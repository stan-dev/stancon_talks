data {
  int N;                                  // number of matches
  int nseason;                            // number of seasons
  int season[N];                          // season index
  int score_home[N];                      // home scores
  int score_away[N];                      // away scores
  int Shots_home[N];                      // home shots on target
  int Shots_away[N];                      // away total shots on target
  int Tot_Shots_home[N];                  // home total shots  
  int Tot_Shots_away[N];                  // away total shots
  int nteams;                             // number of teams  
  int teams_home[N];                      // team home index  
  int teams_away[N];                      // team away index  
  int N_prev;                             // number of predicted matches  
  int teams_home_prev[N_prev];            // team home predicted index
  int teams_away_prev[N_prev];            // team away predicted index
  int season_prev[N_prev];                // season predicted index
  int season_time[nseason];               
  }
 parameters {
  real mu;                                // home effect parameter
  real<lower=0> phi;                      // overdispersion total shots parameter
  vector<lower=0>[nteams] alpha;          // total shots shape
  vector<lower=0>[nteams] beta;           // total shots rate
  vector<lower=0>[nteams] delta;          // beta shape 1
  vector<lower=0>[nteams] epsilon;        // beta shape 2
  matrix[nseason, nteams] att_raw;        // raw attack ability
  matrix[nseason, nteams] def_raw;        // raw defense ability
  vector<lower=0>[N] theta_home;          // home total shots rate
  vector<lower=0>[N] theta_away;          // away total shots rate
}
 transformed parameters{
  matrix[nseason, nteams] att;            // attack abilities
  matrix[nseason, nteams] def;            // defense abilities
  cov_matrix[nseason] Sigma_att;          // Gaussian process attack cov. funct.
  cov_matrix[nseason] Sigma_def;          // Gaussian process defense cov.funct.
  matrix[nseason, nteams] mu_att;         // attack hyperparameter
  matrix[nseason, nteams] mu_def;         // defense hyperparameter
  matrix<lower=0, upper=1>[nteams,nseason-1] p;   // conversion home probabilities
  matrix<lower=0, upper=1>[nteams,nseason-1] q;   // conversion away probabilities
  
  p = rep_matrix(0.2, nteams, nseason-1);
  q = rep_matrix(0.2, nteams, nseason-1);
  
   
  // Gaussian process covariance functions
   for (i in 1:(nseason)){
     for (j in 1:(nseason)){
       Sigma_att[i, j] = exp(-pow(season_time[i] - season_time[j], 2))
       + (i == j ? 0.1 : 0.0);
       Sigma_def[i, j] = exp(-pow(season_time[i] - season_time[j], 2))
                   + (i == j ? 0.1 : 0.0);
     }}
   
  // Sum-to-zero constraint for attack/defense parameters
  att[1]=att_raw[1]-mean(att_raw[1]);
  def[1]=def_raw[1]-mean(def_raw[1]);
   for (t in 2:nseason){
      att[t]=att_raw[t]-mean(att_raw[t]);
      def[t]=def_raw[t]-mean(def_raw[t]);
     }
 
  // Lagged prior mean for attack/defense parameters 
   for (t in 2:(nseason)){
     mu_att[1]=rep_row_vector(0,nteams);
     mu_att[t]=rep_row_vector(0,nteams);
     //att[t-1];
     mu_def[1]=rep_row_vector(0,nteams);
     mu_def[t]=rep_row_vector(0,nteams);
     //def[t-1];
     }
 
  // Conversion probabilities
   for (n in 1:N){
     p[teams_home[n], season[n]]=inv_logit(mu+att[season[n], teams_home[n]]  
                                    +def[season[n], teams_away[n]]);
     q[teams_away[n], season[n]]=inv_logit(att[season[n], teams_away[n]]
                                    +def[season[n], teams_home[n]]);
     }
}
 model{
  // Priors
  target+=normal_lpdf(phi|0,1);
  target+=normal_lpdf(alpha |0,5);
  target+=normal_lpdf(beta |0,5);
  target+=normal_lpdf(delta |0,5);
  target+=normal_lpdf(epsilon |0,5);
   
   for (h in 1:(nteams)){
     att_raw[,h]~multi_normal(mu_att[,h], Sigma_att);
     def_raw[,h]~multi_normal(mu_def[,h], Sigma_def);
   }
  
  // Likelihood
   for (n in 1:N){
     target+= binomial_lpmf(score_home[n] |  Shots_home[n],  
                            p[teams_home[n],season[n]]);
     target+= binomial_lpmf(score_away[n] |  Shots_away[n], 
                            q[teams_away[n], season[n]]);
     target+= beta_binomial_lpmf(Shots_home[n] |  Tot_Shots_home[n],     
                                 delta[teams_home[n]],epsilon[teams_home[n]]);
     target+= beta_binomial_lpmf(Shots_away[n] |  Tot_Shots_away[n], 
                                 delta[teams_away[n]],epsilon[teams_away[n]]);
     target+=gamma_lpdf(theta_home[n]|alpha[teams_home[n]],beta[teams_home[n]]);
     target+=gamma_lpdf(theta_away[n]|alpha[teams_away[n]],beta[teams_away[n]]);
      }
 
  target+= neg_binomial_2_lpmf(Tot_Shots_home|theta_home, phi );
  target+= neg_binomial_2_lpmf(Tot_Shots_away|theta_away, phi);
}
 generated quantities{
  vector[N] log_lik;
  int score_home_rep[N];        // in-sample replication for home scores
  int score_away_rep[N];        // in-sample replications for away scores 
  int Shots_home_rep[N];        // in-sample replications for home shots on t.
  int Shots_away_rep[N];        // in-sample replications for away shots on t.
  int Tot_Shots_home_rep[N];    // in-sample replications for home total shots
  int Tot_Shots_away_rep[N];    // in-sample replications for away total shots
  int score_home_prev[N_prev];  // out-of-sample replications for home scores 
  int score_away_prev[N_prev];  // out-of-sample replications for away scores 
  int Shots_home_prev[N_prev];  // out-of-sample replications for home shots on t. 
  int Shots_away_prev[N_prev];  // out-of-sample replications for away shots on t.
  int Tot_Shots_home_prev[N_prev];  // out-of-sample repl. for home total shots
  int Tot_Shots_away_prev[N_prev];  // out-of-sample repl. for away total shots
  vector<lower=0>[N_prev] theta_home_prev;   // predicted home shooting rate 
  vector<lower=0>[N_prev] theta_away_prev;   // predicted away shooting rate 
  vector<lower=0, upper=1>[nteams] p_prev;   // pred. conversion home probabilities
  vector<lower=0, upper=1>[nteams] q_prev;   // pred. conversion away probabilities
    
    for (n in 1:N){
      Tot_Shots_home_rep[n]=neg_binomial_2_rng(theta_home[n], phi);
      Tot_Shots_away_rep[n]=neg_binomial_2_rng(theta_away[n], phi);
      Shots_home_rep[n]=beta_binomial_rng(Tot_Shots_home_rep[n],
                        delta[teams_home[n]],epsilon[teams_home[n]]);
      Shots_away_rep[n]=beta_binomial_rng(Tot_Shots_away_rep[n],
                        delta[teams_away[n]],epsilon[teams_away[n]]);
      score_home_rep[n]=binomial_rng(Shots_home_rep[n],
                                      p[teams_home[n], season[n]]);
      score_away_rep[n]=binomial_rng(Shots_away_rep[n],
                                      q[teams_away[n], season[n]]);
      log_lik[n] =binomial_lpmf(score_home[n] |  Shots_home[n],  
                                p[teams_home[n],season[n]])+
                  binomial_lpmf(score_away[n] |  Shots_away[n], 
                                q[teams_away[n], season[n]]);
                  //               +
                  // beta_binomial_lpmf(Shots_home[n] |
                  //                   Tot_Shots_home[n], delta[teams_home[n]],
                  //                   epsilon[teams_home[n]])+
                  // beta_binomial_lpmf(Shots_away[n] |
                  //                    Tot_Shots_away[n], delta[teams_away[n]],
                  //                    epsilon[teams_away[n]])+
                  // neg_binomial_2_lpmf(Tot_Shots_home[n]| theta_home[n], phi )+
                  // neg_binomial_2_lpmf(Tot_Shots_away[n]| theta_away[n], phi);
      }
 
  // p_prev initialization
   p_prev=p[,nseason-1];
   q_prev=q[,nseason-1];
 
    for (n in 1:N_prev){
      p_prev[teams_home_prev[n]]=p[teams_home_prev[n], nseason-1];
      //inv_logit(mu+att[nseason,teams_home_prev[n]]+
      //  def[nseason, teams_away_prev[n]]);
      q_prev[teams_away_prev[n]]=q[teams_away_prev[n], nseason-1];
      //inv_logit(att[nseason,teams_away_prev[n]]+
      //def[nseason, teams_home_prev[n]]);
      theta_home_prev[n]=gamma_rng(alpha[teams_home_prev[n]],
                                    beta[teams_home_prev[n]]);
      theta_away_prev[n]=gamma_rng(alpha[teams_away_prev[n]],
                                    beta[teams_away_prev[n]]);
      Tot_Shots_home_prev[n]=neg_binomial_2_rng(theta_home_prev[n], phi);
      Tot_Shots_away_prev[n]=neg_binomial_2_rng(theta_away_prev[n], phi);
      Shots_home_prev[n]=beta_binomial_rng(Tot_Shots_home_prev[n],
                  delta[teams_home_prev[n]],epsilon[teams_home_prev[n]]);
      Shots_away_prev[n]=beta_binomial_rng(Tot_Shots_away_prev[n],                                   delta[teams_away_prev[n]],epsilon[teams_away_prev[n]]);
      score_home_prev[n]=binomial_rng(Shots_home_prev[n],
                                        p_prev[teams_home_prev[n]]);
      score_away_prev[n]=binomial_rng(Shots_away_prev[n],
                                        q_prev[teams_away_prev[n]]);
        }
}
