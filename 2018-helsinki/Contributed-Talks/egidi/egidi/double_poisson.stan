data {
  int N;                                  // number of matches
  int nseason;                            // number of seasons
  int season[N];                          // season index
  int score_home[N];                      // home scores
  int score_away[N];                      // away scores
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
  matrix[nseason, nteams] att_raw;        // raw attack ability
  matrix[nseason, nteams] def_raw;        // raw defense ability
  }
 transformed parameters{
  matrix[nseason, nteams] att;            // attack abilities
  matrix[nseason, nteams] def;            // defense abilities
  cov_matrix[nseason] Sigma_att;          // GP attack cov. funct.
  cov_matrix[nseason] Sigma_def;          // GP defense cov.funct.
  matrix[nseason, nteams] mu_att;         // attack hyperparameter
  matrix[nseason, nteams] mu_def;         // defense hyperparameter
  vector<lower=0>[N] theta_home;          // scoring h. rate
  vector<lower=0>[N] theta_away;          // scoring a. rate
  
  theta_home = rep_vector(0.2, N);
  theta_away = rep_vector(0.2, N);
  
  
   
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
     mu_def[1]=rep_row_vector(0,nteams);
     mu_def[t]=rep_row_vector(0,nteams);
     }
     
     // scoring rates
  for (n in 1:N){
  theta_home[n]=exp(mu+att[season[n], teams_home[n]]                                     +def[season[n], teams_away[n]]);
  theta_away[n]=exp(att[season[n], teams_away[n]]
                      +def[season[n], teams_home[n]]);
     }

}
 model{
  // Priors
 
   
   for (h in 1:(nteams)){
     att_raw[,h]~multi_normal(mu_att[,h], Sigma_att);
     def_raw[,h]~multi_normal(mu_def[,h], Sigma_def);
   }
  
  // Likelihood
   for (n in 1:N){
     target+= poisson_lpmf(score_home[n] |  theta_home[n]);
     target+= poisson_lpmf(score_away[n] |  theta_away[n]);
     }
}
 generated quantities{
  vector[N] log_lik;
  int score_home_rep[N];        // in-sample replication for home scores
  int score_away_rep[N];        // in-sample replications for away scores 
  int score_home_prev[N_prev];  // out-of-sample replications for home scores 
  int score_away_prev[N_prev];  // out-of-sample replications for away scores 
  vector<lower=0>[N_prev] theta_home_prev;   // predicted home shooting rate 
  vector<lower=0>[N_prev] theta_away_prev;   // predicted away shooting rate 
  
    
    for (n in 1:N){
      score_home_rep[n]=poisson_rng(theta_home[n]);
      score_away_rep[n]=poisson_rng(theta_away[n]);
      log_lik[n] = poisson_lpmf(score_home[n] |  theta_home[n])+
                   poisson_lpmf(score_away[n] |  theta_away[n]);
      }
 
  // p_prev initialization
   theta_home_prev=theta_home[1:N_prev];
   theta_away_prev= theta_away[1:N_prev];
 
    for (n in 1:N_prev){
      theta_home_prev[n]=exp(mu+att[season_prev[n], teams_home_prev[n]]                                     +def[season_prev[n], teams_away_prev[n]]);
      theta_away_prev[n]=exp(att[season_prev[n], teams_away_prev[n]]
                              +def[season_prev[n], teams_home_prev[n]]);
      score_home_prev[n]=poisson_rng(theta_home_prev[n]);
      score_away_prev[n]=poisson_rng(theta_away_prev[n]);
        }
}
