functions {
  
  real skellam_lpmf(int k, real mu1, real mu2) {
    real total;
    real log_prob;
    real tmp;

    total = (- mu1 - mu2) + ((log(mu1) - log(mu2)) * k / 2);
    
    tmp = 2 * sqrt(mu1*mu2);
    if(tmp > 700) log_prob = total;
    else
      log_prob = total + log(modified_bessel_first_kind(k, tmp));

    return log_prob;
  }
  
  real zero_inflated_skellam_lpmf(int k, real mu1, real mu2, real p) {
    real base_prob;
    real prob;
    real log_prob;
    
    base_prob = exp(skellam_lpmf(k | mu1, mu2));
    
    if (k == 0) // i.e. a draw
      prob = p + ((1 - p) * base_prob);
    else
      prob = (1 - p) * base_prob;
    
    log_prob = log(prob);
    
    return log_prob;
  }
  real zero_inflated_skellam_rng(real mu1, real mu2, real p){
    int k1_rng;
    int k2_rng;
    real u_rng;
    
    u_rng = uniform_rng(0, 1);
    
    if(u_rng < p)
      return 0;
    else {
      k1_rng = poisson_rng(mu1);
      k2_rng = poisson_rng(mu2);
      return k1_rng - k2_rng;
    }
  }
  real skellam_rng(real mu1, real mu2){
      return poisson_rng(mu1) - poisson_rng(mu2);
  }  
}

data {
  int<lower=1> n_teams;
  int<lower=1> n_games;
  int<lower=1, upper=n_teams> home_team[n_games];
  int<lower=1, upper=n_teams> away_team[n_games];
  int goal_difference[n_games];
  
  int<lower=1> n_weeks; // number of weeks 
  int<lower=1> home_week[n_games]; // week number for the home team
  int<lower=1> away_week[n_games]; // week number for the away team
  
  row_vector[n_teams] prev_perf; // a score between -1 and +1
  
  // out-of-sample predictions for generated quanties block
  int<lower=1> n_games_pred; // number of games
  int<lower=1, upper=n_teams> home_team_pred[n_games_pred];
  int<lower=1, upper=n_teams> away_team_pred[n_games_pred];
}

parameters {
  //real<lower=0, upper=1> mixing_proportion;
  real constant_mu;
  real home_advantage;
  real b_prev_offense;
  real b_prev_defense;// regression coefficient of prev_perf
  real<lower=0> sigma_a0;              // teams ability variation 
  real<lower=0> tau_a;                 // hyper-param for game-to-game variation
  
  row_vector<lower=0>[n_teams] sigma_a_raw; // game-to-game variation
  matrix[n_weeks,n_teams] eta_a_offense;         // random component
  matrix[n_weeks,n_teams] eta_a_defense;         // random component
}

transformed parameters {
  // team abilities
  matrix[n_weeks, n_teams] a_offense;   
  matrix[n_weeks, n_teams] a_defense;   

  row_vector<lower=0>[n_teams] sigma_a; // game-to-game variation
  
  // vectorized over teams
  a_offense[1] = b_prev_offense * prev_perf + sigma_a0 * eta_a_offense[1]; 
  // initial abilities (at week 1)
  a_defense[1] = b_prev_defense * prev_perf + sigma_a0 * eta_a_defense[1]; 
  // initial abilities (at week 1)

  sigma_a = tau_a * sigma_a_raw;
  // evolution of abilities
  for (w in 2:n_weeks) {
    a_offense[w] = a_offense[w-1] + sigma_a .* eta_a_offense[w]; 
    a_defense[w] = a_defense[w-1] + sigma_a .* eta_a_defense[w]; 
  }
}

model {
  vector[n_games] home_expected_goals;
  vector[n_games] away_expected_goals;

  // Priors
  //mixing_proportion ~ uniform(0, 1);
  target += normal_lpdf(home_advantage | 0, 1);
  target += normal_lpdf(constant_mu | 0, 1);

  target += normal_lpdf(b_prev_offense | 0, 2);
  target += normal_lpdf(b_prev_defense | 0, 2);
  target += normal_lpdf(sigma_a0 | 0, 2);
 
  target += normal_lpdf(sigma_a_raw | 0, 1);
  target += cauchy_lpdf(tau_a | 0, 1);
  for (w in 1:n_weeks)
    target += normal_lpdf(eta_a_offense[w] | 0, 1);
  for (w in 1:n_weeks)
    target += normal_lpdf(eta_a_defense[w] | 0, 1);

  // Likelihood
  for (g in 1:n_games) {
    home_expected_goals[g] = exp(
      constant_mu + home_advantage +
        a_offense[home_week[g], home_team[g]] + 
        a_defense[away_week[g], away_team[g]]
    );
    
    away_expected_goals[g] = exp(
      constant_mu + 
      a_offense[away_week[g], away_team[g]] + 
      a_defense[home_week[g], home_team[g]]
    );
    
    target += skellam_lpmf(goal_difference[g] | home_expected_goals[g],
      away_expected_goals[g]);
  }
}
