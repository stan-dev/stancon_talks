functions {
  
  real skellam_lpmf(int k, real mu1, real mu2) {
    real total;
    real log_prob;

    total = (- mu1 - mu2) + ((log(mu1) - log(mu2)) * k / 2);

    log_prob = total + log(modified_bessel_first_kind(k, 2 * sqrt(mu1*mu2)));

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
  real<lower=0, upper=1> mixing_proportion;
  real constant_mu;
  real home_advantage;
  real b_prev;
  real<lower=0> sigma_a0;              // teams ability variation 
  real<lower=0> tau_a;                 // hyper-param for game-to-game variation
  
  row_vector<lower=0>[n_teams] sigma_a_raw; // game-to-game variation
  matrix[n_weeks,n_teams] eta_a;         // random component
}
transformed parameters {
  matrix[n_weeks, n_teams] a;                        // team abilities
  row_vector<lower=0>[n_teams] sigma_a; // game-to-game variation
  a[1] = b_prev * prev_perf + sigma_a0 * eta_a[1]; // initial abilities (at week 1)
  sigma_a = tau_a * sigma_a_raw;
  for (w in 2:n_weeks) {
    a[w] = a[w-1] + sigma_a .* eta_a[w];           // evolution of abilities
  }
}
model {
  vector[n_games] home_expected_goals;
  vector[n_games] away_expected_goals;

  // Priors
  mixing_proportion ~ uniform(0, 1);
  home_advantage ~ normal(0, 1);
  constant_mu ~ normal(0, 1);

  b_prev ~ normal(0,2);
  sigma_a0 ~ normal(0,2);
 
  sigma_a_raw ~ normal(0,1);
  tau_a ~ cauchy(0,1);
  to_vector(eta_a) ~ normal(0,1);
  // Likelihood
  for (g in 1:n_games) {
    home_expected_goals[g] = exp(
      constant_mu + home_advantage +
        a[home_week[g], home_team[g]] - 
        a[away_week[g], away_team[g]]
    );
    
    away_expected_goals[g] = exp(
      constant_mu + 
      a[away_week[g], away_team[g]] - 
      a[home_week[g], home_team[g]]
    );
    
    goal_difference[g] ~ zero_inflated_skellam(
      home_expected_goals[g],
      away_expected_goals[g],
      mixing_proportion
    );
  }
}

generated quantities {
  vector[n_games] goal_difference_rep;
  vector[n_games] home_expected_goals;
  vector[n_games] away_expected_goals;
  vector[n_games_pred] goal_difference_pred_rep;
  
  // Likelihood
  for (g in 1:n_games) {
    home_expected_goals[g] = exp(
      constant_mu + home_advantage +
        a[home_week[g],home_team[g]] - 
        a[away_week[g],away_team[g]]
    );
    
    away_expected_goals[g] = exp(
      constant_mu + 
      a[away_week[g],away_team[g]] - 
      a[home_week[g],home_team[g]]
    );
    
    goal_difference_rep[g] = zero_inflated_skellam_rng(
      home_expected_goals[g],
      away_expected_goals[g],
      mixing_proportion
    );
  }
  // out of sample predictions
  // recycle home_expected_goals & away_expected_goals
  for (g in 1:n_games_pred) {
    home_expected_goals[g] = exp(
      constant_mu + home_advantage +
        a[n_weeks,home_team_pred[g]] - 
        a[n_weeks,away_team_pred[g]]
    );
    
    away_expected_goals[g] = exp(
      constant_mu + 
      a[n_weeks,away_team_pred[g]] - 
      a[n_weeks,home_team_pred[g]]
    );
    
    goal_difference_pred_rep[g] = zero_inflated_skellam_rng(
      home_expected_goals[g],
      away_expected_goals[g],
      mixing_proportion
    );
  }
}
