data {
  int<lower=1> n_teams; // number of teams (20)
  int<lower=1> n_games; // number of games 
  int<lower=1> n_weeks; // number of weeks 
  int<lower=1> home_week[n_games]; // week number for the home team
  int<lower=1> away_week[n_games]; // week number for the away team
  int<lower=1, upper=n_teams> home_team[n_games]; // home team ID (1, ..., 20)
  int<lower=1, upper=n_teams> away_team[n_games]; // away team ID (1, ..., 20)
  vector[n_games] goal_difference;    // home_goals - away_goals
  row_vector[n_teams] prev_perf; // a score between -1 and +1
  // out-of-sample predictions for generated quanties block
  int<lower=1> n_games_pred; // number of games
  int<lower=1, upper=n_teams> home_team_pred[n_games_pred];
  int<lower=1, upper=n_teams> away_team_pred[n_games_pred];
}
parameters {
  real b_home; // the effect of hosting the game in mean of score_diff dist.
  //real b_prev;                         // regression coefficient of prev_perf
  real<lower=0> sigma_a0;              // teams ability variation 
  real<lower=0> tau_a;                 // hyper-param for game-to-game variation
  real<lower=1> nu;                    // t-dist degree of freedom
  real<lower=0> sigma_y;               // score_diff variation
  row_vector<lower=0>[n_teams] sigma_a_raw; // game-to-game variation
  matrix[n_weeks,n_teams] eta_a;         // random component
}
transformed parameters {
  matrix[n_weeks, n_teams] a;                        // team abilities
  row_vector<lower=0>[n_teams] sigma_a; // game-to-game variation
  a[1] = sigma_a0 * eta_a[1]; // initial abilities (at week 1)
  sigma_a = tau_a * sigma_a_raw;
  for (w in 2:n_weeks) {
    a[w] = a[w-1] + sigma_a .* eta_a[w];           // evolution of abilities
  }
}
model {
  vector[n_games] a_diff;
  // Priors
  nu ~ gamma(2,0.1);     
  //b_prev ~ normal(0,1);
  sigma_a0 ~ normal(0,1);
  sigma_y ~ normal(0,5);
  b_home ~ normal(0,1);
  sigma_a_raw ~ normal(0,1);
  tau_a ~ cauchy(0,1);
  to_vector(eta_a) ~ normal(0,1);
  // Likelihood
  for (g in 1:n_games) {
    a_diff[g] = a[home_week[g],home_team[g]] - a[away_week[g],away_team[g]];
  }
  goal_difference ~ student_t(nu, a_diff + b_home, sigma_y);
}
generated quantities {
  vector[n_games] goal_difference_rep;
  vector[n_games_pred] goal_difference_pred_rep;
  
  for (g in 1:n_games)
    goal_difference_rep[g] = student_t_rng(nu, a[home_week[g],home_team[g]] - 
                                        a[away_week[g],away_team[g]]+b_home, sigma_y);
  // predictions for next week's round of games, using final week abilities
  for (g in 1:n_games_pred)
    goal_difference_pred_rep[g] = student_t_rng(nu, a[n_weeks, home_team_pred[g]] - 
                                        a[n_weeks, away_team_pred[g]]+b_home, sigma_y);                                      
}


