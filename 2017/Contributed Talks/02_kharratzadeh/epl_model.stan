data {
  int<lower=1> nteams; // number of teams (20)
  int<lower=1> ngames; // number of games 
  int<lower=1> nweeks; // number of weeks 
  int<lower=1> home_week[ngames]; // week number for the home team
  int<lower=1> away_week[ngames]; // week number for the away team
  int<lower=1, upper=nteams> home_team[ngames]; // home team ID (1, ..., 20)
  int<lower=1, upper=nteams> away_team[ngames]; // away team ID (1, ..., 20)
  vector[ngames] score_diff;    // home_goals - away_goals
  row_vector[nteams] prev_perf; // a score between -1 and +1
}
parameters {
  real b_home; // the effect of hosting the game in mean of score_diff dist.
  real b_prev;                         // regression coefficient of prev_perf
  real<lower=0> sigma_a0;              // teams ability variation 
  real<lower=0> tau_a;                 // hyper-param for game-to-game variation
  real<lower=1> nu;                    // t-dist degree of freedom
  real<lower=0> sigma_y;               // score_diff variation
  row_vector<lower=0>[nteams] sigma_a_raw; // game-to-game variation
  matrix[nweeks,nteams] eta_a;         // random component
}
transformed parameters {
  matrix[nweeks, nteams] a;                        // team abilities
  row_vector<lower=0>[nteams] sigma_a; // game-to-game variation
  a[1] = b_prev * prev_perf + sigma_a0 * eta_a[1]; // initial abilities (at week 1)
  sigma_a = tau_a * sigma_a_raw;
  for (w in 2:nweeks) {
    a[w] = a[w-1] + sigma_a .* eta_a[w];           // evolution of abilities
  }
}
model {
  vector[ngames] a_diff;
  // Priors
  nu ~ gamma(2,0.1);     
  b_prev ~ normal(0,1);
  sigma_a0 ~ normal(0,1);
  sigma_y ~ normal(0,5);
  b_home ~ normal(0,1);
  sigma_a_raw ~ normal(0,1);
  tau_a ~ cauchy(0,1);
  to_vector(eta_a) ~ normal(0,1);
  // Likelihood
  for (g in 1:ngames) {
     a_diff[g] = a[home_week[g],home_team[g]] - a[away_week[g],away_team[g]];
  }
  score_diff ~ student_t(nu, a_diff + b_home, sigma_y);
}
generated quantities {
  vector[ngames] score_diff_rep;
  for (g in 1:ngames)
    score_diff_rep[g] = student_t_rng(nu, a[home_week[g],home_team[g]] - 
      a[away_week[g],away_team[g]]+b_home, sigma_y);
}


