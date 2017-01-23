data {
  int<lower=1> N;
  int<lower=1> N_states;
  int<lower=1> N_regions;
  int<lower=1> N_years_obs;
  int<lower=1> N_years;
  int<lower=1> state_region_ind[N_states];
  int<lower=1,upper=50> state_ind[N];
  int<lower=1,upper=10> region_ind[N];
  int<lower=1> year_ind[N];
  vector<lower=0>[N] turnout_weight;
  vector<lower=0,upper=1>[N] y;
}
transformed data {
  real years[N_years];
  vector[17] counts;

  for (t in 1:N_years)
    years[t] = t;
  for (i in 1:17)
    counts[i] = 2;
}
parameters {
  matrix[N_years,N_regions] GP_region_std;
  matrix[N_years,N_states] GP_state_std;
  vector[N_years_obs] year_std;
  vector[N_states] state_std;
  vector[N_regions] region_std;
  real<lower=0> tot_var;
  simplex[17] prop_var;
  real mu;


  real<lower=0> length_GP_region_long;
  real<lower=0> length_GP_state_long;
  real<lower=0> length_GP_region_short;
  real<lower=0> length_GP_state_short;
}
transformed parameters {
  matrix[N_years,N_regions] GP_region;
  matrix[N_years,N_states] GP_state;

  vector[N_years_obs] year_re;
  vector[N_states] state_re;
  vector[N_regions] region_re;
  vector[17] vars;

  real sigma_year;
  real sigma_region;
  vector[10] sigma_state;

  real sigma_error_state_2;

  real sigma_GP_region_long;
  real sigma_GP_state_long;
  real sigma_GP_region_short;
  real sigma_GP_state_short;

  vars = 17 * prop_var * tot_var;
  sigma_year = sqrt(vars[1]);
  sigma_region = sqrt(vars[2]);
  for (i in 1:10)
    sigma_state[i] = sqrt(vars[i + 2]);

  sigma_GP_region_long = sqrt(vars[13]);
  sigma_GP_state_long = sqrt(vars[14]);
  sigma_GP_region_short = sqrt(vars[15]);
  sigma_GP_state_short = sqrt(vars[16]);
  sigma_error_state_2 = sqrt(vars[17]);

  region_re = sigma_region * region_std;
  year_re = sigma_year * year_std;
  state_re = sigma_state[state_region_ind] .* state_std;
  
  {
    matrix[N_years, N_years] cov_region; 
    matrix[N_years, N_years] cov_state; 
    matrix[N_years, N_years] L_cov_region; 
    matrix[N_years, N_years] L_cov_state; 

    cov_region = cov_exp_quad(years, sigma_GP_region_long, 
                                  length_GP_region_long)
               + cov_exp_quad(years, sigma_GP_region_short, 
                                  length_GP_region_short);
    cov_state = cov_exp_quad(years, sigma_GP_state_long, 
                                  length_GP_state_long)
               + cov_exp_quad(years, sigma_GP_state_short, 
                                  length_GP_state_short);
    for (year in 1:N_years) {
      cov_region[year, year] = cov_region[year, year] + 1e-6;
      cov_state[year, year] = cov_state[year, year] + 1e-6;
    }

    L_cov_region = cholesky_decompose(cov_region);
    L_cov_state = cholesky_decompose(cov_state);
    GP_region = L_cov_region * GP_region_std;
    GP_state = L_cov_state * GP_state_std;
  }
}
model {
  vector[N] obs_mu;

  for (n in 1:N) {
    obs_mu[n] = mu + year_re[year_ind[n]] 
              + state_re[state_ind[n]] 
              + region_re[region_ind[n]]
              + GP_region[year_ind[n],region_ind[n]]
              + GP_state[year_ind[n],state_ind[n]];
  }
  y ~ normal(obs_mu, sigma_error_state_2); #* turnout_weight);

  to_vector(GP_region_std) ~ normal(0, 1);
  to_vector(GP_state_std) ~ normal(0, 1);
  year_std ~ normal(0, 1);
  state_std ~ normal(0, 1);
  region_std ~ normal(0, 1);
  mu ~ normal(.5, .5);
  tot_var ~ gamma(3, 3);
  prop_var ~ dirichlet(counts);
  length_GP_region_long ~ weibull(30,8);
  length_GP_state_long ~ weibull(30,8);
  length_GP_region_short ~ weibull(30,3);
  length_GP_state_short ~ weibull(30,3);
}
generated quantities {
  matrix[N_years,N_states] y_new;
  matrix[N_years,N_states] y_new_pred;

  {
    real level;
    level = normal_rng(0.5, sigma_year);
    for (state in 1:N_states) {
      for (t in 1:N_years) {
        if (t < 12) {
          y_new[t,state] = state_re[state] 
                         + region_re[state_region_ind[state]]
                         + GP_state[t,state]
                         + GP_region[t,state_region_ind[state]]
                         + (mu + year_re[t]);
        } else {
          y_new[t,state] = state_re[state] 
                         + region_re[state_region_ind[state]]
                         + GP_state[t,state]
                         + GP_region[t,state_region_ind[state]]
                         + level;
        }
        y_new_pred[t,state] = normal_rng(y_new[t,state],
                                         sigma_error_state_2);
      }
    }
  }
}
