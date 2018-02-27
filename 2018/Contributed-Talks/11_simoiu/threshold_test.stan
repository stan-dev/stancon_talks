data {
  int<lower=1> N; // number of observations
  int<lower=1> R; // number of suspect races
  int<lower=1> D; // number of counties
  
  int<lower=1,upper=R> r[N]; // race of suspect
  int<lower=1,upper=D> d[N]; // county where stop occurred
  
  int<lower=1> n[N]; // # of stops
    int<lower=0> s[N]; // # of searches
    int<lower=0> h[N]; // # of successful searches (hits)
}

parameters {	
  // hyperparameters
  real<lower=0> sigma_t; #standard deviation for the normal the thresholds are drawn from. 
  
  // search thresholds
  vector[R] t_r;
  vector[N] t_i_raw;
  
  // parameters for signal distribution
  vector[R] phi_r;
  vector[D-1] phi_d_raw;
  real mu_phi;
  
  vector[R] delta_r; 
  vector[D-1] delta_d_raw;
  real mu_delta;
}

transformed parameters {
  vector[D] phi_d;
  vector[D] delta_d;
  vector[N] phi;
  vector[N] delta;
  vector[N] t_i;
  vector<lower=0, upper=1>[N] search_rate;
  vector<lower=0, upper=1>[N] hit_rate;
  real successful_search_rate;
  real unsuccessful_search_rate;

  phi_d[1]      = 0;
  phi_d[2:D]    = phi_d_raw;
  delta_d[1]   = 0;
  delta_d[2:D] = delta_d_raw;
  
  t_i = t_r[r] + t_i_raw * sigma_t;
  
  for (i in 1:N) {	
    
    // phi is the fraction of people of race r, d who are guilty (ie, carrying contraband)
    phi[i]    = inv_logit(phi_r[r[i]] + phi_d[d[i]]);
    
    // mu is the center of the guilty distribution. 
    delta[i] = exp(delta_r[r[i]] + delta_d[d[i]]);
    
    successful_search_rate = phi[i] * (1 - normal_cdf(t_i[i], delta[i], 1));
    unsuccessful_search_rate = (1 - phi[i]) * (1 - normal_cdf(t_i[i], 0, 1)); 
    search_rate[i] = (successful_search_rate + unsuccessful_search_rate);
    hit_rate[i] = successful_search_rate / search_rate[i];
   }
}

model {  
  // Draw threshold hyperparameters
  sigma_t ~ normal(0, 1);
  
  // Draw race parameters. Each is centered at a mu, and we allow for inter-race heterogeneity. 
  mu_phi ~ normal(0, 1);
  mu_delta ~ normal(0, 1);
  
  phi_r    ~ normal(mu_phi, .1);
  delta_r ~ normal(mu_delta, .1);
  t_r ~ normal(0, 1);
  
  // Draw department parameters (for un-pinned departments)
  phi_d_raw    ~ normal(0, .1);   
  delta_d_raw ~ normal(0, .1);    
  //thresholds
  t_i_raw ~ normal(0, 1);

  s ~ binomial(n, search_rate);
  h ~ binomial(s, hit_rate);
}
