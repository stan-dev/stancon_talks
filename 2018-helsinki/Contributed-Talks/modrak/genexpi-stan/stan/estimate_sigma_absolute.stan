data { 
  int num_measurements; 
  int num_targets; 
  vector<lower=0>[num_measurements] expression[num_targets]; 
  real<lower=0> sigma_absolute_prior_sigma; 
  real<lower=0> sigma_relative; 
} 
 
parameters { 
  real<lower=0> sigma_absolute; 
} 
 
model { 
  for(i in 1:num_targets) { 
    expression[i] ~ normal(0, sigma_absolute + sigma_relative * expression[i]); 
  } 
  sigma_absolute_prior_sigma ~ normal(0, sigma_absolute_prior_sigma); 
}