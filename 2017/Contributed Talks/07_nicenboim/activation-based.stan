functions {
    real race(int winner, real RT, real[] alpha, real b, real sigma, real psi){
    real log_lik;
    int N_choices;
    N_choices = num_elements(alpha);
    log_lik = 0;
    for(c in 1:N_choices)
        if(c == winner)
          log_lik = log_lik + lognormal_lpdf(RT - psi|b - alpha[c], sigma);
        else 
          log_lik = log_lik + lognormal_lccdf(RT - psi|b - alpha[c], sigma); 
    return(log_lik);
    }
}
data { 
  int<lower = 0> N_obs; 
  int<lower = 1> N_choices; 
  int<lower = 1, upper = N_choices> winner[N_obs];
  vector<lower = 0>[N_obs] RT;
}
transformed data {
  real  b; //arbitrary threshold
  real min_RT;
  b = 10;
  min_RT = min(RT);
}
parameters{
  real alpha[N_choices]; 
  real<lower=0> sigma;
  real<lower=0,upper=min_RT> psi;
}
model {
  alpha ~ normal(0,10);
  sigma ~ normal(0,2);
  psi ~ normal(0,300); 
  for (n in 1:N_obs) {
     target += race(winner[n], RT[n], alpha, b, sigma, psi);
  }
}
