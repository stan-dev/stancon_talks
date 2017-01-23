functions {
  real da(int winner, real RT, vector beta, real P_b, real mu_da, 
               real mu_b, real sigma, real psi){
    // theta = softmax(beta)
    // log(P(w = 1 | theta, P_b)): 
    real log_P_w1;
    // Prob of direct access given winner = 1 
    real log_P_da_gw1;
    // Prob of backtracking given winner = 1 
    real log_P_b_gw1;
    // Equation (10) in log:
    log_P_w1 = log_sum_exp(categorical_logit_lpmf(1 | beta),
                log(P_b)+ log1m_exp(categorical_logit_lpmf(1|beta)));
    // Equation (14) in log:
    log_P_da_gw1 = categorical_logit_lpmf(1 | beta) - log_P_w1;
    // Equation (15) in log:
    log_P_b_gw1 = log(P_b) + log1m_exp(categorical_logit_lpmf(1 | beta)) -
                  log_P_w1;
    if(winner==1) {
    return (log_P_w1 + // Increment on likelihood due to winner=1
                       // Increment on likelihood due to RT:
            log_sum_exp(log_P_da_gw1 + lognormal_lpdf(RT - psi| mu_da, sigma),
            log_P_b_gw1 + lognormal_lpdf(RT - psi | mu_da + mu_b, sigma) ));
    } else {
      return (log1m(P_b) + categorical_logit_lpmf(winner | beta) + 
              // Increment on likelihood due to RT:
              lognormal_lpdf(RT - psi | mu_da, sigma));
    }
  }
}
data {
  int<lower = 0> N_obs; 
  int<lower = 1> N_choices; 
  int<lower = 1, upper = N_choices> winner[N_obs];
  vector<lower = 0>[N_obs] RT;
}
transformed data {
  real<lower=0> min_RT;
  min_RT = min(RT);
}
parameters{
  real<lower=0,upper=1> P_b;
  real<lower=0> mu_da;
  real<lower=0> mu_b;
  vector[N_choices-2] beta_incorrect; 
  real<lower=0> beta_added;
  real<lower=0> sigma;
  real<lower=0,upper=min_RT> psi;
}
transformed parameters{
  vector[N_choices] beta;
  beta[1] = beta_added + fmax(max(beta_incorrect),0);
  beta[2:N_choices-1] = beta_incorrect;
  beta[N_choices] = 0;
}
model {
  beta_added ~ normal(0,2);
  beta_incorrect ~ normal(0,2);
  P_b ~ beta(1,1);
  mu_da ~ normal(0,10);
  mu_b ~ normal(0,2);
  sigma ~ normal(0,2);
  psi ~ normal(0,300);
  for (n in 1:N_obs) {
    target +=  da(winner[n], RT[n], beta, P_b, mu_da, mu_b, sigma, psi);
  }
}
generated quantities {
  vector[N_choices] theta;
  theta = softmax(beta);
}
