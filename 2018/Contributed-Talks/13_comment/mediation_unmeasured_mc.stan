data {
  // number of observations
  int<lower=0> N;
  // number of columns in design matrix excluding A
  int<lower=0> P;
  // design matrix, excluding A, M, U
  matrix[N, P] X;
  // observed treatment
  vector[N] A;
  // observed mediator
  int<lower=0,upper=1> M[N];
  // outcome
  int<lower=0,upper=1> Y[N];
  // mean of regression priors
  vector[P + 3] alpha_m;
  vector[P + 2] beta_m;
  vector[P] gamma_m;
  // variance-covariance of regression priors
  cov_matrix[P + 3] alpha_vcv;
  cov_matrix[P + 2] beta_vcv;
  cov_matrix[P] gamma_vcv;
}

transformed data {
  // make vector of 1/N for (classical) bootstrapping
  vector[N] boot_probs = rep_vector(1.0/N, N);
  
  // make vector version of M
  vector[N] Mv = to_vector(M);
}

parameters {
  // regression coefficients (confounder model)
  vector[P] gamma;
  
  // regression coefficients (outcome model)
  vector[P + 3] alpha;
  
  // regression coefficients (mediator model)
  vector[P + 2] beta;
}

transformed parameters {
  // partial M coefficient parameters
  vector[P] betaZ = head(beta, P);
  real betaU = beta[P + 1];
  real betaA = beta[P + 2];
  
  // partial Y coefficient parameters
  vector[P] alphaZ = head(alpha, P);
  real alphaU = alpha[P + 1];
  real alphaA = alpha[P + 2];
  real alphaM = alpha[P + 3];
}

model {
  // linear predictors
  // U regression
  vector[N] eta_u;
  // M regression, if U = 0
  vector[N] eta_mu0;
  // Y regression, if U = 0
  vector[N] eta_yu0;
  
  // log-likelihood contributions for U = 0 and U = 1 cases
  real ll_0;
  real ll_1; 
  
  // calculate linear predictors for U = 0 case
  // will selectively add on betaU and alphaU as needed
  eta_u = X * gamma;
  eta_mu0 = X * betaZ + A * betaA;
  eta_yu0 = X * alphaZ + A * alphaA + Mv * alphaM;
  
  // informative priors
  alpha ~ multi_normal(alpha_m, alpha_vcv);
  beta  ~ multi_normal(beta_m, beta_vcv);
  gamma ~ multi_normal(gamma_m, gamma_vcv);
  
  // likelihood
  for (n in 1:N) {
    // contribution if U = 0
    ll_0 = log_inv_logit(eta_yu0[n]) * Y[n] + 
           log1m_inv_logit(eta_yu0[n]) * (1 - Y[n]) + 
           log_inv_logit(eta_mu0[n]) * M[n] + 
           log1m_inv_logit(eta_mu0[n]) * (1 - M[n]) + 
           log1m_inv_logit(eta_u[n]);
          
    // contribution if U = 1
    ll_1 = log_inv_logit(eta_yu0[n] + alphaU) * Y[n] + 
           log1m_inv_logit(eta_yu0[n] + alphaU) * (1 - Y[n]) + 
           log_inv_logit(eta_mu0[n] + betaU) * M[n] + 
           log1m_inv_logit(eta_mu0[n] + betaU) * (1 - M[n]) + 
           log_inv_logit(eta_u[n]);
    
    // contribution is summation over U possibilities
    target += log_sum_exp(ll_0, ll_1);
  }
}

generated quantities {
  // row index to be sampled for bootstrap
  int row_i;
  
  // calculate NDE in the bootstrapped sample
  real NDE = 0;
  vector[N] U;
  vector[N] M_a0;
  vector[N] Y_a1Ma0;
  vector[N] Y_a0Ma0;
  for (n in 1:N) {
    // sample baseline covariates
    row_i = 1;
    //categorical_rng(boot_probs);
    
    // sample U
    U[n] = bernoulli_logit_rng(X[row_i] * gamma);
    
    // sample M_a where a = 0
    M_a0[n] = bernoulli_logit_rng(X[row_i] * betaZ + U[n] * betaU);
    
    // sample Y_(a=0, M=M_0) and Y_(a=1, M=M_0)
    Y_a0Ma0[n] = bernoulli_logit_rng(X[row_i] * alphaZ + M_a0[n] * alphaM + 
                                     U[n] * alphaU);
    Y_a1Ma0[n] = bernoulli_logit_rng(X[row_i] * alphaZ + M_a0[n] * alphaM +
                                     alphaA + U[n] * alphaU);
                                                          
    // add contribution of this observation to the bootstrapped NDE
    NDE = NDE + (Y_a1Ma0[n] - Y_a0Ma0[n])/N;
  }
}

