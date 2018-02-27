data {
  // number of observations
  int<lower=0> N;
  // number of columns in design matrix excluding A (and M)
  int<lower=0> P;
  // design matrix, excluding treatment A
  matrix[N, P] X;
  // observed treatment
  vector[N] A;
  // observed mediator
  int<lower=0,upper=1> M[N];
  // outcome
  int<lower=0,upper=1> Y[N];
  // mean of regression priors
  vector[P + 2] alpha_m;
  vector[P + 1] beta_m;
  // variance-covariance of regression priors
  cov_matrix[P + 2] alpha_vcv;
  cov_matrix[P + 1] beta_vcv;
}

transformed data {
  // make vector of 1/N for (classical) bootstrapping
  vector[N] boot_probs = rep_vector(1.0/N, N);

  // make vector version of M
  vector[N] Mv = to_vector(M);
}

parameters {
  // regression coefficients (outcome model)
  vector[P + 2] alpha;

  // regression coefficients (mediator model)
  vector[P + 1] beta;
}

transformed parameters {
  // partial M coefficient parameters
  vector[P] betaZ = head(beta, P);
  real betaA = beta[P + 1];
  
  // partial Y coefficient parameters
  vector[P] alphaZ = head(alpha, P);
  real alphaA = alpha[P + 1];
  real alphaM = alpha[P + 2];
}

model {
  // priors on causal coefficients weakly informative for binary exposure
  alpha ~ multi_normal(alpha_m, alpha_vcv);
  beta ~ multi_normal(beta_m, beta_vcv);

  // likelihoods
  M ~ bernoulli_logit(X * betaZ + A * betaA);
  Y ~ bernoulli_logit(X * alphaZ + A * alphaA + Mv * alphaM);
}

generated quantities {
  // row index to be sampled for bootstrap
  int row_i;

  // calculate NDE in the bootstrapped sample
  real NDE = 0;
  vector[N] M_a0;
  vector[N] Y_a1Ma0;
  vector[N] Y_a0Ma0;
  for (n in 1:N) {
    // sample baseline covariates
    row_i = categorical_rng(boot_probs);
    
    // sample Ma where a = 0
    M_a0[n] = bernoulli_logit_rng(X[row_i] * betaZ);

    // sample Y_(a=1, M=M_0) and Y_(a=0, M=M_0)
    Y_a1Ma0[n] = bernoulli_logit_rng(X[row_i] * alphaZ + M_a0[n] * alphaM + 
                                     alphaA);
    Y_a0Ma0[n] = bernoulli_logit_rng(X[row_i] * alphaZ + M_a0[n] * alphaM);

    // add contribution of this observation to the bootstrapped NDE
    NDE = NDE + (Y_a1Ma0[n] - Y_a0Ma0[n])/N;
  }
}

