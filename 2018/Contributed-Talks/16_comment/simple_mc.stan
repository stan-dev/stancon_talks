data {
  // number of observations
  int<lower=0> N;
  // number of columns in design matrix excluding A
  int<lower=0> P;
  // design matrix, excluding treatment A
  matrix[N, P] X;
  // observed treatment
  vector[N] A;
  // outcome
  int<lower=0,upper=1> Y[N];
}

transformed data {
  // make vector of 1/N for (classical) bootstrapping
  vector[N] boot_probs = rep_vector(1.0/N, N);
}

parameters {
  // regression coefficients
  vector[P + 1] alpha;
}

transformed parameters {
  vector[P] alphaZ = head(alpha, P);
  real alphaA = alpha[P + 1];
}

model {
  // priors for regression coefficients
  alpha ~ normal(0, 2.5);
  
  // likelihood
  Y ~ bernoulli_logit(X * alphaZ + A * alphaA);
}

generated quantities {
  // row index to be sampled for bootstrap
  int row_i;

  // calculate ATE in the bootstrapped sample
  real ATE = 0;
  vector[N] Y_a1;
  vector[N] Y_a0;
  for (n in 1:N) {
    // sample baseline covariates
    row_i = categorical_rng(boot_probs);
    
    // sample Ya where a = 1 and a = 0
    Y_a1[n] = bernoulli_logit_rng(X[row_i] * alphaZ + alphaA);
    Y_a0[n] = bernoulli_logit_rng(X[row_i] * alphaZ);

    // add contribution of this observation to the bootstrapped ATE
    ATE = ATE + (Y_a1[n] - Y_a0[n])/N;
  }
}

