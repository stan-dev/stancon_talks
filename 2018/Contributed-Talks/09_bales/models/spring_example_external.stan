functions {
  vector eigenvalues_sym_external(matrix K);
}
data {
  int<lower = 1> N;
  real<lower = 0.0> y[N - 1];
  real m;
}
transformed data {
  matrix[N, N] K_unscaled = rep_matrix(0, N, N);
  
  for(n in 1:N) {
    if(n == 1) {
      K_unscaled[n, n] = 1.0 / m;
      K_unscaled[n, n + 1] = -1.0 / m;
    } else if(n == N) {
      K_unscaled[n, n - 1] = -1.0 / m;
      K_unscaled[n, n] = 1.0 / m;
    } else {
      K_unscaled[n, n - 1] = -1.0 / m;
      K_unscaled[n, n] = 2.0 / m;
      K_unscaled[n, n + 1] = -1.0 / m;
    }
  }
}
parameters {
  real<lower = 0.0> k;
  real<lower = 0.0> sigma;
}
transformed parameters {
  vector[N - 1] eigs;
  
  {
    matrix[N, N] K = k * K_unscaled;
    
    eigs = eigenvalues_sym_external(K)[2:N];
  }
}
model {
  k ~ normal(1.0, 1.0);
  sigma ~ normal(0.0, 1.0);
  
  y ~ normal(sqrt(eigs), sigma);
}
