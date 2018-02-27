data{
 int<lower=0> N;
 int<lower=0> K;
 int<lower=0> Kb;
 int<lower=0> Kc;
 int<lower=0, upper=1> yb[N,Kb];
 vector[Kc] yc[N];
}

transformed data {
  matrix[Kc, Kc] I = diag_matrix(rep_vector(1, Kc));
}
parameters {
  vector[Kb] zb[N];
  cholesky_factor_corr[K] L_R;  // first continuous, then binary
  vector<lower=0>[Kc] sigma;
  vector[K] mu;
}

transformed parameters{
  matrix[N, Kb] z;
  vector[Kc] mu_c = head(mu, Kc);
  vector[Kb] mu_b = tail(mu, Kb);
  {
    matrix[Kc, Kc] L_inv = mdivide_left_tri_low(diag_pre_multiply(sigma, L_R[1:Kc, 1:Kc]), I);
    for (n in 1:N){
      vector[Kc] resid = L_inv * (yc[n] - mu_c);
      z[n,] = transpose(mu_b + tail(L_R * append_row(resid, zb[n]), Kb));
    }
  }
}

model{
  mu ~ normal(0,10);
  L_R ~ lkj_corr_cholesky(2);
  sigma~cauchy(0,2.5);
  yc ~ multi_normal_cholesky(mu_c, diag_pre_multiply(sigma, L_R[1:Kc, 1:Kc]));
  for (n in 1:N) zb[n] ~ normal(0,1);
  for (k in 1:Kb) yb[,k] ~ bernoulli_logit(z[,k]);
  
}

generated quantities{
  matrix[K,K] R = multiply_lower_tri_self_transpose(L_R);
  vector[K] full_sigma = append_row(sigma,rep_vector(1, Kb));
  matrix[K,K] Sigma = multiply_lower_tri_self_transpose(diag_pre_multiply(full_sigma,L_R));

}
