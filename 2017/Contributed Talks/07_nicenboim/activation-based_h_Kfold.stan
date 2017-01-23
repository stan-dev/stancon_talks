functions {
    vector mean_of_X_by_Y(vector X, int[] Y){
    int N_rows;
    int N_groups;
    N_rows = num_elements(X);
    if(N_rows!= num_elements(Y)) 
      reject("X and Y don't have the same length");
    N_groups = max(Y);
    { # matrix with a column for each group of Y, 
      #and 1 if Y belong to the group:
      matrix[N_rows, N_groups] matrix_1; 
      for(r in 1:N_rows)
        for(g in 1:N_groups)
        matrix_1[r,g] =  (g ==  Y[r]);  
      return(((X' * matrix_1)   // sum of Xs per group
        // divided by matrix_1^T * matrix_1 (number of times each group appears):
      / crossprod(matrix_1))');     
    }
  }
  real psi_max(vector u_psi, int[] subj, vector RT) {
    real psi_max;
    psi_max = positive_infinity();
    for (i in 1:num_elements(RT))
         psi_max = fmin(psi_max, log(RT[i]) - u_psi[subj[i]]);
    return (psi_max);
  }
    real race(int winner, real RT, vector alpha, real b, real sigma, real psi) {
    real log_lik;
    int N_choices;
    N_choices = num_elements(alpha);
    log_lik = 0;
    for (c in 1:N_choices)
        if (c == winner)
          log_lik = log_lik + lognormal_lpdf(RT-psi|b-alpha[c],sigma);
        else 
          log_lik = log_lik + lognormal_lccdf(RT-psi|b-alpha[c], sigma); 
    return(log_lik);
    }
  vector race_rng(vector alpha, real b, real sigma,
                   real psi) {
    int N_choices;
    real curr_accum_RT; 
    real fastest_accum_RT; 
    vector[2] gen;
    N_choices = num_elements(alpha);
    fastest_accum_RT = positive_infinity();
    for (c in 1:N_choices) {
      curr_accum_RT = psi + lognormal_rng(b - alpha[c], sigma);
      if (curr_accum_RT < fastest_accum_RT){ #this accumulator is faster
        fastest_accum_RT = curr_accum_RT;
        gen[1] = c;
        gen[2] = curr_accum_RT;
      }
    } 
    return(gen);  
  }
}
data { 
  int<lower = 0> N_obs; 
  int<lower = 1> N_choices; 
  int<lower = 1, upper = N_choices> winner[N_obs];
  vector<lower = 0>[N_obs] RT;
  int<lower = 1> N_subj;         
  int<lower = 1> subj[N_obs];    
  int<lower = 1> item[N_obs];    
  int<lower = 1> N_item;        
  vector[N_obs] holdout; 
}
transformed data {
  real b; 
  real min_RT;
  real logmean_RT;
  vector[N_choices] logmean_RT_w;
  real N_holdout;
  b = 10;
  min_RT = min(RT);
  logmean_RT = log(mean(RT));
  logmean_RT_w = log(mean_of_X_by_Y(RT, winner));
}
parameters{
  vector[N_choices] alpha_0raw; 
  vector<lower = 0> [N_choices]  tau_u;     
  cholesky_factor_corr[N_choices] L_u;  
  matrix[N_choices, N_subj] z_u;
  vector<lower = 0> [N_choices]  tau_w;     
  cholesky_factor_corr[N_choices] L_w;  
  matrix[N_choices, N_item] z_w;
  real<lower = 0> sigma;
  vector[N_subj] u_psi;
  real<lower = 0> tau_psi; 
  real<upper = psi_max(u_psi, subj, RT) / logmean_RT> psi_0raw;
}
transformed parameters {
  real psi_0;
  matrix[N_choices, N_subj] u; 
  matrix[N_choices, N_item] w; 
  vector[N_choices] alpha_0; 

  // Optimization through Cholesky Fact:
  u = diag_pre_multiply(tau_u, L_u) //matrix[N_choices,N_choices]
      * z_u;  
  w = diag_pre_multiply(tau_w, L_w) //matrix[N_choices,N_choices]
      * z_w;
  // Unit-scaling
  alpha_0 = b - alpha_0raw .* logmean_RT_w;
  psi_0 = psi_0raw * logmean_RT;   
}
model {
  alpha_0raw ~ normal(0, 1);
  tau_u ~ normal(0, 1);
  L_u ~ lkj_corr_cholesky(2.0);
  to_vector(z_u) ~ normal(0, 1); 
  tau_w ~ normal(0, 1);
  L_w ~ lkj_corr_cholesky(2.0);
  to_vector(z_w) ~ normal(0, 1); 
  sigma ~ normal(0, 2);
  psi_0raw ~ normal(0, 1);
  tau_psi ~ normal(0, 1);
  u_psi ~ normal(0, tau_psi);
  for (n in 1:N_obs) {
    if(holdout[n]==0){
      vector[N_choices] alpha; 
      real psi;
      alpha = alpha_0 + u[, subj[n]] + w[, item[n]];
      psi = exp(psi_0 + u_psi[subj[n]]);
      target += race(winner[n], RT[n], alpha, b, sigma, psi);
    }
  }
}
generated quantities {
  matrix[N_choices, N_choices] Cor_u;
  matrix[N_choices, N_choices] Cor_w;
  vector[N_obs] log_lik;
  Cor_u = tcrossprod(L_u);  
  Cor_w = tcrossprod(L_w);  
  for (n in 1:N_obs) {
    vector[N_choices] alpha; 
    real psi;
    alpha = alpha_0 + u[, subj[n]] + w[, item[n]];
    psi = exp(psi_0 + u_psi[subj[n]]);
    log_lik[n] = race(winner[n], RT[n], alpha, b, sigma, psi);
  }
}
