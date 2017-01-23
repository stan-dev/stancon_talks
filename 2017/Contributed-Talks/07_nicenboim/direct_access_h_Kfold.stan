functions {
  real psi_max(vector u_psi, int[] subj, vector RT) {
    real psi_max;
    psi_max = positive_infinity();
    for (i in 1:num_elements(RT))
         psi_max = fmin(psi_max, log(RT[i]) - u_psi[subj[i]]);
    return (psi_max);
  }
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
  vector da_rng(vector theta, real P_b, real mu_da, real mu_b, real sigma, 
                real psi) {
    int orig_choice;
    int backtracking;
    vector[2] gen;
    orig_choice = categorical_rng(theta);
    backtracking = 0;
    if (orig_choice!=1) backtracking = bernoulli_rng(P_b);
    # Change the answer to 1 if there was backtracking:   
    gen[1] = backtracking ? 1 : orig_choice;
    { real mu; # it adds the mu_b if there is backtracking:
      mu = mu_da + (backtracking ? mu_b : 0);
      gen[2] = psi + lognormal_rng(mu, sigma);
    } 
    return(gen);  
  }
}
data {
  int<lower=0> N_obs; 
  int<lower=1> N_choices; 
  vector<lower=0>[N_obs] RT;
  int<lower=1,upper=N_choices> winner[N_obs];
  int<lower = 1> subj[N_obs];    
  int<lower = 1> N_subj;         
  int<lower = 1> item[N_obs];    
  int<lower = 1> N_item;  
  vector[N_obs] holdout; 
}
transformed data {
  real<lower=0> min_RT;
  real logmean_RT;
  min_RT = min(RT);
  logmean_RT = log(mean(RT));
}
parameters{
  real<lower=0> sigma;
  real<lower=0> mu_da_0raw;
  real<lower=0> mu_b_0;
  vector[N_choices-2] beta_incorrect; 
  real<lower=0> beta_added;
  vector<lower = 0> [N_choices - 1]  tau_u;     
  cholesky_factor_corr[N_choices - 1] L_u;  
  matrix[N_choices - 1, N_subj] z_u;
  vector<lower = 0> [2]  tau_u_RT;     
  cholesky_factor_corr[2] L_u_RT;  
  matrix[2, N_subj] z_u_RT;
  vector<lower = 0> [N_choices - 1]  tau_w;     
  cholesky_factor_corr[N_choices - 1] L_w;  
  matrix[N_choices - 1, N_item] z_w;
  vector<lower = 0> [2]  tau_w_RT;     
  cholesky_factor_corr[2] L_w_RT;  
  matrix[2, N_item] z_w_RT;
  real<lower=0,upper=1> P_b;
  vector[N_subj] u_psi;
  real<lower = 0> tau_psi; 
  real<upper = psi_max(u_psi, subj, RT) / logmean_RT> psi_0raw;
}
transformed parameters{
  real<lower=0> mu_da_0;
  vector[N_choices] beta_0;
  matrix[2, N_subj] u_RT; 
  matrix[N_choices, N_subj] u; 
  matrix[2, N_item] w_RT; 
  matrix[N_choices, N_item] w;
  real psi_0;
  u_RT = diag_pre_multiply(tau_u_RT, L_u_RT) * z_u_RT;   
  u[1:N_choices-1] = diag_pre_multiply(tau_u, L_u) * z_u;
  u[N_choices] = rep_row_vector(0,N_subj); 
  w_RT = diag_pre_multiply(tau_w_RT, L_w_RT) * z_w_RT;   
  w[1:N_choices-1] = diag_pre_multiply(tau_w, L_w) * z_w;
  w[N_choices] = rep_row_vector(0,N_item); 
  beta_0[1] = beta_added + fmax(max(beta_incorrect),0);
  beta_0[2:N_choices-1] = beta_incorrect;
  beta_0[N_choices] = 0;
  mu_da_0 = mu_da_0raw * logmean_RT;
  psi_0 = psi_0raw * logmean_RT;   
}
model {
  sigma ~ normal(0,2);
  beta_added ~ normal(0,2);
  beta_incorrect ~ normal(0,2);
  psi_0raw ~ normal(0, 1);
  tau_psi ~ normal(0, 1);
  u_psi ~ normal(0, tau_psi);
  to_vector(z_u_RT) ~ normal(0, 1); 
  to_vector(z_u) ~ normal(0, 1); 
  tau_u_RT ~ normal(0, 1);
  tau_u ~ normal(0, 1);
  L_u_RT ~ lkj_corr_cholesky(2.0);
  L_u ~ lkj_corr_cholesky(2.0);
  to_vector(z_w_RT) ~ normal(0, 1); 
  to_vector(z_w) ~ normal(0, 1); 
  tau_w_RT ~ normal(0, 1);
  tau_w ~ normal(0, 1);
  L_w_RT ~ lkj_corr_cholesky(2.0);
  L_w ~ lkj_corr_cholesky(2.0);
  P_b ~ beta(1,1);
  mu_da_0raw ~ normal(0,1);
  mu_b_0 ~ normal(0,2);
  for (n in 1:N_obs) {
    if(holdout[n]==0){
      real mu_da;
      real mu_b;
      vector[N_choices] beta;
      real psi;
      mu_da = mu_da_0 + u_RT[1,subj[n]] + w_RT[1,item[n]];
      mu_b = mu_b_0 + u_RT[2,subj[n]] + w_RT[2,item[n]];
      beta = beta_0 + u[,subj[n]] + w[,item[n]];
      psi = exp(psi_0 + u_psi[subj[n]]);
      target += da(winner[n], RT[n], beta, P_b, mu_da, mu_b, sigma, psi);
    }
  }
}
generated quantities {
  vector[N_choices] theta_0;
  matrix[N_choices-1, N_choices-1] Cor_u;
  matrix[N_choices-1, N_choices-1] Cor_w;  
  matrix[2, 2] Cor_u_RT;
  matrix[2, 2] Cor_w_RT;
  vector[N_obs] log_lik;
  theta_0 = softmax(beta_0);
  Cor_u = tcrossprod(L_u);  
  Cor_w = tcrossprod(L_w);  
  Cor_u_RT = tcrossprod(L_u_RT);  
  Cor_w_RT = tcrossprod(L_w_RT);  
  for (n in 1:N_obs) {
    real mu_da;
    real mu_b;
    vector[N_choices] beta;
    real psi;
    mu_da = mu_da_0 + u_RT[1,subj[n]] + w_RT[1,item[n]];
    mu_b = mu_b_0 + u_RT[2,subj[n]] + w_RT[2,item[n]];
    beta = beta_0 + u[,subj[n]] + w[,item[n]];
    psi = exp(psi_0 + u_psi[subj[n]]);
    log_lik[n] = da(winner[n], RT[n], beta, P_b, mu_da, mu_b, sigma, psi);
  }
}
