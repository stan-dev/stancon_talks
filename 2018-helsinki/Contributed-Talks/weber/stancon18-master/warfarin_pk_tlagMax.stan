functions {
  #include "utils.stan"
  #include "gamma2_overdisp.stan"
  #include "models.stan"
}
data {
  int J;
  real<lower=0> tlag_max;
  vector<lower=0>[J] dose;
  int N[J];
  vector<lower=0>[sum(N)] dv;
  vector<lower=tlag_max>[sum(N)] time;
  vector<lower=0>[J] weight;
  real<lower=0> ref;
}
transformed data {
  int Ndv = sum(N);
  row_vector[J] LSweight = to_row_vector(log(weight) - log(70.));
  vector[J] ldose = log(dose);
  real ltlag_max = log(tlag_max);
  int sidx[J+1] = make_slice_index(N);
  real refSq = square(ref);
}
parameters {
  // log population parameters of 1-tlag, 2-ka, 3-Cl, 4-V
  vector[4] theta;

  // subject specific random effects
  matrix[4,J] Eta;
  vector<lower=0>[4] sigma_eta;

  real<lower=0> sigma_y;
  real<lower=0> kappa;
}
transformed parameters {
  vector[Ndv] mu;
  matrix[4,J] Eta_cov;

  // lag time is fitted as fraction in logit space relative to the
  // maximal tlag
  Eta_cov[1] = log_inv_logit(theta[1] + sigma_eta[1] * Eta[1]) + ltlag_max;
  Eta_cov[2] = theta[2] + sigma_eta[2] * Eta[2];
  // apply to CL and V the weight covariate effects
  Eta_cov[3] = theta[3] + sigma_eta[3] * Eta[3] +  0.75 * LSweight;
  Eta_cov[4] = theta[4] + sigma_eta[4] * Eta[4] +         LSweight;
  
  // loop over each subject and evaluate PK model
  for(j in 1:J) {
    int start = sidx[j];
    int end   = sidx[j+1]-1;
    mu[start:end] = exp(pk_1cmt_oral_tlagMax(time[start:end], ldose[j], Eta_cov[1,j], Eta_cov[2,j], Eta_cov[3,j], Eta_cov[4,j])) + 1E-5;
  }
}
model {
  theta[1] ~ normal(0, 2);
  theta[2] ~ normal(log(1),   log(2)/1.96);
  theta[3] ~ normal(log(0.1), log(10)/1.96);
  theta[4] ~ normal(log(10),  log(10)/1.96);

  sigma_eta ~ normal(0, 0.5);
 
  to_vector(Eta) ~ normal(0, 1);

  sigma_y ~ normal(0, 2);
  kappa ~ gamma(0.2, 0.2);

  dv ~ gamma2_overdisp(mu, sigma_y, kappa * refSq);
}
generated quantities {
  real tlag = inv_logit(theta[1]) * tlag_max;
  real ka = exp(theta[2]);
  real CL = exp(theta[3]);
  real V = exp(theta[4]);
  vector[Ndv] log_lik;
  vector[Ndv] ypred;
  vector[Ndv] ypred_cond;
  vector[J] log_lik_patient;
  real sigma_y_ref = sqrt(square(sigma_y) + 1.0/kappa);

  for(o in 1:Ndv) {
    log_lik[o] = gamma2_overdisp_lpdf( dv[o:o] | mu[o:o], sigma_y, kappa * refSq);
    ypred_cond[o] = gamma2_overdisp_rng(mu[o], sigma_y, kappa * refSq);
  }

  for(j in 1:J) {
    int start = sidx[j];
    int end   = sidx[j+1]-1;
    vector[4] eta = multi_normal_cholesky_rng(theta, diag_matrix(sigma_eta));
    vector[4] eta_cov;
    vector[N[j]] mu_pred;

    eta_cov[1] = log_inv_logit(eta[1]) + ltlag_max;
    eta_cov[2] = eta[2];
    // apply to CL and V the weight covariate effects
    eta_cov[3] = eta[3] +  0.75 * LSweight[j];
    eta_cov[4] = eta[4] +         LSweight[j];

    mu_pred = exp(pk_1cmt_oral_tlagMax(time[start:end], ldose[j], eta_cov[1], eta_cov[2], eta_cov[3], eta_cov[4])) + 1E-5;

    log_lik_patient[j] = gamma2_overdisp_lpdf(dv[start:end]| mu_pred, sigma_y, kappa * refSq);
    for(k in 1:N[j]){
      ypred[start + k - 1] = gamma2_overdisp_rng(mu_pred[k], sigma_y, kappa * refSq);
    }
  }
}
