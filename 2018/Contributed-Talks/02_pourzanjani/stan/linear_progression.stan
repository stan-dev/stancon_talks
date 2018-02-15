functions {
  real[] f(real[] s, real[] a, real[] b, real[] c, real[] d) {

    int tot_obs;
    real ret[dims(s)[1]];
    tot_obs = dims(s)[1];

    //apply logistic function to each observation
    for(obs in 1:tot_obs) {
      ret[obs] = a[obs]/(1 + exp(-b[obs]*s[obs]+c[obs])) + d[obs];
    }
    
    return ret;
  }
}
data {
  int<lower=0> N;             //number of patients
  int<lower=0> K;             //number biomarkers
  int<lower=0> tot_obs;       //total number of observations
  
  int patient_idx[tot_obs];   //identify which patient the obs came from
  int biomarker_idx[tot_obs]; //identify which type of biomarker the obs is
  
  vector[tot_obs] age;        //age at the time of observation
  real y[tot_obs];            //obs value
}
parameters {
  
  //biomarker specific parameters
  real<upper=0> a_abeta;
  real<upper=0> a_hippo;
  real<lower=0> a_tau;
  
  real<lower=0> d_abeta;
  real<lower=0> d_hippo;
  real<lower=0> d_tau;
  
  real<lower=0> b_hippo;
  real<lower=0> b_mmse;
  real<lower=0> b_tau;
  
  real c_hippo;
  real c_mmse;
  real c_tau;
  
  real<lower=0> sigma_abeta;
  real<lower=0> sigma_hippo;
  real<lower=0> sigma_mmse;
  real<lower=0> sigma_tau;
  
  //individual specific disease paramters
  real<lower=0> gamma;
  vector<lower=0>[N] alpha;
  vector[N] beta;
}
transformed parameters {
  //use single arrays for simple indexing
  real a[K];
  real d[K];
  real<lower=0, upper=10> b[K];
  real c[K];
  real<lower=0> sigma[K];
  
  //MMSE does not need a and d because we know
  //it is designed to go from a scale of 0-30.
  //Because of the individual progressions, one
  //of the biomarker b,c parameters must be fixed
  //for identifiability. We choose ABETA
  a[1] = a_abeta;
  a[2] = a_hippo;
  a[3] = -30;
  a[4] = a_tau;
  
  d[1] = d_abeta;
  d[2] = d_hippo;
  d[3] = 30;
  d[4] = d_tau;
  
  b[1] = 1;
  b[2] = b_hippo;
  b[3] = b_mmse;
  b[4] = b_tau;
  
  c[1] = 0;
  c[2] = c_hippo;
  c[3] = c_mmse;
  c[4] = c_tau;
  
  sigma[1] = sigma_abeta;
  sigma[2] = sigma_hippo;
  sigma[3] = sigma_mmse;
  sigma[4] = sigma_tau;
}
model {
  
  real s[tot_obs] = to_array_1d(alpha[patient_idx] .* age + beta[patient_idx]);
  y ~ normal(f(s, a[biomarker_idx], b[biomarker_idx], c[biomarker_idx], d[biomarker_idx]), sigma[biomarker_idx]);
  
  a_abeta ~ normal(-110, 1);
  a_hippo ~ normal(-0.19,0.1);
  a_tau ~ normal(50,10);
  
  d_abeta ~ normal(245,1);
  d_hippo ~ normal(0.72,0.1);
  d_tau ~ normal(50,10);
  
  b_hippo ~ normal(1,1);
  b_mmse ~ normal(1,1);
  b_tau ~ normal(1,1);
  
  c_hippo ~ normal(0,20);
  c_mmse ~ normal(0,20);
  c_tau ~ normal(0,20);
  
  gamma ~ normal(0, 10);
  alpha ~ normal(1, gamma);
  beta ~ normal(0, 10);
}
generated quantities {
  
  //generate value of each of the biomarkers over time for each patient
  vector[81] grid;
  real shat[N,81];
  real fhat[N,K,81];
  
  real muhat[tot_obs];
  real yhat[tot_obs];
  real s[tot_obs] = to_array_1d(alpha[patient_idx] .* age + beta[patient_idx]);
  
  for(i in 1:81) grid[i] = -20 + 0.5*(i-1);
  for(n in 1:N) {
    shat[n,] = to_array_1d(alpha[n]*grid + rep_vector(beta[n],81));
    for(k in 1:K)
      fhat[n,k,] = f(to_array_1d(alpha[n]*grid+beta[n]),rep_array(a[k], 81),rep_array(b[k], 81),rep_array(c[k], 81),rep_array(d[k], 81));
  }
  
  muhat = f(s, a[biomarker_idx], b[biomarker_idx], c[biomarker_idx], d[biomarker_idx]);
  for(n in 1:tot_obs) yhat[n] = normal_rng(muhat[n], sigma[biomarker_idx][n]);
}
