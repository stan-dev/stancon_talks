functions {
  real psi_max(vector z_psi, real sd_psi, vector RT, int[] subj) {
    real psi_max = positive_infinity();
    vector[num_elements(z_psi)] r_psi = sd_psi * z_psi;
    vector[num_elements(RT)] logRT = log(RT);
    for (i in 1:num_elements(RT))
         psi_max = fmin(psi_max, logRT[i] - r_psi[subj[i]]);
    return psi_max;
  } 
  real nlcdf_approx(real x) {
    return log_inv_logit(0.07056 * x^3 + 1.5976* x);
  }  
  real nlpdf(real x) {
    return normal_lpdf(x | 0, 1);
  }
  real log_calc_exp(real[,] terms) {
    // Each term in the calculation is represented by an array of length two, t,
    // t[1] > 0 indicates that exp(t[2]) should be summed 
    // t[1] < 0 indicates that exp(t[2]) should be substracted 
    // t[1] == 0 indicates that exp(t[2]) should be ignored 
    real out;
    real c = max(terms[,2]);
    int nterms = dims(terms)[1];
    real expterms[nterms];
    if(dims(terms)[2] > 2) print("The extra dimensions in log_calc_exp will be ignored.")
    for(t in 1:nterms){
      real sign = terms[t, 1];
      if(sign == 0)
        expterms[t] = 0;
      else if(sign > 0 )
        expterms[t] = exp(terms[t, 2] - c);
      else if(sign < 0 )
        expterms[t] = -exp(terms[t, 2] - c);
    }
    out = c + log(sum(expterms));
    return out;
  }
  real tLBA_lpdf(real t, real A, real b, real nu, real s){
    real bAt = (b - A - t * nu)/(t * s);
    real bt = (b - t * nu)/(t * s);
    real prod_term = - log(A) - nlcdf_approx(nu/s); 
    real terms[4,2] = {{nu, log(fabs(nu)) + nlcdf_approx(bt)},
                      {1, log(s) + nlpdf(bAt)},
                      {-nu, log(fabs(nu)) + nlcdf_approx(bAt)},
                      {-1, log(s) + nlpdf(bt)}};

    real lpdf = log_calc_exp(terms) + prod_term; 
    if(is_nan(lpdf)) reject("Negative PDF!!")
    return lpdf;
  }    
  real tLBA_lcdf(real t, real A, real b, real nu, real s){
    real bmte = b - t * nu;
    real bmAte = bmte - A;
    real bt = bmte/(t * s);
    real bAt = bmAte/(t * s);
    real logs = log(s);
    real logt = log(t);
    real logA = log(A);
    real logb = log(b);
    real prod_term =  -nlcdf_approx(nu/s); 
    real lcdf;
    real terms[4, 2] = {{-bmAte, log(fabs(bmAte)) + nlcdf_approx(bAt)},
                        {-1, logt + logs + nlpdf(bAt)},
                        {bmte,log(fabs(bmte)) + nlcdf_approx(bt)},
                        {1, logt + logs + nlpdf(bt)}};

    real calc = log_calc_exp(terms) - logA ;
    // log(calc) should be between 0 and 1
    // We need lcdf =< 0  (exp(lcdf) is between 0 and 1)
    if(is_nan(calc))  
      // This happens when calc = log(x < 0)
      // This is because of approximation errors, and in this case
      // lcdf should can be set as 0
      lcdf = 0;
    else if(calc >= 0 )  
      // If calc on the raw scale is larger than 1,
      // also because of approx errors log(1 - x), with x > 1
      // then cdf should be 0 and lcdf is -Inf.
      lcdf = negative_infinity();
    else 
      lcdf = log1m_exp(calc) + prod_term;
      // lcdf can't be larger than zero
      lcdf = lcdf > 0 ?  0 : lcdf;
    return lcdf;
  } 
  real tLBA_lccdf(real t, real A, real b, real nu, real s){
    real lcdf = tLBA_lcdf(t|A, b, nu, s);
    real lccdf;
    if(lcdf == 0) {
      lccdf =  negative_infinity() ;
    } else if(lcdf == negative_infinity()) {
      lccdf =  0;
    } else{
      lccdf = log1m_exp(tLBA_lcdf(t | A, b, nu, s));
    }
                     
    if(is_nan(lccdf)) reject("lccdf is NaN, Parameters t:",t,", A:",A,", b:",b,", nu:",nu,", s:",s);                      
    return lccdf;
  }  
  real tLBA(int response, real RT, vector A, vector b, 
            vector nu, vector s, real t_0){
    
    real log_lik = 0;
    int N_acc = num_elements(nu);
    for(c in 1:N_acc){
      // Warnings:
      if (s[c]<=0)
        reject("s[",c,"] <= 0; found s[",c,"] = ", s[c])
      if (A[c]<=0)
        reject("A[",c,"] <= 0; found A[",c,"] = ", A[c])
      if (b[c]<=A[c])
        reject("b[",c,"] <= A[",c,"]; found A[",c,"] = " , A[c], "b[",c,"] = ",b[c])
      if (RT<=t_0)
        reject("RT <= t_0; found RT = ", RT, "t_0 = ",t_0)
      if (t_0<=0)
        reject("t_0 <= 0; found t_0 = ", t_0)
      if (RT<=0)
        reject("RT <= 0; found RT = ", RT)
      // end of Warnings
      if(c == response)
        log_lik = log_lik + tLBA_lpdf(RT - t_0|A[c], b[c], nu[c], s[c]);
        else 
        log_lik = log_lik + tLBA_lccdf(RT - t_0|A[c], b[c], nu[c], s[c]); 
    }
    return log_lik;
    }
matrix helmert(int levels) {
  // Matrix of helmert contrasts:
  matrix[levels, levels - 1] m = rep_matrix(-1, levels, levels - 1);
  real k = 0;
  for(j in 1:(levels - 1))
    for(i in 1:levels)
      if(j <= i - 2) 
        m[i, j] = 0;
      else if(j == i - 1) {
        k = k + 1;
        m[i, j] = k;
      }
  return m;    
  }   
}

data { 
  int N_acc;
  int<lower = 0> N_obs; 
  int<lower = 0, upper = N_acc> response[N_obs]; // response selected
  vector[N_obs] RT;
  int N_pred;
  int<lower=1> J_subj[N_obs]; 
  int<lower=1> N_subj; 
  matrix[N_obs, N_pred] X[N_acc];
}
transformed data {
  real minRT = min(RT); 
  matrix[N_acc, N_acc - 1] helmert_m = helmert(N_acc);
  vector<lower = 0>[N_acc] s = rep_vector(1,N_acc); 
}
parameters {
  // "random effects"
  vector<lower = 0>[N_pred] sd_nu_subj;  // group-level standard deviations 
  cholesky_factor_corr[N_pred] L_nu_subj;  
  matrix[N_pred, N_subj] z_nu_subj;  // unscaled group-level effects
  vector<lower = 0>[N_acc * 2] sd_Ak_subj;  // group-level standard deviations 
  cholesky_factor_corr[N_acc * 2] L_Ak_subj;  // unscaled group-level effects
  matrix[N_acc * 2, N_subj] z_Ak_subj;
  real<lower = 0> sd_psi_subj; 
  vector[N_subj] z_psi_subj;

  // "fixed effects"
  real<upper = psi_max(z_psi_subj, sd_psi_subj, RT, J_subj)> psi_0;
  vector[N_pred] beta;  
  vector[N_acc - 1]  A_diff; //  difference between drifts
  real<lower = -min(helmert_m * A_diff)>  A_b; //  baseline A
  // The baseline shouldn't be smaller than the most negative value,
  // which is the first line of the helmert contrast (every beta * -1)
  vector[N_acc - 1]  k_diff; //  difference between drifts
  real<lower = -min(helmert_m * k_diff)>  k_b; //  baseline drift 
//  baseline drift   
}
transformed parameters {
  vector<lower = 0>[N_acc] A = A_b + helmert_m * A_diff;  
  vector<lower = 0>[N_acc] k = k_b + helmert_m * k_diff; 
  vector<lower = 0>[N_acc] b = k + A; //distance to the threshold
  // For convienience, I leave r_nu_subj here: 
  matrix[N_pred, N_subj] r_nu_subj = 
              diag_pre_multiply(sd_nu_subj, L_nu_subj) * z_nu_subj;   
}
model {
  matrix[N_acc*2, N_subj] r_Ak_subj  = 
      diag_pre_multiply(sd_Ak_subj, L_Ak_subj) 
                * z_Ak_subj;
  vector[N_subj] r_psi_subj = sd_psi_subj * z_psi_subj;
  real log_lik[N_obs];
  vector[N_obs] t_0;
  matrix[N_obs, N_acc] t_nu;

  A_b ~ normal(0, 2);
  A_diff ~ normal(0, .5);
  k_b ~ normal(0, 2);
  k_diff ~ normal(0, .5);
  beta ~ normal(0, 10);
  psi_0 ~ normal(-1.7, .5); 

  sd_nu_subj ~ normal(0, 1);
  L_nu_subj ~ lkj_corr_cholesky(2.0);
  to_vector(z_nu_subj) ~ normal(0, 1);
  sd_Ak_subj ~ normal(0, 1);
  L_Ak_subj ~ lkj_corr_cholesky(2.0);
  to_vector(z_Ak_subj) ~ normal(0, 1);
  sd_psi_subj ~ normal(0, 1);
  z_psi_subj ~ normal(0, 1);


  for(a in 1:N_acc)
    t_nu[, a] =  X[a] * beta;

  t_0 = exp(r_psi_subj[J_subj] + psi_0);


  for (n in 1:N_obs) {
    vector[N_acc] A_n;
    vector[N_acc] b_n;
    vector[N_acc] nu_n;                
    for(a in 1:N_acc) {
      real  r_acc = 0;
      A_n[a] = A[a] * exp(r_Ak_subj[a, J_subj[n]]);
      b_n[a] = k[a] * exp(r_Ak_subj[N_acc + a, J_subj[n]]) + A_n[a];
      nu_n[a] = t_nu[n, a];
      for(p in 1:(N_pred))
        r_acc = r_acc + X[a][n,p] * r_nu_subj[p, J_subj[n]];
      nu_n[a] = nu_n[a] + r_acc;   
      }
    
    log_lik[n] = tLBA(response[n], RT[n], A_n, b_n, nu_n, s, t_0[n]);
  }
  target += sum(log_lik);
}
generated quantities {
  matrix[N_pred, N_pred] Corr_nu = L_nu_subj * L_nu_subj';
  matrix[N_acc * 2,N_acc * 2] Corr_Ak = L_Ak_subj * L_Ak_subj';
  real t_0 = exp(psi_0);
}

