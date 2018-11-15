functions {
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
      else if(sign > 0)
        expterms[t] = exp(terms[t, 2] - c);
      else if(sign < 0)
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
  matrix[N_obs, N_pred] X[N_acc];

}
transformed data {
  real minRT = min(RT); 
  matrix[N_acc, N_acc - 1] helmert_m = helmert(N_acc);
  vector<lower = 0>[N_acc] s = rep_vector(1,N_acc); 
}
parameters {
  real<lower = 0, upper = minRT> t_0; // extra time for non-decision
  vector[N_pred] beta;  
  vector[N_acc - 1]  A_diff; //  difference between drifts
  real<lower = -min(helmert_m * A_diff)>  A_b; //  baseline A
  // The baseline shouldn't be smaller than the most negative value,
  // which is the first line of the helmert contrast (every beta * -1)
  vector[N_acc - 1]  k_diff; //  difference between drifts
  real<lower = -min(helmert_m * k_diff)>  k_b; //  baseline drift 
}
transformed parameters {
  vector<lower = 0>[N_acc] A = A_b + helmert_m * A_diff;  
  vector<lower = 0>[N_acc] k = k_b + helmert_m * k_diff; 
  vector<lower = 0>[N_acc] b = k + A; //distance to the threshold
}
model {
  matrix[N_acc, N_obs] nu;
  real log_lik[N_obs];

  A_b ~ normal(0, 2);
  A_diff ~ normal(0, .5);
  k_b ~ normal(0, 2);
  k_diff ~ normal(0, .5);
  beta ~ normal(0, 10);
  t_0 ~ lognormal(-1.7, .5); 

  for(a in 1:N_acc)
    nu[a, ] =  (X[a] * beta)';

  for (n in 1:N_obs) 
     log_lik[n] = tLBA(response[n], RT[n], A, b, nu[,n], s, t_0);
  
  target += sum(log_lik);
}

