functions {
  real I1(real x) { 
	 real ret = 0;
	 if(x < -20) return(0);

	 if(x >= -20 && x < -20) return(ret + 0 - 0);
	 else ret = ret + 0;

	 if(x >= -20 && x < -20) return(ret + 0 - 0);
	 else ret = ret + 0;

	 if(x >= -20 && x < -20) return(ret + 0 - 0);
	 else ret = ret + 0;

	 if(x >= -20 && x < -16) return(ret + -0.9999999999 * x^4/256 + -47.99999999 * x^3/192 + -767.9999999 * x^2/128 + -4095.9999998976 * x/64 - 254.999999926333);
	 else ret = ret + 1.00000006033338;

	 if(x >= -16) return(1);

	 return(ret);
}

real I2(real x) { 
	 real ret = 0;
	 if(x < -20) return(0);

	 if(x >= -20 && x < -20) return(ret + 0 - 0);
	 else ret = ret + 0;

	 if(x >= -20 && x < -20) return(ret + 0 - 0);
	 else ret = ret + 0;

	 if(x >= -20 && x < -16) return(ret + 448 * x^4/131072 + 22272 * x^3/98304 + 365568 * x^2/65536 + 1986560 * x/32768 - -246.875);
	 else ret = ret + 0.875;

	 if(x >= -16 && x < -12) return(ret + -0.9999999999 * x^4/2048 + -36 * x^3/1536 + -431.9999999 * x^2/1024 + -1727.9999999568 * x/512 - 10.00000002685);
	 else ret = ret + 0.124999987212497;

	 if(x >= -12) return(1);

	 return(ret);
}

real I3(real x) { 
	 real ret = 0;
	 if(x < -20) return(0);

	 if(x >= -20 && x < -20) return(ret + 0 - 0);
	 else ret = ret + 0;

	 if(x >= -20 && x < -16) return(ret + -1408 * x^4/1179648 + -75264 * x^3/884736 + -1320960 * x^2/589824 + -7577600 * x/294912 - 107.638888888889);
	 else ret = ret + 0.361111111111086;

	 if(x >= -16 && x < -12) return(ret + 224 * x^4/294912 + 8832 * x^3/221184 + 112128 * x^2/147456 + 464896 * x/73728 - -20);
	 else ret = ret + 0.583333333333329;

	 if(x >= -12 && x < -8) return(ret + -0.9999999999 * x^4/4608 + -24 * x^3/3456 + -192 * x^2/2304 + -511.9999999872 * x/1152 - 0.83333333365);
	 else ret = ret + 0.0555555552388891;

	 if(x >= -8) return(1);

	 return(ret);
}

real I4(real x) { 
	 real ret = 0;
	 if(x < -20) return(0);

	 if(x >= -20 && x < -16) return(ret + 1 * x^4/6144 + 60 * x^3/4608 + 1200 * x^2/3072 + 8000 * x/1536 - -26.0416666666667);
	 else ret = ret + 0.0416666666666714;

	 if(x >= -16 && x < -12) return(ret + -3 * x^4/6144 + -132 * x^3/4608 + -1872 * x^2/3072 + -8383.9999999999 * x/1536 - 16.6666666666656);
	 else ret = ret + 0.458333333333599;

	 if(x >= -12 && x < -8) return(ret + 3 * x^4/6144 + 84 * x^3/4608 + 720 * x^2/3072 + 1984 * x/1536 - -3.125);
	 else ret = ret + 0.458333333333332;

	 if(x >= -8 && x < -4) return(ret + -0.9999999999 * x^4/6144 + -12 * x^3/4608 + -47.99999999 * x^2/3072 + -63.9999999984 * x/1536 - 2.66666633219614e-10);
	 else ret = ret + 0.0416666664520834;

	 if(x >= -4) return(1);

	 return(ret);
}

real I5(real x) { 
	 real ret = 0;
	 if(x < -16) return(0);

	 if(x >= -16 && x < -12) return(ret + 1 * x^4/6144 + 48 * x^3/4608 + 768 * x^2/3072 + 4096 * x/1536 - -10.6666666666667);
	 else ret = ret + 0.0416666666666643;

	 if(x >= -12 && x < -8) return(ret + -3 * x^4/6144 + -96 * x^3/4608 + -960 * x^2/3072 + -2816 * x/1536 - 2.875);
	 else ret = ret + 0.458333333333332;

	 if(x >= -8 && x < -4) return(ret + 3 * x^4/6144 + 48 * x^3/4608 + 192 * x^2/3072 + 256 * x/1536 - -0.666666666666666);
	 else ret = ret + 0.458333333333333;

	 if(x >= -4 && x < 0) return(ret + -0.9999999999 * x^4/6144 - -0.0416666666625);
	 else ret = ret + 0.0416666666625;

	 if(x >= 0) return(1);

	 return(ret);
}

real I6(real x) { 
	 real ret = 0;
	 if(x < -12) return(0);

	 if(x >= -12 && x < -8) return(ret + 1 * x^4/6144 + 36 * x^3/4608 + 432 * x^2/3072 + 1728 * x/1536 - -3.375);
	 else ret = ret + 0.0416666666666661;

	 if(x >= -8 && x < -4) return(ret + -3 * x^4/6144 + -60 * x^3/4608 + -336 * x^2/3072 + -30720 * x/147456 - -0.666666666666666);
	 else ret = ret + 0.458333333333333;

	 if(x >= -4 && x < 0) return(ret + 3 * x^4/6144 + 12 * x^3/4608 + -48 * x^2/3072 + 64 * x/1536 - -0.458333333333333);
	 else ret = ret + 0.458333333333333;

	 if(x >= 0 && x < 4) return(ret + -1 * x^4/6144 + 12 * x^3/4608 + -48 * x^2/3072 + 64 * x/1536 - 0);
	 else ret = ret + 0.0416666666666667;

	 if(x >= 4) return(1);

	 return(ret);
}

real I7(real x) { 
	 real ret = 0;
	 if(x < -8) return(0);

	 if(x >= -8 && x < -4) return(ret + 1 * x^4/6144 + 24 * x^3/4608 + 192 * x^2/3072 + 512 * x/1536 - -0.666666666666667);
	 else ret = ret + 0.0416666666666666;

	 if(x >= -4 && x < 0) return(ret + -3 * x^4/6144 + -24 * x^3/4608 + 256 * x/1536 - -0.458333333333333);
	 else ret = ret + 0.458333333333333;

	 if(x >= 0 && x < 4) return(ret + 3 * x^4/6144 + -24 * x^3/4608 + 256 * x/1536 - 0);
	 else ret = ret + 0.458333333333333;

	 if(x >= 4 && x < 8) return(ret + -1 * x^4/6144 + 24 * x^3/4608 + -192 * x^2/3072 + 512 * x/1536 - 0.625);
	 else ret = ret + 0.0416666666666666;

	 if(x >= 8) return(1);

	 return(ret);
}

real I8(real x) { 
	 real ret = 0;
	 if(x < -4) return(0);

	 if(x >= -4 && x < 0) return(ret + 1 * x^4/6144 + 12 * x^3/4608 + 48 * x^2/3072 + 64 * x/1536 - -0.0416666666666667);
	 else ret = ret + 0.0416666666666667;

	 if(x >= 0 && x < 4) return(ret + -3 * x^4/6144 + 12 * x^3/4608 + 48 * x^2/3072 + 64 * x/1536 - 0);
	 else ret = ret + 0.458333333333333;

	 if(x >= 4 && x < 8) return(ret + 3 * x^4/6144 + -60 * x^3/4608 + 336 * x^2/3072 + -320 * x/1536 - 0.208333333333333);
	 else ret = ret + 0.458333333333333;

	 if(x >= 8 && x < 12) return(ret + -1 * x^4/6144 + 36 * x^3/4608 + -432 * x^2/3072 + 1728 * x/1536 - 3.33333333333333);
	 else ret = ret + 0.0416666666666661;

	 if(x >= 12) return(1);

	 return(ret);
}

real I9(real x) { 
	 real ret = 0;
	 if(x < 0) return(0);

	 if(x >= 0 && x < 4) return(ret + 1 * x^4/6144 - 0);
	 else ret = ret + 0.0416666666666667;

	 if(x >= 4 && x < 8) return(ret + -3 * x^4/6144 + 48 * x^3/4608 + -192 * x^2/3072 + 256 * x/1536 - 0.208333333333333);
	 else ret = ret + 0.458333333333333;

	 if(x >= 8 && x < 12) return(ret + 3 * x^4/6144 + -96 * x^3/4608 + 960 * x^2/3072 + -2816 * x/1536 - -3.33333333333333);
	 else ret = ret + 0.458333333333332;

	 if(x >= 12 && x < 16) return(ret + -1 * x^4/6144 + 48 * x^3/4608 + -768 * x^2/3072 + 4096 * x/1536 - 10.625);
	 else ret = ret + 0.0416666666666643;

	 if(x >= 16) return(1);

	 return(ret);
}

real I10(real x) { 
	 real ret = 0;
	 if(x < 4) return(0);

	 if(x >= 4 && x < 8) return(ret + 1 * x^4/6144 + -12 * x^3/4608 + 48 * x^2/3072 + -64 * x/1536 - -0.0416666666666667);
	 else ret = ret + 0.0416666666666667;

	 if(x >= 8 && x < 12) return(ret + -3 * x^4/6144 + 84 * x^3/4608 + -720 * x^2/3072 + 1984 * x/1536 - 2.66666666666667);
	 else ret = ret + 0.458333333333332;

	 if(x >= 12 && x < 16) return(ret + 3 * x^4/6144 + -132 * x^3/4608 + 1872 * x^2/3072 + -8383.9999999999 * x/1536 - -17.1249999999992);
	 else ret = ret + 0.458333333333599;

	 if(x >= 16 && x < 20) return(ret + -1 * x^4/6144 + 60 * x^3/4608 + -1200 * x^2/3072 + 8000 * x/1536 - 26);
	 else ret = ret + 0.0416666666666714;

	 if(x >= 20) return(1);

	 return(ret);
}

real I11(real x) { 
	 real ret = 0;
	 if(x < 8) return(0);

	 if(x >= 8 && x < 12) return(ret + 1 * x^4/4608 + -24 * x^3/3456 + 192 * x^2/2304 + -512 * x/1152 - -0.888888888888889);
	 else ret = ret + 0.0555555555555558;

	 if(x >= 12 && x < 16) return(ret + -224 * x^4/294912 + 8832 * x^3/221184 + -112128 * x^2/147456 + 464896 * x/73728 - 19.4166666666667);
	 else ret = ret + 0.583333333333329;

	 if(x >= 16 && x < 20) return(ret + 1408 * x^4/1179648 + -75264 * x^3/884736 + 1320960 * x^2/589824 + -7577600 * x/294912 - -108);
	 else ret = ret + 0.361111111111086;

	 if(x >= 20 && x < 20) return(ret + 0 - 0);
	 else ret = ret + 0;

	 if(x >= 20) return(1);

	 return(ret);
}

real I12(real x) { 
	 real ret = 0;
	 if(x < 12) return(0);

	 if(x >= 12 && x < 16) return(ret + 1 * x^4/2048 + -36 * x^3/1536 + 432 * x^2/1024 + -1728 * x/512 - -10.125);
	 else ret = ret + 0.125;

	 if(x >= 16 && x < 20) return(ret + -448 * x^4/131072 + 22272 * x^3/98304 + -365568 * x^2/65536 + 1986560 * x/32768 - 246);
	 else ret = ret + 0.875;

	 if(x >= 20 && x < 20) return(ret + 0 - 0);
	 else ret = ret + 0;

	 if(x >= 20 && x < 20) return(ret + 0 - 0);
	 else ret = ret + 0;

	 if(x >= 20) return(1);

	 return(ret);
}

real I13(real x) { 
	 real ret = 0;
	 if(x < 16) return(0);

	 if(x >= 16 && x < 20) return(ret + 1 * x^4/256 + -48 * x^3/192 + 768 * x^2/128 + -4095.9999999999 * x/64 - -255.999999999975);
	 else ret = ret + 1.00000000000637;

	 if(x >= 20 && x < 20) return(ret + 0 - 0);
	 else ret = ret + 0;

	 if(x >= 20 && x < 20) return(ret + 0 - 0);
	 else ret = ret + 0;

	 if(x >= 20 && x < 20) return(ret + 0 - 0);
	 else ret = ret + 0;

	 if(x >= 20) return(1);

	 return(ret);
}

real[] ispline(vector x, matrix c) {

	 int N;
	 real ret[dims(x)[1]];
	 N = dims(x)[1];

	 for(n in 1:N) ret[n] = c[n][1]*I1(x[n]) + c[n][2]*I2(x[n]) + c[n][3]*I3(x[n]) + c[n][4]*I4(x[n]) + c[n][5]*I5(x[n]) + c[n][6]*I6(x[n]) + c[n][7]*I7(x[n]) + c[n][8]*I8(x[n]) + c[n][9]*I9(x[n]) + c[n][10]*I10(x[n]) + c[n][11]*I11(x[n]) + c[n][12]*I12(x[n]) + c[n][13]*I13(x[n]);

	 return(ret);
}

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
  int<lower=0> N; //number of patients
  int<lower=0> K; //number biomarkers
  int<lower=0> tot_obs;
  
  int patient_idx[tot_obs];
  int biomarker_idx[tot_obs];
  
  vector[tot_obs] age;
  real y[tot_obs];
}
parameters {
  
  //biomarker specific parameters
  real<upper=0> a_abeta;
  real<upper=0> a_hippo;
  real<lower=0> a_tau;
  
  real<lower=0> d_abeta;
  real<lower=0> d_hippo;
  real<lower=0> d_tau;
  
  real<lower=0> b_abeta;
  real<lower=0> b_hippo;
  real<lower=0> b_mmse;
  real<lower=0> b_tau;
  
  real c_abeta;
  real c_hippo;
  real c_mmse;
  real c_tau;
  
  //we now have a patient specific sigma_abeta
  real<lower=0> sigma_abeta[N];
  real<lower=0> sigma_hippo;
  real<lower=0> sigma_mmse;
  real<lower=0> sigma_tau;
  
  //each patient has a positive linear combination, beta,
  //of the I-splines. Each of these coefficients, respectively
  //follow a hierarchical prior
  matrix<lower=0>[N,13] beta_raw;
  real<lower=0> gamma[13];
}
transformed parameters {
  real a[K];
  real d[K];
  real<lower=0, upper=10> b[K];
  real c[K];
  real<lower=0> sigma[K];
  matrix<lower=0>[N,13] beta;
  
  //each observation will now have a sigma
  //that depends on the biomarker being measured
  //and the specific patient being measured
  real<lower=0> sigma_obs[tot_obs];
  
  a[1] = a_abeta;
  a[2] = a_hippo;
  a[3] = -30;
  a[4] = a_tau;
  
  d[1] = d_abeta;
  d[2] = d_hippo;
  d[3] = 30;
  d[4] = d_tau;
  
  //with s constrained to be positive c_abeta and b_abeta identifiable
  //so we no longer set them to fixed constants
  b[1] = b_abeta;
  b[2] = b_hippo;
  b[3] = b_mmse;
  b[4] = b_tau;
  
  c[1] = c_abeta;
  c[2] = c_hippo;
  c[3] = c_mmse;
  c[4] = c_tau;
  
  sigma[1] = 1;
  sigma[2] = sigma_hippo;
  sigma[3] = sigma_mmse;
  sigma[4] = sigma_tau;
  
  //fill sigma_obs
  for(obs in 1:tot_obs) {
    if(biomarker_idx[obs] == 1) sigma_obs[obs] = sigma_abeta[patient_idx[obs]];
    else sigma_obs[obs] = sigma[biomarker_idx[obs]]; 
  }
  
  for(i in 1:13) beta[,i] = beta_raw[,i]*gamma[i];
  
}
model {
  real s[tot_obs] = ispline(age, beta[patient_idx]);
  y ~ normal(f(s, a[biomarker_idx], b[biomarker_idx], c[biomarker_idx], d[biomarker_idx]), sigma_obs);

  sigma_abeta ~ normal(0,30);
  sigma_hippo ~ normal(0,0.01);
  sigma_mmse ~ normal(0,1);

  //hierarchical prior
  gamma ~ normal(0,10);
  for(i in 1:13) beta_raw[,i] ~ normal(0,1);
  
  a_abeta ~ normal(-110, 20);
  a_hippo ~ normal(-0.19,0.1);
  a_tau ~ normal(50,20);
  
  d_abeta ~ normal(245,20);
  d_hippo ~ normal(0.72,0.1);
  d_tau ~ normal(50,20);
  
  b_abeta ~ normal(1,1);
  b_hippo ~ normal(1,1);
  b_mmse ~ normal(1,1);
  b_tau ~ normal(1,1);
  
  c_abeta ~ normal(1,1);
  c_hippo ~ normal(0,20);
  c_mmse ~ normal(0,20);
  c_tau ~ normal(0,20);
}
generated quantities {
  vector[81] grid;
  real fhat[N,K,81];
  real shat[N,81];
  
  real muhat[tot_obs];
  real yhat[tot_obs];
  real s[tot_obs] = ispline(age, beta[patient_idx]);
  
  for(i in 1:81) grid[i] = -20 + 0.5*(i-1);
  for(n in 1:N) {
    shat[n,] = ispline(grid, beta[rep_array(n, 81)]);
    for(k in 1:K)
      fhat[n,k,] = f(ispline(grid, beta[rep_array(n, 81)]),rep_array(a[k], 81),rep_array(b[k], 81),rep_array(c[k], 81),rep_array(d[k], 81));
  }
  
  muhat = f(s, a[biomarker_idx], b[biomarker_idx], c[biomarker_idx], d[biomarker_idx]);
  for(n in 1:tot_obs) yhat[n] = normal_rng(muhat[n], sigma_obs[n]);
}