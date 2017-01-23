functions{

  real[] effCptModel1(real t0, real t, real[] init, real amt, int cmt, int evid,
			  real CL, real Q, real V1, real V2, real ka, real ke0) {
    real k10;
    real k12;
    real k21;
    matrix[4, 4] K;

    real x[4];

    x = rep_array(0, 4);

    k10 = CL / V1;
    k12 = Q / V1;
    k21 = Q / V2;


    K = rep_matrix(0, 4, 4);
    
    K[1, 1] = -ka;
    K[2, 1] = ka;
    K[2, 2] = -(k10 + k12);
    K[2, 3] = k21;
    K[3, 2] = k12;
    K[3, 3] = -k21;
    K[4, 2] = ke0;
    K[4, 4] = -ke0;

    x = to_array_1d(matrix_exp((t - t0) * K) * to_vector(init));
   
    if (evid == 1) x[cmt] = x[cmt] + amt;

    return x;
  }

  matrix effCptModel(real[] time, real[] amt, int[] cmt, int[] evid, 
		     real CL, real Q, real V1, real V2, real ka, real ke0) {
    real init[4];
    matrix[size(time), 4] result;
    int nt;
    real t0;

    nt = size(time);

    init = rep_array(0, 4);
    t0 = time[1];
    for (i in 1:nt){
      init = effCptModel1(time[max(1, i - 1)], time[i], init, amt[i], cmt[i], evid[i],
			       CL, Q, V1, V2, ka, ke0);
      for (j in 1:4) result[i, j] = init[j];
      t0 = time[i];
    }
    return result;
  }
  
}

data {
  int<lower = 1> nt;
  int<lower = 1> nObs;
  int<lower = 1> iObs[nObs];
  real<lower = 0> amt[nt];
  int<lower = 1> cmt[nt];
  int<lower = 0> evid[nt];
  real<lower = 0> time[nt];
  vector<lower = 0>[nObs] cObs;
  vector[nObs] respObs;
}

transformed data {
  vector[nObs] logCObs;
  logCObs = log(cObs);
}

parameters {
  real<lower = 0> CL;
  real<lower = 0> Q;
  real<lower = 0> V1;
  real<lower = 0> V2;
  real<lower = 0> ka;
  real<lower = 0> ke0;
  real<lower = 0> EC50;
  real<lower = 0> sigma;
  real<lower = 0> sigmaResp;
}

transformed parameters {
  vector<lower = 0>[nt] cHat;
  vector<lower = 0>[nObs] cHatObs;
  vector<lower = 0>[nt] respHat;
  vector<lower = 0>[nObs] respHatObs;
  vector<lower = 0>[nt] ceHat;
  matrix[nt, 4] x;

  x = effCptModel(time, amt, cmt, evid,
                  CL, Q, V1, V2, ka, ke0);
                  
  cHat = 1000 * x[ ,2] ./ V1;
  ceHat = 1000 * x[ ,4] ./ V1;
  respHat = 100 * ceHat ./ (EC50 + ceHat);

  cHatObs = cHat[iObs];
  respHatObs = respHat[iObs];
  
  // for(i in 1:nObs){
  //   cHatObs[i] = cHat[iObs[i]];
  //   respHatObs[i] = respHat[iObs[i]];
  // }

}

model {
    CL ~ lognormal(log(10), 0.25);
    Q ~ lognormal(log(15), 0.5);
    V1 ~ lognormal(log(35), 0.25);
    V2 ~ lognormal(log(105), 0.5);
    ka ~ normal(0, 5);
    ke0 ~ normal(0, 2);
    EC50 ~ normal(0, 200);
    sigma ~ cauchy(0, 2);
    sigmaResp ~ cauchy(0, 5);

    logCObs ~ normal(log(cHatObs), sigma); 
    respObs ~ normal(respHatObs, sigmaResp); 
}

generated quantities {
  real cObsPred[nObs];
  real respObsPred[nObs];
  
  for(i in 1:nObs) {
    cObsPred[i] = exp(normal_rng(log(cHatObs[i]), sigma));
    respObsPred[i] = exp(normal_rng(log(respHatObs[i]), sigma));
  }

}

