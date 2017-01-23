functions {
  vector twoCptModel1(real dt, vector init, real amt, int cmt, int evid,
		    real CL, real Q, real V1, real V2, real ka) {
    int nCmt = 3;
    vector[nCmt] x;
    real k10;
    real k12;
    real k21;
    real ksum;
    vector[nCmt] a;
    vector[nCmt] alpha;

    k10 = CL / V1;
    k12 = Q / V1;
    k21 = Q / V2;
    ksum = k10 + k12 + k21;
    alpha[1] = (ksum + sqrt(ksum * ksum - 4.0 * k10 * k21))/2.0;
    alpha[2] = (ksum - sqrt(ksum * ksum - 4.0 * k10 * k21))/2.0;
    alpha[3] = ka;

    x = rep_vector(0.0, nCmt);

    if (init[1] != 0.0) {
      x[1] = init[1] * exp(-alpha[3] * dt);
      a[1] = ka * (k21 - alpha[1]) / ((ka - alpha[1]) * (alpha[2] - alpha[1]));
      a[2] = ka * (k21 - alpha[2]) / ((ka - alpha[2]) * (alpha[1] - alpha[2]));
      a[3] = -(a[1] + a[2]);
      x[2] = init[1] * sum(a .* exp(-alpha * dt));
      a[1] = ka * k12 / ((ka - alpha[1]) * (alpha[2] - alpha[1]));
      a[2] = ka * k12 / ((ka - alpha[2]) * (alpha[1] - alpha[2]));
      a[3] = -(a[1] + a[2]);
      x[3] = init[1] * sum(a .* exp(-alpha * dt));
    }
    
    if (init[2] != 0) {
      a[1] = (k21 - alpha[1]) / (alpha[2] - alpha[1]);
      a[2] = (k21 - alpha[2]) / (alpha[1] - alpha[2]);
      x[2] = x[2] + init[2] * sum(segment(a, 1, 2) .* exp(-segment(alpha, 1, 2) * dt));
      a[1] = k12 / (alpha[2] - alpha[1]);
      a[2] = -a[1];
      x[3] = x[3] + init[2] * sum(segment(a, 1, 2) .* exp(-segment(alpha, 1, 2) * dt));
    }

    if (init[3] != 0) {
      a[1] = k21 / (alpha[2] - alpha[1]);
      a[2] = -a[1];
      x[2] = x[2] + init[3] * sum(segment(a, 1, 2) .* exp(-segment(alpha, 1, 2) * dt));
      a[1] = (k10 + k12 - alpha[1]) / (alpha[2] - alpha[1]);
      a[2] = (k10 + k12 - alpha[2]) / (alpha[1] - alpha[2]);
      x[3] = x[3] + init[3] * sum(segment(a, 1, 2) .* exp(-segment(alpha, 1, 2) * dt));
    }

    if (evid == 1) x[cmt] = x[cmt] + amt;

    return x;
  }

  matrix twoCptModel(real[] time, real[] amt, int[] cmt, int[] evid, 
		     real CL, real Q, real V1, real V2, real ka) {
    int nCmt = 3;
    vector[nCmt] init;
    real dt;
    real t0;
    matrix[size(time), nCmt] result;
    int nt;

    nt = size(time);

    init = rep_vector(0, nCmt);
    t0 = time[1];
    for (i in 1:nt) {
      dt = time[i] - t0;
      init = twoCptModel1(dt, init, amt[i], cmt[i], evid[i],
			   CL, Q, V1, V2, ka);
      for(j in 1:nCmt) result[i, j] = init[j];
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
}

transformed data {
  vector[nObs] logCObs;
  int nCmt;

  logCObs = log(cObs);
  nCmt = 3; ## Fixed. The code specifically handles models with 3 compartments.
}

parameters {
  real<lower = 0> CL;
  real<lower = 0> Q;
  real<lower = 0> V1;
  real<lower = 0> V2;
  real<lower = 0> ka;
  real<lower = 0> sigma;
}

transformed parameters {
  vector<lower = 0>[nt] cHat;
  vector<lower = 0>[nObs] cHatObs;
  matrix<lower = 0>[nt, nCmt] x;

  x = twoCptModel(time, amt, cmt, evid,
	            	   CL, Q, V1, V2, ka);

  cHat = col(x, 2) ./ V1;

  cHatObs = cHat[iObs]; ## predictions for observed data records
}

model {
  CL ~ lognormal(log(10), 0.25);
  Q ~ lognormal(log(15), 0.5);
  V1 ~ lognormal(log(35), 0.25);
  V2 ~ lognormal(log(105), 0.5);
  sigma ~ cauchy(0, 1);

  logCObs ~ normal(log(cHatObs), sigma);
}

generated quantities{
    real cObsPred[nObs];

    for (i in 1:nObs){
      cObsPred[i] = exp(normal_rng(log(cHatObs[i]), sigma));
    }
}
