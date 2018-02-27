// Model to calculate drug concentration after patient has reached a
// steady state.

functions {
  // Analytical solution to the ODEs
  vector twoCptModel1(real dt, vector init, vector theta,
                      real amt, int cmt, int evid) {
    real CL = theta[1];
    real Q = theta[2];
    real V1 = theta[3];
    real V2 = theta[4];
    real ka = theta[5];
    real k10 = CL / V1;
    real k12 = Q / V1;
    real k21 = Q / V2;
    real ksum = k10 + k12 + k21;
    vector[3] alpha;
    vector[3] a;
    vector[3] x = rep_vector(0.0, 3);

    alpha[1] = (ksum + sqrt(ksum * ksum - 4.0 * k10 * k21))/2.0;
    alpha[2] = (ksum - sqrt(ksum * ksum - 4.0 * k10 * k21))/2.0;
    alpha[3] = ka;

    if(init[1] != 0.0){
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

    if(init[2] != 0){
      a[1] = (k21 - alpha[1]) / (alpha[2] - alpha[1]);
      a[2] = (k21 - alpha[2]) / (alpha[1] - alpha[2]);
      x[2] = x[2] + init[2] * sum(segment(a, 1, 2) .* exp(-segment(alpha, 1, 2) * dt));
      a[1] = k12 / (alpha[2] - alpha[1]);
      a[2] = -a[1];
      x[3] = x[3] + init[2] * sum(segment(a, 1, 2) .* exp(-segment(alpha, 1, 2) * dt));
    }

    if(init[3] != 0){
      a[1] = k21 / (alpha[2] - alpha[1]);
      a[2] = -a[1];
      x[2] = x[2] + init[3] * sum(segment(a, 1, 2) .* exp(-segment(alpha, 1, 2) * dt));
      a[1] = (k10 + k12 - alpha[1]) / (alpha[2] - alpha[1]);
      a[2] = (k10 + k12 - alpha[2]) / (alpha[1] - alpha[2]);
      x[3] = x[3] + init[3] * sum(segment(a, 1, 2) .* exp(-segment(alpha, 1, 2) * dt));
    }

    if(evid == 1) x[cmt] = x[cmt] + amt;

    return x;
  }
  
  // Functor to pass to the algebraic solver
  vector f(vector y, vector theta, real[] x_r, int[] x_i) {
    real ii = x_r[1];
    real amt = x_r[2];
    int cmt = x_i[1];
    int evid = 1;

    // return the difference between evolved and initial state
    return twoCptModel1(ii, y, theta, amt, cmt, evid) - y;
  }
  
  // Function to solve the ODEs within the event schedule of the clinical trial
  matrix twoCptModel(vector init, real[] time, real[] amt, real[] ii,
                     int[] cmt, int[] evid, vector theta) {
    int nCmt = 3;
    vector[nCmt] y;
    real dt;
    real t0;
    matrix[size(time), nCmt] result;
    int nt;
    
    for (j in 1:nCmt) result[1, j] = init[j];

    nt = size(time);
    t0 = time[1];
    y = init;
    for (i in 2:nt) {
      dt = time[i] - t0;
      y = twoCptModel1(dt, y, theta, amt[i], cmt[i], evid[i]);

      for (j in 1:nCmt) result[i, j] = y[j];
      t0 = time[i];
    }
    return result;
  }
}

data {
  int<lower = 1> nt;
  int<lower = 1> nObs;
  int<lower = 1> iObs[nObs];
  vector<lower = 0>[nObs] cObs;
  real<lower = 0> time[nt];
  real<lower = 0> amt[nt];
  real<lower = 0> ii[nt];
  int<lower = 1> cmt[nt];
  int<lower = 0> evid[nt];
}

transformed data {
  vector[nObs] logCObs = log(cObs);
  int nCmt = 3;
  int nParms = 5;
  vector[nCmt] init_guess;
  real x_r[2] = {ii[1], amt[1]};
  int x_i[1] = {cmt[1]};
  
  // guess for steady state. Use amt as a scale.
  init_guess[1] = amt[1];
  init_guess[2] = amt[1] * 0.5;
  init_guess[3] = amt[1] * 0.35;
}

parameters {
  real<lower = 0> CL;
  real<lower = 0> Q;
  real<lower = 0> VC;
  real<lower = 0> VP;
  real<lower = 0> ka;
  real<lower = 0> sigma;
}

transformed parameters {
  vector[nParms] theta;
  vector<lower = 0>[nt] cHat;
  vector<lower = 0>[nObs] cHatObs;
  matrix<lower = 0>[nt, nCmt] x;
  vector[nCmt] init;

  theta[1] = CL;
  theta[2] = Q;
  theta[3] = VC;
  theta[4] = VP;
  theta[5] = ka;

  // start the system at a steady state
  init = algebra_solver(f, init_guess, theta, x_r, x_i);

  x = twoCptModel(init, time, amt, ii, cmt, evid, theta);

  cHat = col(x, 2) ./ VC;

  cHatObs = cHat[iObs];  // predictions for observed data records
}

model {
  CL ~ lognormal(log(10), 0.5);
  Q ~ lognormal(log(15), 0.5);
  VC ~ lognormal(log(35), 0.5);
  VP ~ lognormal(log(105), 0.5);
  sigma ~ cauchy(0, 1);

  logCObs ~ normal(log(cHatObs), sigma);
}

generated quantities{
    real cObsPred[nObs];
    for (i in 1:nObs){
      cObsPred[i] = exp(normal_rng(log(cHatObs[i]), sigma));
    }
}
