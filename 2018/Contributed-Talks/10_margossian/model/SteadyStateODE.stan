// Model to calculate drug concentration after patient has reached a
// steady state.

functions {
  // Function specifying the ODE for a Two Cpt model
  real[] twoCptModelODE(real t,
			real[] x,
			real[] parms,
			real[] x_r,
			int[] x_i){
    real CL = parms[1];
    real Q = parms[2];
    real V1 = parms[3];
    real V2 = parms[4];
    real ka = parms[5];
    
    real k10 = CL / V1;
    real k12 = Q / V1;
    real k21 = Q / V2;
    
    real y[3];

    y[1] = -ka*x[1];
    y[2] = ka*x[1] - (k10 + k12)*x[2] + k21*x[3];
    y[3] = k12*x[2] - k21*x[3];

    return y;
  }
  
  // Functor to pass to the algebraic solver
  vector f(vector y, vector theta, real[] x_r, int[] x_i) {
    real amt = x_r[2];
    int cmt = x_i[1];
    real y_ii[3] = to_array_1d(y);

    y_ii[cmt] = y_ii[cmt] + amt;
    y_ii = integrate_ode_rk45(twoCptModelODE, y_ii, 0, rep_array(x_r[1], 1), 
                              to_array_1d(theta), rep_array(0.0, 1),
                              rep_array(0, 1))[1];

    // return the difference between evolved and initial state
    return to_vector(y_ii) - y;
  }

  real[] twoCptModel1(real t0, real[] t, real[] init, real amt, int cmt, 
                      int evid, vector theta, real[] x_r, int[] x_i) {
    int nCmt = 3;
    real x[nCmt];
    real temp[1, nCmt];
    x = rep_array(0, nCmt);

    if (t0 == t[1]) {
      x = init;
    } else {
      temp = integrate_ode_rk45(twoCptModelODE, init, t0, t, to_array_1d(theta),
                                x_r, x_i);
      x = to_array_1d(temp);
    }

    if (evid == 1) x[cmt] = x[cmt] + amt;

    return x;
  }
  
  // Function to solve the ODEs within the event schedule of the clinical trial
  matrix twoCptModel(vector init, real[] time, real[] amt, real[] ii,
                     int[] cmt, int[] evid, vector theta, real[] x_r,
                     int[] x_i) {
    int nCmt = 3;
    real y[nCmt];
    matrix[size(time), nCmt] result;
    int nt;
    
    for (j in 1:nCmt) result[1, j] = init[j];

    nt = size(time);
    y = to_array_1d(init);
    for (i in 2:nt) {
      y = twoCptModel1(time[i - 1], time[i:i], y, amt[i], cmt[i], evid[i],
                       theta, x_r, x_i);

      for (j in 1:nCmt) result[i, j] = y[j];
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
  real x_r2[0];
  int x_i2[0];
  real rel_tol = 1e-10;  // default
  real f_tol = 1e-5;  // adjusted empirically
  int max_steps = 1000;  // default
  
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
  init = algebra_solver(f, init_guess, theta, x_r, x_i, 
                        rel_tol, f_tol, max_steps);
  init[cmt[1]] = init[cmt[1]] + amt[1];

  x = twoCptModel(init, time, amt, ii, cmt, evid, theta, x_r2, x_i2);

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
