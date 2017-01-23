functions{

    real[] twoCptNeutModelODE(real t,
			real[] x,
			real[] parms,
			real[] rdummy,
			int[] idummy) {
    real k10;
    real k12;
    real k21;
    real CL;
    real Q;
    real V1;
    real V2;
    real ka;
    real mtt;
    real circ0;
    real gamma;
    real alpha;
    real ktr;
    real dxdt[8];
    real conc;
    real EDrug;
    real transit1;
    real transit2;
    real transit3;
    real circ;
    real prol;

    CL = parms[1];
    Q = parms[2];
    V1 = parms[3];
    V2 = parms[4];
    ka = parms[5];
    mtt = parms[6];	
    circ0 = parms[7];
    gamma = parms[8];
    alpha = parms[9];

    k10 = CL / V1;
    k12 = Q / V1;
    k21 = Q / V2;

    ktr = 4 / mtt;
  
    dxdt[1] = -ka * x[1];
    dxdt[2] = ka * x[1] - (k10 + k12) * x[2] + k21 * x[3];
    dxdt[3] = k12 * x[2] - k21 * x[3];
    conc = x[1]/V1;
    EDrug = alpha * conc;
    // x[4], x[5], x[6], x[7] and x[8] are differences from circ0.
    prol = x[4] + circ0;
    transit1 = x[5] + circ0;
    transit2 = x[6] + circ0;
    transit3 = x[7] + circ0;
    circ = fmax(machine_precision(), x[8] + circ0); // Device for implementing a modeled initial condition
    dxdt[4] = ktr * prol * ((1 - EDrug) * ((circ0 / circ)^gamma) - 1);
    dxdt[5] = ktr * (prol - transit1);
    dxdt[6] = ktr * (transit1 - transit2);
    dxdt[7] = ktr * (transit2 - transit3);
    dxdt[8] = ktr * (transit3 - circ);

    return dxdt;

  }

  real[] twoCptNeutModel1(real t0, real[] t, real[] init, real amt, int cmt, int evid,
			  real CL, real Q, real V1, real V2, real ka,
			  real mtt, real circ0, real gamma, real alpha,
			  real[] rdummy, int[] idummy) {
    int nCmt = 8;
    real x[nCmt];
    real parms[9];
    real temp[1, nCmt];

    parms[1]= CL;
    parms[2] = Q;
    parms[3] = V1;
    parms[4] = V2;
    parms[5] = ka;
    parms[6] = mtt;
    parms[7] = circ0;
    parms[8] = gamma;
    parms[9] = alpha;

    x = rep_array(0, nCmt);

    if (t0 == t[1]) {
      x = init;
    } else {
      temp = integrate_ode_rk45(twoCptNeutModelODE, init, t0, t, parms, rdummy, idummy,
				  1.0E-6, 1.0E-6, 1.0E8);
      x = to_array_1d(temp);
    }
   
    if (evid == 1) x[cmt] = x[cmt] + amt;

    return x;
  }

  matrix twoCptNeutModel(real[] time, real[] amt, int[] cmt, int[] evid, 
		     real CL, real Q, real V1, real V2, real ka,
		     real mtt, real circ0, real gamma, real alpha,
		     real[] rdummy, int[] idummy) {
    int nCmt = 8;
    real init[nCmt];
    matrix[size(time), nCmt] result;
    int nt;
    real t0;

    nt = size(time);

    init = rep_array(0, nCmt);
    t0 = time[1];
    for(i in 1:nt){
      init = twoCptNeutModel1(time[max(1, i - 1)], time[i:i], init, amt[i], cmt[i], evid[i],
			       CL, Q, V1, V2, ka, mtt, circ0, gamma, alpha, rdummy, idummy);
      for(j in 1:nCmt) result[i, j] = init[j];
      t0 = time[i];
    }
    return result;
  }
  
}

data {
  int<lower = 1> nt;
  int<lower = 1> nObsPK;
  int<lower = 1> nObsPD;
  int<lower = 1> iObsPK[nObsPK];
  int<lower = 1> iObsPD[nObsPD];
  real<lower = 0> amt[nt];
  int<lower = 1> cmt[nt];
  int<lower = 0> evid[nt];
  real<lower = 0> time[nt];
  vector<lower = 0>[nObsPK] cObs;
  vector<lower = 0>[nObsPD] neutObs;
  real<lower = 0> CLPrior;
  real<lower = 0> QPrior;
  real<lower = 0> V1Prior;
  real<lower = 0> V2Prior;
  real<lower = 0> kaPrior;
  real<lower = 0> CLPriorCV;
  real<lower = 0> QPriorCV;
  real<lower = 0> V1PriorCV;
  real<lower = 0> V2PriorCV;
  real<lower = 0> kaPriorCV;
  real<lower = 0> circ0Prior;
  real<lower = 0> circ0PriorCV;
  real<lower = 0> mttPrior;
  real<lower = 0> mttPriorCV;
  real<lower = 0> gammaPrior;
  real<lower = 0> gammaPriorCV;
  real<lower = 0> alphaPrior;
  real<lower = 0> alphaPriorCV;
}

transformed data {
  vector[nObsPK] logCObs;
  vector[nObsPD] logNeutObs;
  int idummy[0];
  real rdummy[0];
  int nCmt = 8;

  logCObs = log(cObs);
  logNeutObs = log(neutObs);
}

parameters {
  real<lower = 0> CL;
  real<lower = 0> Q;
  real<lower = 0> V1;
  real<lower = 0> V2;
  real<lower = 0> ka;
  real<lower = 0> mtt;
  real<lower = 0> circ0;
  real<lower = 0> alpha;
  real<lower = 0> gamma;
  real<lower = 0> sigma;
  real<lower = 0> sigmaNeut;
}

transformed parameters {
  vector[nt] cHat;
  vector[nObsPK] cHatObs;
  vector[nt] neutHat;
  vector[nObsPD] neutHatObs;

  matrix[nt, nCmt] x;

  x = twoCptNeutModel(time,
		       amt,
		       cmt,
		       evid,
		       CL, Q, V1, V2, ka,
		       mtt, circ0, gamma, alpha,
		       rdummy, idummy);
  cHat = x[, 2] / V1;
  neutHat = x[, nCmt] + circ0;

  cHatObs = cHat[iObsPK]; ## predictions for observed data records  }
  neutHatObs = neutHat[iObsPD]; ## predictions for observed data records
}

model {
  
  // CL ~ lognormal(log(CLPrior), CLPriorCV);
  // Q ~ lognormal(log(QPrior), QPriorCV);
  // V1 ~ lognormal(log(V1Prior), V1PriorCV);
  // V2 ~ lognormal(log(V2Prior), V2PriorCV);
  // ka ~ lognormal(log(kaPrior), kaPriorCV);
  // sigma ~ cauchy(0, 1);
  CL ~ lognormal(log(10), 0.25);
  Q ~ lognormal(log(15), 0.5);
  V1 ~ lognormal(log(35), 0.25);
  V2 ~ lognormal(log(105), 0.5);
  sigma ~ cauchy(0, 1);
  mtt ~ lognormal(log(mttPrior), mttPriorCV);
  circ0 ~ lognormal(log(circ0Prior), circ0PriorCV);
  alpha ~ lognormal(log(alphaPrior), alphaPriorCV);
  gamma ~ lognormal(log(gammaPrior), gammaPriorCV);
  sigmaNeut ~ cauchy(0, 1);

  logCObs ~ normal(log(cHatObs), sigma); ## observed data likelihood

  logNeutObs ~ normal(log(neutHatObs), sigmaNeut);
}

generated quantities {
  real cPred[nt];
  real neutPred[nt];

  for (i in 1:nt) {
    if(time[i] == 0) {
      cPred[i] = 0;
    } else {
      cPred[i] = exp(normal_rng(log(fmax(machine_precision(), cHat[i])), sigma));
    }
    neutPred[i] = exp(normal_rng(log(fmax(machine_precision(), neutHat[i])), sigmaNeut));
  }

}

