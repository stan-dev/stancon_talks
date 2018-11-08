functions {
  
  // declare C++ functions
  real getCubicHittingTime(real[] theta, real t0, real t1, real threshold, int NMAX, real TOL);
  real getMax(real[] theta, real t0, real t1, int NMAX, real TOL);
  
  // coag system
  real[] dz_dt(real t, real[] z, real[] theta, real[] x_r, int[] x_i) {
    
    // state values
    real FII = z[1];
    real FIIa = z[2];
    real AT = z[3];
    real Fg = z[4];
    real Fn = z[5];
    real tPA = z[6];
    real PAI = z[7];

    // rate constants
    real g_FII = x_r[1];
    real g_AT = x_r[2];
    
    real k_FIIa = x_r[3];
    real c_FIIa = x_r[4];
    real K_FIIa = x_r[5];
    
    real k_AT = x_r[6];
    
    real k_clot = x_r[7];
    real c_clot = x_r[8];
    real K_clot = x_r[9];
    
    real k_lys = x_r[10];
    real c_lys = x_r[11];
    real K_lys = x_r[12];
    
    real k_PAI = x_r[13];
    
    real c_tfpi = x_r[14];
    real b_tfpi = x_r[15];
    real tfpi_min = x_r[16];
    
    // patient parameters
    real c_cascade = theta[1];
    real b_cascade = theta[2];
    
    real cascade_delay = 1/(1+exp(c_cascade-b_cascade*t));
    real tfpi = (1-tfpi_min)/(1+exp(-c_tfpi+b_tfpi*t)) +  tfpi_min;

    // reaction rates
    real r_FIIa = k_FIIa*cascade_delay*tfpi*((FII/g_FII)^c_FIIa)/(K_FIIa + (FII/g_FII)^c_FIIa);
    real r_AT = k_AT*FIIa*(AT/g_AT);
    real r_clot = k_clot*FIIa*((Fg/9e-6)^c_clot)/(K_clot + (Fg/9e-6)^c_clot);
    real r_lys = k_lys*tPA*(Fn^c_lys)/(K_lys + Fn^c_lys);
    real r_PAI = k_PAI*tPA*PAI;

    // system
    real dFII = g_FII*(-r_FIIa);
    real dFIIa = r_FIIa -r_AT;
    real dAT = -g_AT*r_AT;
    real dFg = -r_clot;
    real dFn = 1e7*(r_clot - r_lys);
    real dtPA = -r_PAI;
    real dPAI = -r_PAI;

    return { dFII, dFIIa, dAT, dFg, dFn, dtPA, dPAI };
  }
  
   // solves tridiagonal matrix Ax=f with 2 across
  // across the diagonal, a on the lower diag,
  // and c on the upper diag
  vector thomas_solve(vector a, vector c, vector f) {
    
    int N = rows(f);
    
    real alpha[N];
    real beta[N-1];
    real y[N];
    vector[N] x;
    
    alpha[1] = 2;
    y[1] = f[1];
    
    for(i in 2:N) {
      beta[i-1] = a[i-1]/alpha[i-1];
      alpha[i] = 2 - beta[i-1]*c[i-1];
  
      y[i] = f[i] - beta[i-1]*y[i-1];
    }
  
    x[N] = y[N]/alpha[N];
    for(i in 1:(N-1)) {
      x[N-i] = (y[N-i]-c[N-i]*x[N-i+1])/alpha[N-i];
    }
    
    return x;
  } 
  
  vector diff(vector x) {
    
    int N = rows(x);
    vector[N-1] d = x[2:N]-x[1:(N-1)];
    return d;
    
  }
  
  // takes in points x, and function evaluated at those points, f
  vector get_cubic_spline_coefficients(vector x, vector f) {
    
    int N = rows(x)-1;
    vector[N] h = diff(x);
  
    vector[N-1] mu_ = h[1:(N-1)] ./ (h[1:(N-1)] + h[2:N]);
    vector[N-1] lambda_ = h[2:N] ./ (h[1:(N-1)] + h[2:N]);
    vector[N-1] d_ = 6.0 ./ (h[1:(N-1)]+h[2:N]) .* ((f[3:(N+1)]-f[2:N])./h[2:N]-(f[2:N]-f[1:(N-1)])./h[1:(N-1)]);
  
    // extend vectors because there are still unknowns. see pg. 358 of Quarteroni Numerical Analysis book
    vector[N] mu = append_row(mu_, 1);
    vector[N] lambda = append_row(1,lambda_);
    vector[N+1] d = append_row(append_row(d_[1], d_), d_[N-1]);
    
    // solve for coefficients
    vector[N+1] M = thomas_solve(mu, lambda, d);
    vector[N] C = (f[2:(N+1)]-f[1:N])./ h - (h/6).*(M[2:(N+1)]-M[1:N]);
    vector[N] Ct = f[1:N] - M[1:N].*((h .* h)/6);
    
    // package up cubic spline parameters M, C, Ct into single vector
    return(append_row(append_row(append_row(M,C),Ct),x));
    
  }
  
  // return the index of the array that is just after
  // the first element of the array that is greater than a
  int first_interval_index(vector x, real a) {
    
    int N = rows(x);
    for(idx in 1:N) if(x[idx] > a) return idx;
    
    return -1;
  }
  
}

data {
  
  int<lower = 0> Nt; // total num time points to capture ODE
  real ts[Nt];       // time points to capture ODE
  
  // protein and TEG data
  real proteins[5]; // FII, AT, Fg, tPA, PAI
  real teg[4];      // K, R, MA, LY30

  // vector indicating which of the protein measurements are missing
  int<lower = 0> num_missing[5];
}

transformed data {

  real g_FII = 100/1.4e-6;
  real g_AT = 100/3.4e-6;
  //real c_cascade = theta[3];
  //real b_cascade = theta[4];
  //real c_tfpi = theta[5];
  //real b_tfpi = theta[6];

  real k_FIIa = 3.5e-9;
  real c_FIIa = 1.0;
  real K_FIIa = 1.4e-6;
    
  real k_AT = 1.6e4;
    
  real k_clot = 2.900221;
  real c_clot = 1.0;
  real K_clot = 0.7527255;
    
  real k_lys = 1.0;
  real c_lys = 1.0;
  real K_lys = 0.5;
    
  real k_PAI = 4.5e5;
  
  real c_tfpi = 9.931251;
  real b_tfpi = 0.01513794;
  real tfpi_min = 0.2;
  
  real rate_constants[16] = {g_FII, g_AT, k_FIIa, c_FIIa, K_FIIa, k_AT, k_clot, c_clot, K_clot, k_lys, c_lys, K_lys,k_PAI, c_tfpi, b_tfpi, tfpi_min};
}

parameters { 

  // model missing observations as random variables (unknown parameters)
  real<lower=0> FII_missing[num_missing[1]];
  real<lower=0> AT_missing[num_missing[2]];
  real<lower=0> Fg_missing[num_missing[3]];
  real<lower=0> tPA_missing[num_missing[4]];
  real<lower=0> PAI_missing[num_missing[5]];
  
  // individual's parameters governing the onset and slow down of their clotting cascade 
  real<lower=0> theta[2];
}

transformed parameters {
  
  // forward simulated TEG values for every patient
  real teg_sim[4];
  

  // get TEG readings by forward simulating ODE
  {

    // z0, z hold soln curves for the patients clotting system
    // Clot holds the value of the TEG curve and is a function of the state Fn
    real z0[7];
    real z[Nt-1, 7];
    vector[Nt] Fn;
    vector[Nt] FnSq;
    vector[Nt] Clot;
    
    // spline parameters are used to interpolate the discrete ODE soln
    vector[Nt+2*(Nt-1)+Nt] spline_params;
    vector[Nt] M;
    vector[Nt-1] C;
    vector[Nt-1] Ct;
    int i;

    // coalesce all initial conditions in to a single data structure whether
    // missing or observed if they're missing fill with parameter, otherwise w/ data
    // dFII, dFIIa, dAT, dFg, dFn, dtPA, dPAI
    z0[1]  = (num_missing[1] != 0) ? FII_missing[1] : proteins[1];
    z0[2]  = 0.0;
    z0[3]  = (num_missing[2] != 0) ? AT_missing[1] : proteins[2];
    z0[4]  = (num_missing[3] != 0) ? Fg_missing[1] : proteins[3];
    z0[5]  = 0.0;
    z0[6]  = (num_missing[4] != 0) ? tPA_missing[1] : proteins[4];
    z0[7]  = (num_missing[5] != 0) ? PAI_missing[1] : proteins[5];

    // integrate ODE and get Fn reprenting fibrin
    z = integrate_ode_bdf(dz_dt, z0, 0, ts[2:Nt], theta, rate_constants,
                                  rep_array(0, 0), 1e-6, 1e-5, 1e3);
                                  
    Fn = append_row(z0[5], to_vector(z[,5]));
    
    // convert fibrin to "clot"
    FnSq = Fn .* Fn;
    Clot = 64.0*FnSq ./ (100.0 + FnSq);
    
    // use ODE soln to get cubic spline coefficients
    spline_params = get_cubic_spline_coefficients(to_vector(ts),Clot);
    M = spline_params[1:Nt];
    C = spline_params[(Nt+1):(2*Nt-1)];
    Ct = spline_params[(2*Nt):(3*Nt-2)];
    
    // to get R and K time first find which interval they're in, then get hitting
    // time in seconds. then convert to minutes
    i = first_interval_index(Clot, 2.0);
    
    teg_sim[1] = getCubicHittingTime({M[i-1], M[i], C[i-1], Ct[i-1]}, ts[i-1], ts[i], 2.0, 100, 1e-6);
    i = first_interval_index(Clot, 20.0);
    
    teg_sim[2] = getCubicHittingTime({M[i-1], M[i], C[i-1], Ct[i-1]}, ts[i-1], ts[i], 20.0, 100, 1e-6);
      
    teg_sim[2] = (teg_sim[2]-teg_sim[1])/60.0;
    teg_sim[1] = teg_sim[1]/60.0;
      
    // for MA use different function to find max. first find where the discrete soln
    // is at its maximum. then the continuous max will either be to the left or
    // right depending on if the deriv. there is positive or negative
    i = sort_indices_desc(Clot)[1];
    
    if(i == Nt) {
      teg_sim[3] = Clot[Nt];
    }
    else {
      if(-(M[i]/(2*(ts[i+1]-ts[i])))*pow(ts[i+1]-ts[i], 2) + C[i] > 0) {
        i = i + 1;
      }
      teg_sim[3] = getMax({M[i-1], M[i], C[i-1], Ct[i-1]}, ts[i-1], ts[i], 100, 1e-6);
    }

    teg_sim[4] = 100*(1-Clot[Nt]/teg_sim[3]);

  }

}

model {
  
  // priors for protein concentrations
  FII_missing ~ normal(50, 30);
  AT_missing ~ normal(80, 30);
  Fg_missing ~ normal(200, 50);
  tPA_missing ~ exponential(1/(4e-10));
  PAI_missing ~ exponential(1/(9.304e-10));
  
  // priors for individual patient clotting parameters
  theta[1] ~ normal(3.0, 1.0); // c_cascade
  theta[2] ~ normal(0.03, 2e-2); // b_cascade
  
  // likelihood of forward simulated TEG values
  teg[1] ~ normal(teg_sim[1], 0.1);//R
  teg[2] ~ normal(teg_sim[2], 0.2);//K
  teg[3] ~ normal(teg_sim[3], 5.0);//MA
  teg[4] ~ normal(teg_sim[4], 0.2);//Ly30
}
