functions {
  
  // declare C++ functions
  real getCubicHittingTime(real[] theta, real t0, real t1, real threshold, int NMAX, real TOL);
  real getMax(real[] theta, real t0, real t1, int NMAX, real TOL);
  
  // coag system
  real[] dz_dt(real t, real[] z, real[] theta, real[] x_r, int[] x_i) {
    
    real lambda = theta[1];
    real dz = -lambda*z[1];

    return { dz };
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
    for(idx in 1:N) if(x[idx] < a) return idx;
    
    return -1;
  }
  
}

data {
  
  int<lower = 0> Nt; // total num time points to capture ODE
  real ts[Nt];       // time points to capture ODE
  
  real<lower=0> teg;
}

parameters { 
  real<lower=0> y0[1];
}

transformed parameters {
  
  real<lower=0> teg_sim;
  
  {
    // temp variables to get soln curves from forward sim, which in turn are used
    // to get continuous version of Fn soln via splines, which is used to get TEG
    real z[Nt-1, 1];
    vector[Nt] Clot;
    
    vector[Nt+2*(Nt-1)+Nt] spline_params;
    vector[Nt] M;
    vector[Nt-1] C;
    vector[Nt-1] Ct;
    
    int i;
    real lambda[1] = {1.0};
  
    // integrate ODE and get Fn reprenting fibrin
    z = integrate_ode_bdf(dz_dt, y0, 0, ts[2:Nt], lambda, rep_array(0.0, 0),
                                    rep_array(0, 0), 1e-15, 1e-15, 1e3);
    
    Clot = append_row(y0[1], to_vector(z[,1]));
      
    // use ODE soln to get cubic spline coefficients
    spline_params = get_cubic_spline_coefficients(to_vector(ts),Clot);
    M = spline_params[1:Nt];
    C = spline_params[(Nt+1):(2*Nt-1)];
    Ct = spline_params[(2*Nt):(3*Nt-2)];
  
    
    i = first_interval_index(Clot, 0.6);
    teg_sim = getCubicHittingTime({M[i-1], M[i], C[i-1], Ct[i-1]}, ts[i-1], ts[i], 0.6, 100, 1e-14);
  }
}

model {
  
  y0 ~ normal(1, 0.1);
  teg ~ normal(teg_sim, 1e-2);
}
