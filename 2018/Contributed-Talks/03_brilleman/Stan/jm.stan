// Note that this file contains extracts from the rstanarm package.
// Copyright (C) 2015, 2016, 2017 Trustees of Columbia University
// Copyright (C) 2016, 2017 Sam Brilleman
  
functions {

  /** 
  * Evaluate the linear predictor for the glmer submodel
  *
  * @param X Design matrix for fe
  * @param Z Design matrix for re, for a single grouping factor
  * @param Z_id Group indexing for Z
  * @param gamma The intercept parameter
  * @param beta Vector of population level parameters
  * @param bMat Matrix of group level params
  * @param shift Number of columns in bMat
  *   that correpond to group level params from prior glmer submodels
  * @return A vector containing the linear predictor for the glmer submodel
  */  
  vector evaluate_eta(matrix X, vector[] Z, int[] Z_id, real gamma, 
                      vector beta, matrix bMat, int shift) {
    int N = rows(X);    // num rows in design matrix
    int K = rows(beta); // num predictors
    int p = size(Z);    // num group level params
    vector[N] eta;
    
    if (K > 0) eta = X * beta;
    else eta = rep_vector(0.0, N);
    
    for (k in 1:p)
      for (n in 1:N)
        eta[n] = eta[n] + (bMat[Z_id[n], k + shift]) * Z[k,n];
    
    return eta;
  } 

  /** 
  * Get the indices corresponding to the lower tri of a square matrix
  *
  * @param dim The number of rows in the square matrix
  * @return A vector of indices
  */ 
  int[] lower_tri_indices(int dim) {
    int indices[dim + choose(dim, 2)];
    int mark = 1;
    for (r in 1:dim) {
      for (c in r:dim) {
        indices[mark] = (r - 1) * dim + c;
        mark = mark + 1;
      }
    }
    return indices;
  }
  
}
data {
  
  //----- Longitudinal submodels

  // number of long. submodels
  // NB this is fixed equal to 2 for this simplified jm.stan 
  // file. See the jm.stan file in rstanarm for more general 
  // code that allows for 1, 2 or 3 longitudinal submodels.
  int<lower=2,upper=2> M; 
  
  // population level dimensions
  int<lower=0> y_N[2]; // num observations
  int<lower=0> y_K[2]; // num predictors

  // population level data
  // NB these design matrices are evaluated at the observation times
  vector[y_N[1]] y1; // response vectors
  vector[y_N[2]] y2;  
  matrix[y_N[1],y_K[1]] y1_X; // fe design matrix
  matrix[y_N[2],y_K[2]] y2_X; 
  vector[y_K[1]] y1_Xbar; // predictor means
  vector[y_K[2]] y2_Xbar;
  
  // group level dimensions
  int<lower=0> b_N;     // num groups
  int<lower=0> b_K;     // total num params
  int<lower=0> b_KM[2]; // num params in each submodel

  // group level data
  vector[y_N[1]] y1_Z[b_KM[1]]; // re design matrix
  vector[y_N[2]] y2_Z[b_KM[2]];
  int<lower=0> y1_Z_id[y_N[1]]; // group ids for y*_Z
  int<lower=0> y2_Z_id[y_N[2]];
  
  //----- Event submodel
  
  // data for calculating event submodel linear predictor in GK quadrature
  // NB these design matrices are evaluated AT the event time and
  // the (unstandardised) quadrature points
  int<lower=0> e_K;           // num. of predictors in event submodel
  int<lower=0> a_K;           // num. of association parameters
  int<lower=0> Npat;          // num. individuals
  int<lower=0> Nevents;       // num. events (ie. not censored)  
  int<lower=0> qnodes;        // num. of nodes for GK quadrature 
  int<lower=0> Npat_times_qnodes; 
  int<lower=0> nrow_e_Xq;     // num. rows in event submodel predictor matrix
  vector[nrow_e_Xq] e_times;  // event times and quadrature points
  matrix[nrow_e_Xq,e_K] e_Xq; // predictor matrix (event submodel)
  vector[e_K] e_xbar;         // predictor means (event submodel)
  int<lower=0> basehaz_df;    // df for B-splines baseline hazard
  matrix[nrow_e_Xq,basehaz_df] basehaz_X; // design matrix (basis terms) for baseline hazard
  vector[Npat_times_qnodes] qwts; // GK quadrature weights with (b-a)/2 scaling 

  // data for calculating long. submodel linear predictor in GK quadrature
  // NB these design matrices are evaluated AT the event time and
  // the (unstandardised) quadrature points
  int<lower=0> nrow_y_Xq[2]; // num. rows in long. predictor matrix at quadpoints
  matrix[nrow_y_Xq[1],y_K[1]] y1_Xq; // fe design matrix at quadpoints
  matrix[nrow_y_Xq[2],y_K[2]] y2_Xq; 
  vector[nrow_y_Xq[1]] y1_Zq[b_KM[1]]; // re design matrix at quadpoints
  vector[nrow_y_Xq[2]] y2_Zq[b_KM[2]];
  int<lower=0> y1_Zq_id[nrow_y_Xq[1]]; // group indexing for re design matrix
  int<lower=0> y2_Zq_id[nrow_y_Xq[2]]; 

  //----- Hyperparameters for prior distributions
  
  // means for priors
  // coefficients
  vector[y_K[1]]       y1_prior_mean;
  vector[y_K[2]]       y2_prior_mean;
  vector[e_K]          e_prior_mean;
  vector[a_K]          a_prior_mean;
  vector[M]            y_prior_mean_for_intercept;
  vector<lower=0>[M]   y_prior_mean_for_aux;
  vector[basehaz_df]   e_prior_mean_for_aux;
  
  // scale for priors
  vector<lower=0>[y_K[1]] y1_prior_scale;
  vector<lower=0>[y_K[2]] y2_prior_scale;
  vector<lower=0>[e_K]    e_prior_scale;
  vector<lower=0>[a_K]    a_prior_scale;
  vector<lower=0>[M]      y_prior_scale_for_intercept;
  vector<lower=0>[M]      y_prior_scale_for_aux;
  vector<lower=0>[basehaz_df] e_prior_scale_for_aux;
  
  // lkj prior stuff
  vector<lower=0>[b_K] b_prior_scale;
  vector<lower=0>[b_K] b_prior_df;
  real<lower=0> b_prior_regularization;
}
transformed data {
  // indexing used to extract lower tri of RE covariance matrix
  int b_cov_idx[b_K + choose(b_K, 2)];
  if (b_K > 0) 
    b_cov_idx = lower_tri_indices(b_K);
}
parameters {
  real y1_gamma; // intercepts in long. submodels
  real y2_gamma;
  vector[y_K[1]] y1_z_beta; // primitive coefs in long. submodels
  vector[y_K[2]] y2_z_beta;
  vector[e_K] e_z_beta; // primitive coefs in event submodel (log hazard ratios)
  vector[a_K] a_z_beta; // primitive assoc params (log hazard ratios)
  real<lower=0> y1_aux_unscaled; // unscaled residual error SDs 
  real<lower=0> y2_aux_unscaled; 
  vector[basehaz_df] e_aux_unscaled; // unscaled coefs for baseline hazard      

  // group level params   
  vector<lower=0>[b_K] b_sd; // group level sds  
  matrix[b_K,b_N] z_b_mat;   // unscaled group level params 
  cholesky_factor_corr[b_K > 1 ? b_K : 0] 
    b_cholesky;              // cholesky factor of corr matrix
}
transformed parameters {
  vector[y_K[1]] y1_beta;     // population level params for long. submodels
  vector[y_K[2]] y2_beta;
  real y1_aux; // residual error SDs for long. submodels
  real y2_aux;  
  vector[e_K] e_beta;       // coefs in event submodel (log hazard ratios)
  vector[a_K] a_beta;       // assoc params in event submodel (log hazard ratios) 
  vector[basehaz_df] e_aux; // b-spline coefs for baseline hazard
  matrix[b_N,b_K] b_mat;    // group level params
  
  // coefs for long. submodels
  y1_beta = y1_z_beta .* y1_prior_scale + y1_prior_mean;
  y2_beta = y2_z_beta .* y2_prior_scale + y2_prior_mean;

  // coefs for event submodel (incl. association parameters)
  e_beta = e_z_beta .* e_prior_scale + e_prior_mean;
  a_beta = a_z_beta .* a_prior_scale + a_prior_mean;
  
  // residual error SDs for long. submodels
  y1_aux = y1_aux_unscaled * y_prior_scale_for_aux[1] + y_prior_mean_for_aux[1];
  y2_aux = y2_aux_unscaled * y_prior_scale_for_aux[2] + y_prior_mean_for_aux[2];

  // b-spline coefs for baseline hazard
  e_aux = e_aux_unscaled .* e_prior_scale_for_aux + e_prior_mean_for_aux;
 
  // group level params
  if (b_K == 1) 
    b_mat = (b_sd[1] * z_b_mat)'; 
  else if (b_K > 1) 
    b_mat = (diag_pre_multiply(b_sd, b_cholesky) * z_b_mat)';
}
model {

  //---- Log-lik for longitudinal submodels

  {
    // declare linear predictors
    vector[y_N[1]] y1_eta; 
    vector[y_N[2]] y2_eta;

    // evaluate linear predictor for each long. submodel
    y1_eta = evaluate_eta(y1_X, y1_Z, y1_Z_id, y1_gamma, y1_beta, b_mat, 0);
    y2_eta = evaluate_eta(y2_X, y2_Z, y2_Z_id, y2_gamma, y2_beta, b_mat, b_KM[1]);
    
    // increment the target with the log-lik
    target += normal_lpdf(y1 | y1_eta, y1_aux);
    target += normal_lpdf(y2 | y2_eta, y2_aux);
  }
  
  //----- Log-lik for event submodel (Gauss-Kronrod quadrature)
  
  {
    vector[nrow_y_Xq[1]] y1_eta_q; 
    vector[nrow_y_Xq[2]] y2_eta_q;
    vector[nrow_e_Xq] e_eta_q; 
    vector[nrow_e_Xq] log_basehaz;  // log baseline hazard AT event time and quadrature points
    vector[nrow_e_Xq] log_haz_q;    // log hazard AT event time and quadrature points
    vector[Nevents] log_haz_etimes; // log hazard AT the event time only
    vector[Npat_times_qnodes] log_haz_qtimes; // log hazard AT the quadrature points
    
    // Event submodel: linear predictor at event time and quadrature points
    e_eta_q = e_Xq * e_beta;
    
    // Long. submodel: linear predictor at event time and quadrature points
    y1_eta_q = evaluate_eta(y1_Xq, y1_Zq, y1_Zq_id, y1_gamma, y1_beta, b_mat, 0);
    y2_eta_q = evaluate_eta(y2_Xq, y2_Zq, y2_Zq_id, y2_gamma, y2_beta, b_mat, b_KM[1]);
  
    // Event submodel: add on contribution from association structure to
    // the linear predictor at event time and quadrature points
    e_eta_q = e_eta_q + a_beta[1] * y1_eta_q + a_beta[2] * y2_eta_q;
    
    // Log baseline hazard at event time and quadrature points
    log_basehaz = basehaz_X * e_aux; 
    
    // Log hazard at event time and quadrature points
    log_haz_q = log_basehaz + e_eta_q;
  
    // Log hazard at event times only
    log_haz_etimes = head(log_haz_q, Nevents);
  
    // Log hazard at quadrature points only
    log_haz_qtimes = tail(log_haz_q, Npat_times_qnodes);
 
    // Log likelihood for event submodel
    // NB The first term is the log hazard contribution to the log  
    // likelihood for the event submodel. The second term is the log  
    // survival contribution to the log likelihood for the event submodel.  
    // The latter is obtained by summing over the quadrature points to get 
    // the approximate integral (i.e. cumulative hazard). Note that the
    // 'qwts' vector already incorporates (b-a)/2 scaling such that the
    // integral is evaluated over limits (a,b) rather than (-1,+1), where
    // 'a' is baseline, i.e. time 0, and 'b' is the event or censoring
    // time for the individual.
    target += sum(log_haz_etimes) - dot_product(qwts, exp(log_haz_qtimes));  
  }    
    
  //----- Log-priors
    
    // intercepts for long. submodels
    target += normal_lpdf(y1_gamma | 
      y_prior_mean_for_intercept[1], y_prior_scale_for_intercept[1]);    
    target += normal_lpdf(y2_gamma | 
      y_prior_mean_for_intercept[2], y_prior_scale_for_intercept[2]);    

    // coefficients for long. submodels   
    target += normal_lpdf(y1_z_beta | 0, 1);
    target += normal_lpdf(y2_z_beta | 0, 1);
    
    // coefficients for event submodel
    target += normal_lpdf(e_z_beta | 0, 1);
    target += normal_lpdf(a_z_beta | 0, 1);
    
    // residual error SDs for long. submodels
    target += normal_lpdf(y1_aux_unscaled | 0, 1);
    target += normal_lpdf(y2_aux_unscaled | 0, 1);
    
    // b-spline coefs for baseline hazard
    target += normal_lpdf(e_aux_unscaled | 0, 1);

    // group level terms
      // sds
      target += student_t_lpdf(b_sd | b_prior_df, 0, b_prior_scale);
      // primitive coefs
      target += normal_lpdf(to_vector(z_b_mat) | 0, 1); 
      // corr matrix
      if (b_K > 1) 
        target += lkj_corr_cholesky_lpdf(b_cholesky | b_prior_regularization);
}
generated quantities {
  real y1_alpha; // transformed intercepts for long. submodels
  real y2_alpha;
  real e_alpha; // transformed intercept for event submodel 
  vector[size(b_cov_idx)] b_cov; // var-cov for REs
    
  // Transformed intercepts for long. submodels
  y1_alpha = y1_gamma - dot_product(y1_Xbar, y1_beta);
  y2_alpha = y2_gamma - dot_product(y2_Xbar, y2_beta);
  
  // Transformed intercept for event submodel 
  e_alpha = 0.0 - dot_product(e_xbar, e_beta);
  
  // Transform variance-covariance matrix for REs
  if (b_K == 1)
    b_cov[1] = b_sd[1] * b_sd[1];
  else
    b_cov = to_vector(quad_form_diag(
      multiply_lower_tri_self_transpose(b_cholesky), b_sd))[b_cov_idx];
}
