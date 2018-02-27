// our Stan model, saved as vsb.stan
// first we define the function that takes data and parameters and returns predicted market shares
functions {
  // calculates shares for a given market
  row_vector shares(real alpha, vector beta, matrix bigX, matrix Sigma, row_vector xi, matrix z) {
    matrix[rows(z), cols(xi)+1] utilities;
    matrix[rows(z), cols(xi)+1] probs;
    row_vector[cols(xi)+1] shares;
    // 1. Rather than pass in p and x separately, we'll pass in bigX = append_col(p, X)
    // 2. append alpha_shock, beta_shock
    {
      matrix[rows(z), cols(xi)] tmp;
      
      tmp = rep_matrix((bigX*append_row(alpha, beta) + xi')', rows(z));
      
      // replace the append_col wing single values (might speed things up)
      utilities[1:rows(z), 1:cols(xi)] = tmp + z * cholesky_decompose(Sigma)' * bigX';
      utilities[1:rows(z),cols(utilities)] = rep_vector(0.0, rows(z));
      
      for(i in 1:rows(z)) {
         probs[i] = softmax(utilities[i]')';
      }
      
    }
    
    for(j in 1:cols(probs)) {
      shares[j] = mean(col(probs, j));
    }
    
    return(shares);
  }
}
// next define our data
data {
  int NS; // number of individuals in integration
  int J; // number of products
  int T; // number of markets
  int P; // number of features
  matrix[NS, P+1] z; // normal(0,1) draws of the shocks
  matrix[T, J] price; // price for each unit
  int sales[T, J+1]; // unit sales across T markets for J+1 products (inc outside good)
  matrix[T*J, P] X_repeat; // T Xs stacked on top of each other. This format allows characteristics to vary across markets.
  real nu;
}
// next join the product data together into single matrices
transformed data {
  matrix[T*J, P+1] bigX;
  bigX = append_col(to_vector(price'), X_repeat);
}
// define parameters
parameters {
  real alpha; 
  vector[P] beta;
  vector[P] gamma;
  real gamma0;
  real<lower = 0> price_scale;
  matrix[T, J] xi;
  vector<lower = 0>[P+1] scales;
  corr_matrix[P+1] Omega;
  real<lower = 0> lambda;
}

transformed parameters {
  cov_matrix[P+1] Sigma;
  Sigma = quad_form_diag(Omega, scales);
}
// and the model
model {
  // priors
  alpha ~ normal(0, 1);
  beta ~ normal(0, 1);
  gamma0 ~ normal(0, 3);
  gamma ~ normal(0, 3);
  price_scale ~ normal(0, 1);
  lambda ~ normal(0, 1);
  to_vector(xi) ~ normal(0, 1);
  scales ~ normal(0, 1);
  Omega ~ lkj_corr(nu);
  
  // model of price -- this helps pin down xi
  for (i in 1:rows(to_vector(price'))) {
    log(to_vector(price')[i]) ~ normal(gamma0 + X_repeat[i] * gamma + lambda * to_vector(xi')[i], price_scale);
  }
  
  
  // model of sales 
  {
    matrix[T, J+1] pred_shares;
    for(t in 1:T) {
      // predicted market shares given data and parameters
      pred_shares[t] = shares(alpha, beta, bigX[(t*J-J+1):(t*J)], Sigma, xi[t], z);
      // sales are measured with multinomial measurement error
      sales[t] ~ multinomial(pred_shares[t]');
    }
    
  }
}