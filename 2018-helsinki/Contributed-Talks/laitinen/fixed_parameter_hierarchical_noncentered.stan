
data{
  int <lower=0> T;               // total number of samples
  int <lower=0> n_series;        // number of series
  int <lower=0> samples_per_series[n_series]; 
  vector<lower=0> [T] observations; // observation concatenated
  vector [T] time;               // observation times concatenated
  
  vector[n_series] kappa_log;
  vector[n_series] mu;
}
transformed data{
  
  // vector[T] observation_vec = to_vector(observations);
  vector[T] time_vec = to_vector(time);
  
  int sample_start[n_series];
  sample_start[1] = 1;
  if(n_series > 1) {
    for(i in 2:n_series)
      sample_start[i] = sample_start[i-1] + samples_per_series[i-1];
  }
  

  
}

parameters{
  vector[n_series] lambda_log;
  // vector[n_series] kappa_log;
  // vector[n_series] mu;
  real <lower=2> student_df;
  vector[T] epsilon; // t-process 
}

transformed parameters{
  vector[n_series] sigma_log = 0.5*(kappa_log + lambda_log + log2());
  vector[n_series] sigma = exp(sigma_log);
  vector[n_series] lambda = exp(lambda_log);
  vector[n_series] kappa = exp(kappa_log);
  vector[n_series] kappa_inv = exp(-kappa_log);
  vector[n_series] kappa_sqrt = exp(0.5*kappa_log);
  vector[T] latent_observation;
  
  // For each raw latent observation, transform it using the
  // conditional expression for the OU process.
  for(i in 1:n_series){
    
    int n = samples_per_series[i];
    int offset = sample_start[i];
    vector [n-1] time_diff = time_vec[offset+1:offset+n-1] -time_vec[offset:offset+n-2];
    real cum_squares = 0;
    real last_t = -1;
    real lv;
    for(k in 0:n-1){
      real lv_raw = epsilon[offset + k];
      if(k == 0){
        //For the first latent observation use the stationary distribution.
        lv = mu[i] + lv_raw * kappa_sqrt[i];
      }else{
        real t = time_diff[k];
        real exp_neg_lambda_t = exp(-t*lambda[i]);
        real sd_scale = kappa_sqrt[i] .* sqrt(1-square(exp_neg_lambda_t));
        lv = mu[i] - (mu[i] - lv) .* exp_neg_lambda_t + lv_raw .* sd_scale;
        last_t = t;
      }
      latent_observation[offset+k] = lv;
    }
    
  }
}
model {
  
  target += sum(lambda_log);
  target += sum(kappa_log);
  
  
  // Increment the log probability according to the conditional expression for
  // the unit multivariate t-distribution.
  for(i in 1:n_series){
    int n = samples_per_series[i];
    int offset = sample_start[i];
    vector[n] sq_lv = square(epsilon[offset:(offset + n - 1)]);
    vector[n] cum_squares = cumulative_sum(append_row(0,sq_lv[1:n-1]));  
    
    // log1p_term = log1p(sq_lv[k] / (student_df + cum_squares[k] - 2));
    
    for(k in 1:samples_per_series[i]){
      target +=  (lgamma((student_df + k) * 0.5) - lgamma((student_df+ k - 1 )* 0.5));     
      target += -0.5 * (student_df + k) * log1p(sq_lv[k] / (student_df + cum_squares[k] - 2));
      target += -0.5 * log(student_df + cum_squares[k] - 2);
    }
    
  }
  
  // Add the log probability for the observations given the latent observations
  observations ~ normal(latent_observation, 0.01);
  
  
  // Prior probabilities (not the original ones).
  lambda ~ normal(0,5);
  kappa  ~ normal(0,5);
  mu     ~ normal(0,5);
  student_df ~ gamma(2,.1);
}
