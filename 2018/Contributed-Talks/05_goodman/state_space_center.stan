//Copyright 2017 Aaron Goodman <aaronjg@stanford.edu>. Licensed under the GPLv3 or later.
data{
  int <lower=0> NSMPL;
  int <lower=0> replicates;
  int <lower=0> replicate_samples[replicates];
  int <lower=0> value[NSMPL];
  vector [NSMPL] time;
}
transformed data{
  vector[NSMPL] value_vec = to_vector(value);
  vector[NSMPL] time_vec = to_vector(time);
  int sample_start[replicates];
  vector [NSMPL - replicates] time_diffs;
  
  sample_start[1] = 1;
  if(replicates > 1)
    for(i in 2:replicates)
      sample_start[i] = sample_start[i-1] + replicate_samples[i-1];

  for(i in 1:replicates){
    int smpl = replicate_samples[i];
    int start = sample_start[i];
    time_diffs[(start - (i-1)):(start - (i - 1) + smpl-2)]  = segment(time_vec, start + 1, smpl - 1) -
      segment(time_vec, start, smpl - 1);
  }
}

parameters{
  real lambda_log;
  real kappa_log;
  real mu;
  real <lower=2> student_df;
  vector[NSMPL] latent_value;
}

transformed parameters{
  real sigma_log = 0.5*(kappa_log + lambda_log + log2());
  real sigma = exp(sigma_log);
  real lambda = exp(lambda_log);
  real kappa = exp(kappa_log);
  real kappa_inv = exp(-kappa_log);
  real kappa_sqrt = exp(0.5*kappa_log);
  vector[NSMPL] latent_value_raw;
  
  for(i in 1:replicates){
    int n = replicate_samples[i];
    int offset = sample_start[i];
    vector [n-1] time_diff = segment(time_diffs, offset - (i-1), n - 1);
    real cum_squares = 0;
    real last_t = -1;
    real lv_raw;
    for(k in 0:n-1){
      real lv = latent_value[offset + k];
      if(k == 0){
        lv_raw = (lv - mu)/kappa_sqrt;
      }else{
        real t = time_diff[k];
        real exp_neg_lambda_t = exp(-t*lambda);
        real sd_scale = kappa_sqrt .* sqrt(1-square(exp_neg_lambda_t));
        lv_raw = (lv - mu - (latent_value[offset + k - 1] - mu) * exp_neg_lambda_t)/sd_scale;
      }
      latent_value_raw[offset+k] = lv_raw;
    }
  }
}
model {
  target += lambda_log;
  target += kappa_log;
  student_df ~ gamma(2,.1);
  for(i in 1:replicates){
    int n = replicate_samples[i];
    int offset = sample_start[i];
    vector [n] sq_lv = square(segment(latent_value_raw, offset, n));
    vector [n] cum_squares = cumulative_sum(append_row(0, sq_lv[1:n-1]));
    for(k in 1:NSMPL){
      target +=  (lgamma((student_df + k) * 0.5) - lgamma((student_df+ k - 1 )* 0.5));     
      target += -0.5 * (student_df + k) * log1p(sq_lv[k] / (student_df + cum_squares[k] - 2));
      target += -0.5 * log(student_df + cum_squares[k] - 2);
    }
  }

  //jacobian of the transformation
  target += -NSMPL * (0.5 * kappa_log);
  target += -0.5 * log1m_exp(-2*lambda*time_diffs);
  
  value ~ poisson_log(latent_value);
  lambda ~ gamma(2,2);
  kappa ~ gamma(2,2);
  mu ~ normal(0,5);
}

generated quantities{
  real carrying_capacity = mu + kappa;
}
