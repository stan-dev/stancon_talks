//Copyright 2017 Aaron Goodman <aaronjg@stanford.edu>. Licensed under the GPLv3 or later.
functions{

  // helper functions for calculating the precion matrix
  real exp_logit(real x){
    return x / (1 - x);
  }
  real exp_half_logit(real x){
    return  x / (1 - square(x));
  }
  real on_diagonal_diff(real sq_exp_lambda_t){
    return(1+exp_logit(sq_exp_lambda_t));
  }
  real on_diagonal_middle_diff(real sq_exp_lambda_t1, real sq_exp_lambda_t2){
    return(1+exp_logit(square(sq_exp_lambda_t1)) + exp_logit(square( sq_exp_lambda_t2)));
  }
  real off_diagonal_diff(real exp_lambda_t){
    return(- exp_half_logit(exp_lambda_t));
  }

  // calculate the determinant of the a symmetric, tridiagonal matrix with
  // diagonal elements in m[1,1:n] and off diagonals in m[2,1:n-1]
  real tridiag_det(real [,] m){
    real f0 = 1;
    real f1 = m[1,1];
    real ftemp;
    for(i in 2:dims(m)[2]){
      ftemp = f1;
      f1 = m[1,i] .* f1 - square(m[2,i-1]) .* f0;
      f0 = ftemp;
    }
    return f1;
  }

  // calculate the precision matrix of an OU covariance kernel
  // with mean reversion parameter lambda and time intervals week_num
  real [,] ou_precision(real []week_num,real lambda){
    int n = size(week_num);
    vector [n] weeks = to_vector(week_num);
    real  ou_p [2,n];
    real min_diff = min(weeks[2:] - weeks[:n-1]);
    ou_p[1,1] = on_diagonal_diff(exp((week_num[1] - week_num[2])*lambda));
    ou_p[1,n] = on_diagonal_diff(exp((week_num[n-1] - week_num[n])*lambda));
    for(k in 1:(n-1)){
      real difference = week_num[k+1] - week_num[k];
      real exp_lambda_t = exp(lambda *(-difference));
      real off_diag = off_diagonal_diff(exp_lambda_t);
      ou_p[2,k] = off_diag;
      if(k != 1){
        real back_difference = week_num[k] - week_num[k-1];
        real exp_lambda_t0 = exp(lambda *(-back_difference));      
        real on_diag = on_diagonal_middle_diff(exp_lambda_t0,exp_lambda_t);
        ou_p[1,k] = on_diag;
      }
    }
    return(ou_p);
  }
}
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
  real mean_temporal_difference = (time[NSMPL] - time[1]) / (rows(time) - 1);
  int sample_start[replicates];
  sample_start[1] = 1;
  if(replicates > 1)
    for(i in 2:replicates)
      sample_start[i] = sample_start[i-1] + replicate_samples[i-1];
}

parameters{
  real lambda_log;
  real kappa_log;
  real mu;
  real <lower=2> student_df;
  vector[NSMPL] latent_value_raw;
}

transformed parameters{
  real sigma_log = 0.5*(kappa_log + lambda_log + log2());
  real sigma = exp(sigma_log);
  real lambda = exp(lambda_log);
  real kappa = exp(kappa_log);
  real kappa_inv = exp(-kappa_log);
  real kappa_sqrt = exp(0.5*kappa_log);
  vector[NSMPL] latent_value = mu + latent_value_raw*kappa_sqrt;
}
model {
  target += lambda_log;
  target += kappa_log;
  {
    int sample_offset = 1;
    int taxa_offset = 1;
    for(i in 1:replicates){
      int n = replicate_samples[i];
      int offset = sample_start[i];
      real next_lv;
      real prev_lv;
      vector[n] tot;
      real time_sample [n] = to_array_1d(segment(time,offset,n));
      real ou_p [2,n] = ou_precision(time_sample,lambda);
      for(k in 0:n){
        real normalized;
        real lv = next_lv;
        if(k > 0){
          normalized = lv * ou_p[1,k];
        }
        if(k > 1){
          normalized = normalized + prev_lv * ou_p[2,k-1];
        }
        if(k < n){
          prev_lv = lv;
          next_lv = latent_value_raw[offset+k];
          if(k > 0){
            normalized = normalized + next_lv * ou_p[2,k];
          }
        }
        if(k > 0){
          tot[k] = normalized .* lv;
        }
      }
      target += -0.5*(student_df + n) * log1p(sum(tot) / 
                                                  (student_df - 2));
      target +=  -n * 0.5 * log(student_df - 2) +
                  lgamma(0.5*(student_df + n)) -
                  lgamma(0.5*student_df);
      target += 0.5 * log(tridiag_det(ou_p));
    }
  }
  
  student_df ~ gamma(2,.1);  
  lambda ~ gamma(2,2);
  kappa ~ gamma(2,2);
  mu ~ normal(0,5);

  value ~ poisson_log(latent_value);
}
generated quantities{
  real carrying_capacity = exp(mu + kappa);
}
