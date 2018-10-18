
data {
  int num_time;
  int num_measurements;
  int num_regulators;
  int num_targets;
  int<lower=0,upper=1> regulators_measured;

  int<lower=1,upper=num_time> measurement_times[num_measurements];

  matrix<lower=0>[regulators_measured ? num_measurements : 0, num_regulators] regulator_expression;
  matrix<lower=0>[num_measurements, num_targets] expression;
  int<lower=-1,upper=1> regulation_signs[num_regulators, num_targets];

  int num_spline_basis;
  matrix[num_time, num_spline_basis] spline_basis;  // matrix of B-splines


  int<lower=0,upper=1> coeffs_prior_given;
  vector[num_spline_basis] coeffs_prior_mean[coeffs_prior_given ? num_regulators : 0];
  cov_matrix[num_spline_basis] coeffs_prior_cov[coeffs_prior_given ? num_regulators : 0];

  real<lower = 0> initial_condition_prior_sigma;
  real<lower = 0> asymptotic_normalized_state_prior_sigma;
  real<lower = 0> mean_regulatory_input_prior_sigma;
  real<lower = 0> sd_regulatory_input_prior_sigma;
  real degradation_prior_mean;
  real<lower = 0> degradation_prior_sigma;

  //Intercept is only important when regulator is measured, otherwise, it is hidden in b
  real<lower = 0> intercept_prior_sigma[regulators_measured ? 1 : 0];


  int<lower=0,upper=1> sigma_given;
  real<lower = 0> sigma_relative_prior_mean[sigma_given ? 0 : 1];
  real<lower = 0> sigma_absolute_prior_mean[sigma_given ? 0 : 1];
  real<lower = 0> sigma_relative_prior_sigma[sigma_given ? 0 : 1];
  real<lower = 0> sigma_absolute_prior_sigma[sigma_given ? 0 : 1];

  real<lower = 0> sigma_relative_data[sigma_given ? 1 : 0];
  real<lower = 0> sigma_absolute_data[sigma_given ? 1 : 0];


  real<lower=0> scale;
}

transformed data {
  matrix[num_time, num_spline_basis] spline_basis_scaled;  // matrix of B-splines
  real time_real[num_time];
  vector[num_targets] regulation_signs_real[num_regulators];
  vector[num_targets] max_target_expression;
  vector<lower=0>[num_targets] basal_transcription;

  cholesky_factor_cov[num_spline_basis] coeffs_prior_chol_cov[coeffs_prior_given ? num_regulators : 0];



  if(regulators_measured != 0 && coeffs_prior_given != 0) {
    reject("You are not supposed to give coeffs prior if regulator is measured. Use scale parameter instead.");
  }

  if (coeffs_prior_given) {
    for(r in 1:num_regulators) {
      coeffs_prior_chol_cov[r] = cholesky_decompose(coeffs_prior_cov[r]);
    }
  }

  for(i in 1:num_time) {
    time_real[i] = i;
  }

  for(r in 1:num_regulators) {
    for(t in 1:num_targets) {
      regulation_signs_real[r, t] = regulation_signs[r, t];
    }
  }

  for(t in 1:num_targets) {
    max_target_expression[t] = max(expression[,t]);
  }

  {
    for(m in 2:num_measurements) {
      if(measurement_times[m - 1] >= measurement_times[m]) {
        reject("Measurement times have to be increasing");
      }
    }
  }


  spline_basis_scaled = spline_basis * scale;

  for(t in 1:num_targets) {
    basal_transcription[t] = 0;
  }
}


parameters {
  matrix[num_spline_basis, num_regulators] coeffs;
  vector<lower=0>[num_targets] initial_condition;
  vector<lower=0>[num_targets] asymptotic_normalized_state;
  simplex[num_regulators] relative_weights[num_targets];
  vector[num_targets] mean_regulatory_input_raw;
  vector<lower=0>[num_targets] sd_regulatory_input_raw;
  vector<lower=0>[num_targets] degradation;
  real<lower = 0> intercept_raw[regulators_measured ? num_regulators : 0];

  real<lower=0> sigma_relative_param[sigma_given ? 0 : 1];
  real<lower=0> sigma_absolute_param[sigma_given ? 0 : 1];
}


transformed parameters {
  vector[num_targets] mean_regulatory_input;
  vector<lower=0>[num_targets] sd_regulatory_input;
  vector<lower=0>[num_targets] mean_synthesis;
  matrix[num_time, num_regulators] regulator_profile;
  matrix<lower=0>[num_time, num_targets] predicted_expression;
  matrix<lower=0>[num_time, num_regulators] predicted_regulator_expression;
  matrix[num_targets, num_regulators] w;
  vector[num_targets] b_centered;
  vector<lower=0>[num_targets] sensitivity;
  vector[num_regulators] intercept;

  real sigma_relative = sigma_given ? sigma_relative_data[1] : sigma_relative_param[1];
  real sigma_absolute = sigma_given ? sigma_absolute_data[1] : sigma_absolute_param[1];


  for(r in 1:num_regulators) {
    if(coeffs_prior_given == 0) {
      regulator_profile[,r] = spline_basis_scaled * coeffs[,r]; //Note the transpose
    } else {
      vector[num_spline_basis] coeffs_transformed = coeffs_prior_chol_cov[r] * coeffs[,r]  + coeffs_prior_mean[r];

      regulator_profile[,r] = spline_basis_scaled * coeffs_transformed;
    }

    {
      real min_intercept = -min(regulator_profile[,r]);
      if(regulators_measured == 0) {
        intercept[r] = min_intercept;
      } else {
        intercept[r] = min_intercept + intercept_raw[r] * intercept_prior_sigma[1];
      }
      predicted_regulator_expression[,r] = regulator_profile[,r] + intercept[r];
    }
  }


  if(num_targets > 0) {
    mean_synthesis = asymptotic_normalized_state  .* degradation .* max_target_expression;
    mean_regulatory_input = mean_regulatory_input_raw * mean_regulatory_input_prior_sigma;
    sd_regulatory_input = sd_regulatory_input_raw * sd_regulatory_input_prior_sigma;
    {
      vector[num_regulators] mean_regulator_profile;
      for(r in 1:num_regulators)
      {
        real sd_regulator_profile = sd(regulator_profile[,r]);
        mean_regulator_profile[r] = mean(regulator_profile[,r]);
        w[,r] = regulation_signs_real[r,] .* (sd_regulatory_input ./ sd_regulator_profile) .* to_vector(relative_weights[,r]); //Note: proper handling would be to take into account the covariance of the regulatory profiles, but this has issues of its own (especially in that there might be multiple solutions)
      }
      b_centered = mean_regulatory_input - w * mean_regulator_profile;
    }
  }

  for(t in 1:num_targets)
  {
    real basal_over_degradation = basal_transcription[t] / degradation[t];
    real degradation_per_step = exp(-degradation[t]);

    real residual;
    vector[num_time] synthesis;

    synthesis = inv(1 + exp( (-regulator_profile * (w[t,]')) - b_centered[t]));
    sensitivity[t] = mean_synthesis[t] / mean(synthesis);

    predicted_expression[1,t] = initial_condition[t];

    //Calculating the integral by trapezoid rule in a single pass for all values
    residual = -0.5 * synthesis[1];
    for (measurement in 2:num_measurements)
    {
      for (time in (measurement_times[measurement - 1] + 1):measurement_times[measurement])
      {
        residual = (residual + synthesis[time - 1]) * degradation_per_step;

        { //new block to let me define new vars
          real integral_value = residual + 0.5 * synthesis[time];
          predicted_expression[time,t] = basal_over_degradation + (initial_condition[t] - basal_over_degradation) * exp(-degradation[t] * (time - 1)) + sensitivity[t] * integral_value;
        }

      }
    }
  }
}

model {

    //Observation model
    for(m in 1:num_measurements) {
      row_vector[num_targets] sigma = sigma_absolute + sigma_relative * predicted_expression[measurement_times[m],];
      for(t in 1:num_targets) {
        expression[m, t] ~ normal(predicted_expression[measurement_times[m],t], sigma[t]) T[0,];
      }

      if(regulators_measured != 0) {
        row_vector[num_regulators] sigma_regulator = sigma_absolute + sigma_relative * predicted_regulator_expression[measurement_times[m],];
        for(r in 1:num_regulators) {
          regulator_expression[m, r] ~ normal(predicted_regulator_expression[measurement_times[m], r], sigma_regulator[r])  T[0,];
        }
      }
    }

    initial_condition ~ normal(0, initial_condition_prior_sigma);
    to_vector(coeffs) ~ normal(0,1); //coeffs are rescaled by the scale parameter
    intercept_raw ~ normal(0, 1);
    asymptotic_normalized_state ~ normal(1, asymptotic_normalized_state_prior_sigma);

    degradation ~ lognormal(degradation_prior_mean, degradation_prior_sigma);
    mean_regulatory_input_raw ~ normal(0, 1);
    sd_regulatory_input_raw ~ normal(0, 1);

    if(sigma_given == 0) {
      sigma_relative_param[1] ~ normal(sigma_relative_prior_mean[1], sigma_relative_prior_sigma[1]);
      sigma_absolute_param[1] ~ normal(sigma_absolute_prior_mean[1], sigma_absolute_prior_sigma[1]);
    }
}

generated quantities {
  matrix<lower=0>[num_measurements, num_targets] expression_replicates;
  matrix[num_measurements, num_targets] log_likelihood;
  vector[num_targets] b;

  if(num_targets > 0) {
    b = b_centered + (w * intercept);

    for(m in 1:num_measurements) {
      row_vector[num_targets] sigma = sigma_absolute + sigma_relative * predicted_expression[measurement_times[m],];
      for(t in 1:num_targets) {
        //Draw from the truncated normal
        real lower_bound = 0;
        real p = normal_cdf(lower_bound, predicted_expression[measurement_times[m],t], sigma[t]);
        real u = uniform_rng(p, 1);
        expression_replicates[m, t] = inv_Phi(u) * sigma[t] + predicted_expression[measurement_times[m],t];

        log_likelihood[m, t] = normal_lpdf(expression[m, t]| predicted_expression[measurement_times[m], t], sigma[t]) - log_diff_exp(1, normal_lcdf(lower_bound | predicted_expression[measurement_times[m],t], sigma[t]));
      }
    }
  }
}
