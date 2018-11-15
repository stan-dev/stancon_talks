target_ODE <- function(t, state, parameters, protein_transform = identity)
{
  with(as.list(c(state, parameters)), {
    regulatoryInput = bias + weight * protein_transform(protein(t));
    dX = basal_transcription + (sensitivity/(1 + exp(-regulatoryInput))) - degradation * x;
    list(dX)
  })
}

constant_synthesis_ODE <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dX = synthesis - degradation * x;
    list(dX)
  })
}

gp_covariance <- function(distance, gp_scale, gp_length, periodic = FALSE, period = 1) {
  if(periodic) {
    return((gp_scale^2) * exp(-2*(sin(3.141592 * abs(distance) / period)^2) / (gp_length ^ 2)))
  } else {
    return((gp_scale^2) * exp( (-0.5 / (gp_length ^ 2)) * (distance ^ 2) ));
  }

}


generate_random_profile <- function(time, scale, length, mean_func = 0, periodic = FALSE, period = 1, positive_transform = TRUE) {
  # Construct the squared exponential covariance matrix
  cov_matrix = array(0,c(length(time), length(time)));
  maxTime = max(time)
  minTime = min(time)
  for(i in 1:length(time))  {
    cov_matrix[i,i] = scale^2 + 0.00000001; #Adding jitter to ensure positive semi-definiteness
    if (i < length(time)) {
      for(j in (i+1):length(time)) {
        distance = (time[i] - time[j])
        covariance = gp_covariance(distance, scale, length, periodic, period)
        cov_matrix[i,j] = covariance
        cov_matrix[j,i] = covariance
      }
    }
  }
  # Draw from the multinormal distribution using the cholesky decomposition of the cov. matrix
  chol_cov = t(chol(cov_matrix));
  raw_profile = chol_cov %*% rnorm(length(time)) + mean_func;
  # Transform to strictly positive values
  if(positive_transform) {
    positive_profile = log1p(exp(raw_profile))
    return(t(positive_profile))
  } else {
    return(t(raw_profile))
  }
}

test_constant_as_in_genexpi <- function(profile) {
  errors <- profile * 0.2;
  errors[errors < 0.5] <- 0.5
  minUpperBound = min(profile + errors)
  maxLowerBound = max(profile - errors)
  return(minUpperBound > maxLowerBound)

}

#Eliminates profiles that are either too flat or too similar (corelation > 0.9) to the regulator profile
generate_useful_random_profile <- function(time, scale, length, original_profile, original_profile_time) {
  numTries = 0;
  repeat { #Rejection sampling to get a non-flat, non-similar profile
    random_profile = generate_random_profile(time, scale, length);
    numTries = numTries + 1
    if(!test_constant_as_in_genexpi(random_profile)) {
      if(!is.null(original_profile)) {
        na_mask = !is.na(original_profile)
        if(numTries > 100) {
          warning("Could not find different profile - scale: ", scale, " length: ", length, " time: ", time);
          break;
        }
        else if(cor(original_profile[na_mask], (random_profile[1,original_profile_time])[na_mask]) < 0.9) {
          break;
        }
      }
      break;
    }
    if(numTries > 200) {
      warning("Could not create non-flat profile. - scale: ", scale, " length: ", length, " time: ", time)
      break;
    }
  }
  return(random_profile)
}


plot_random_profiles <- function(n, time, scale, length, true_time = time, true_profile = NULL, mean_func = 0, periodic = FALSE, period = 1, positive_transform = TRUE) {
  profiles = array(0, c(n, length(time)));
  for(i in 1:n) {
    profiles[i,] = generate_random_profile(time, scale, length, mean_func = mean_func, periodic = periodic, period = period, positive_transform = positive_transform);
  }


  result = ggmatplot(time, t(profiles))
  if(!is.null(true_profile)) {
    result = result + geom_point(data = data.frame(x = true_time, y = true_profile), aes(x=x, y=y))
  }
  return(result)
}


simulate_regulated_spline <- function(num_targets, num_time, measurement_times,
                                             num_df = NULL,
                                             spline_basis = NULL,
                                             measurement_sigma_absolute_prior_sigma = NULL,
                                             measurement_sigma_relative_prior_sigma = NULL,
                                             measurement_sigma_absolute = NULL,
                                             measurement_sigma_relative = NULL,
                                             regulator_profile = NULL,
                                             regulation_signs = NULL,
                                             regulator_measured = TRUE, integrate_ode45 = TRUE) {
  time <- 1:num_time;

  if(is.null(spline_basis)) {
    spline_degree = 3
    spline_basis <- bs(1:num_time, df = num_df, degree = spline_degree)

    num_coeff <- dim(spline_basis)[2]
    if(all(spline_basis[,num_coeff] == 0)) {
      num_coeff <- num_coeff - 1
      spline_basis = spline_basis[,1:num_coeff]
    }
  } else {
    num_coeff <- dim(spline_basis)[2]
  }

  scale <- 5

  #Rejection sampling to get interesting profile
  n_rejections <- 0
  if(is.null(regulator_profile)) {
    repeat {
      coeffs <- rnorm(num_coeff, 0, 1)
      regulator_profile <- array((spline_basis %*% coeffs)  * scale, num_time)
      if(sd(regulator_profile) > 1
         && sd(regulator_profile[1:floor(num_time / 2)]) > 1
         && sd(regulator_profile[ceiling(num_time / 2):num_time]) > 1
         && mean(diff(regulator_profile) > 1e-2) > 0.2 #Avoid monotonicity
         && mean(diff(regulator_profile) < 1e-2) > 0.2
      ) {
        break;
      }
      n_rejections <- n_rejections + 1
    }
    cat(n_rejections, " rejections regulator\n")
  } else {
    coeffs <- array(0, dim(spline_basis)[2])
  }

  min_intercept <-  -min(regulator_profile)
  if(regulator_measured) {
    intercept_prior_sigma <- 1#scale
    intercept <- min_intercept + abs(rnorm(1, 0, intercept_prior_sigma))
  } else {
    intercept <- min_intercept
    intercept_prior_sigma <- 0.01
  }


  if(is.null(measurement_sigma_absolute)) {
    if(!is.null(measurement_sigma_relative)) {
      stop("Both sigma have to be null")
    }
    sigma_given <- FALSE
    measurement_sigma_absolute <- abs(rnorm(1,0, measurement_sigma_absolute_prior_sigma))
    measurement_sigma_relative <- abs(rnorm(1,0, measurement_sigma_relative_prior_sigma))
    measurement_sigma <- measurement_sigma_prior(0, measurement_sigma_absolute_prior_sigma, 0, measurement_sigma_relative_prior_sigma)
  } else {
    sigma_given <- TRUE
    measurement_sigma <- measurement_sigma_given(sigma_absolute = measurement_sigma_absolute, sigma_relative = measurement_sigma_relative)
  }

  expression_true <- array(-1, c(num_targets,num_time))
  expression_observed <- array(-1, c(num_targets,length(measurement_times)))

  asymptotic_normalized_state_prior_sigma <-  0.5
  degradation_prior_mean <- -2
  degradation_prior_sigma <- 1
  mean_regulatory_input_prior_sigma <- 5
  sd_regulatory_input_prior_sigma <- 3
  initial_condition_prior_sigma <- 2

  if(is.null(regulation_signs)) {
    regulation_signs <- rbernoulli(num_targets) * 2 - 1 #Generates randomly -1 or 1
  } else if(length(regulation_signs) == 1) {
    regulation_signs <- array(regulation_signs, num_targets)
  }

  b_centered <- array(-Inf, num_targets)
  b <- array(-Inf, num_targets)
  w <- array(-Inf, num_targets)
  sensitivity <- array(-Inf, num_targets)
  initial_condition <- array(-Inf, num_targets)
  degradation <- array(-Inf, num_targets)
  mean_synthesis <- array(-Inf, num_targets)
  mean_regulatory_input <- array(-Inf, num_targets)
  sd_regulatory_input <- array(-Inf, num_targets)
  asymptotic_normalized_state <- array(-Inf, num_targets)

  n_rejections <- 0
  for(t in 1:num_targets) {
    #Rejection sampling to have interesting profiles
    repeat {
      #The actual parameters
      mean_regulatory_input[t] <- rnorm(1, 0, mean_regulatory_input_prior_sigma)
      sd_regulatory_input[t] <- abs(rnorm(1, 0, sd_regulatory_input_prior_sigma))
      asymptotic_normalized_state[t] <- abs(rnorm(1, 0, asymptotic_normalized_state_prior_sigma))
      degradation[t] <- rlnorm(1, degradation_prior_mean, degradation_prior_sigma)


      #The derived parameters
      w[t] <- regulation_signs[t] * sd_regulatory_input[t] / sd(regulator_profile)
      b_centered[t] <- mean_regulatory_input[t] - w[t] * mean(regulator_profile)
      b[t] <- b_centered[t] + (intercept * w[t])

      synthesis <- 1 / (1 + exp( -regulator_profile * w[t] - b_centered[t]))

      #In the actual model, there is max(expression) instead of scale
      mean_synthesis[t] <- asymptotic_normalized_state[t] * degradation[t] * scale
      sensitivity[t] <- mean_synthesis[t] / mean(synthesis)

      #degradation_over_sensitivity[t] <- rlnorm(1,0,1)
      #degradation[t] <- degradation_over_sensitivity[t] * sensitivity[t];
      #degradation[t] <- abs(rnorm(1, 0, degradation_prior_sigma))

      initial_condition[t] <- abs(rnorm(1, 0,initial_condition_prior_sigma))

      if(integrate_ode45){
        params <- c(degradation = degradation[t], bias = b_centered[t], sensitivity = sensitivity[t], weight = w[t], basal_transcription = 0, protein = approxfun(time, regulator_profile, rule=2));
        expression_true[t,] <-  ode( y = c(x = initial_condition[t]), times = time, func = target_ODE, parms = params, method = "ode45")[,"x"];
      } else {
        expression_true[t,] <- numerical_integration(0, degradation[t], initial_condition[t], sensitivity[t], weight = w[t], bias = b_centered[t], regulator_profile, num_time)
      }

      if(all(expression_true[t,] > 0)
         #&& sd(expression_true[t,]) > 0.5
         #&& any(expression_true[t,] > 2 * measurement_sigma_absolute)
         #&& mean(diff(expression_true[t,]) > 1e-2) > 0.2 #Avoid monotonicity
         #&& mean(diff(expression_true[t,]) < 1e-2) > 0.2
      ) {
        break;
      }
      n_rejections <- n_rejections + 1
    }
    expression_observed[t,] <- rtruncnorm(length(measurement_times),
                                          mean = expression_true[t,measurement_times],
                                          sd =  c(measurement_sigma_absolute) + c(measurement_sigma_relative) * expression_true[t,measurement_times],
                                          a = 0)

    #Hack around the fact that the actual model uses max(expression) to determine parameters (the model is non-generative in this sense)
    asymptotic_normalized_state[t] <- mean_synthesis[t] / (degradation[t] * max(expression_observed[t,]))
  }

  cat(n_rejections, " rejections targets\n")

  regulator_expression_true <- regulator_profile + intercept
  regulator_expression_observed <- rtruncnorm(length(measurement_times),
                                              mean = regulator_expression_true[measurement_times],
                                              sd =  c(measurement_sigma_absolute) + c(measurement_sigma_relative) * regulator_expression_true[measurement_times],
                                              a = 0)

  return(list(
    true = list (
      coeffs = coeffs,
      regulator_profile = regulator_profile,
      regulator_expression = regulator_expression_true,
      w = w,
      b = b,
      initial_condition = initial_condition,
      sensitivity = sensitivity,
      expression = expression_true,
      intercept = intercept,
      degradation = degradation,
      sigma_relative_param = if(sigma_given) { array(measurement_sigma_relative,1) } else { numeric(0) },
      sigma_absolute_param = if(sigma_given) { array(measurement_sigma_absolute,1) } else { numeric(0) }
    ), observed =

      regulated_model_params(
        measurement_times = measurement_times,
        regulator_expression = regulator_expression_observed,
        target_expression = t(expression_observed),
        measurement_sigma = measurement_sigma,
        spline_params = spline_params(spline_basis, scale),
        params_prior = params_prior(
          initial_condition_prior_sigma = initial_condition_prior_sigma,
          asymptotic_normalized_state_prior_sigma = asymptotic_normalized_state_prior_sigma,
          degradation_prior_mean = degradation_prior_mean,
          degradation_prior_sigma = degradation_prior_sigma,
          mean_regulatory_input_prior_sigma = mean_regulatory_input_prior_sigma,
          sd_regulatory_input_prior_sigma = sd_regulatory_input_prior_sigma,
          intercept_prior_sigma = if (regulator_measured) { intercept_prior_sigma } else { NULL }
        ),
        regulation_signs = matrix(regulation_signs, 1, num_targets)
      )
  ))
}


simulate_constant_synthesis <- function(measurement_times, measurement_sigma_absolute, measurement_sigma_relative, integrate_ode45 = TRUE) {



  expression_true <- array(-1, length(measurement_times))
  expression_observed <- array(-1, length(measurement_times))

  synthesis_over_degradation_prior_mean <- 3
  synthesis_over_degradation_prior_sigma <-  3
  degradation_prior_mean <- -3
  degradation_prior_sigma <- 1
  initial_condition_prior_sigma <- 1

  n_rejections <- 0
  #Rejection sampling to have interesting profiles
  repeat {
    synthesis_over_degradation <- rnorm(1,synthesis_over_degradation_prior_mean,synthesis_over_degradation_prior_sigma)
    degradation <- rlnorm(1, degradation_prior_mean, degradation_prior_sigma)
    synthesis = synthesis_over_degradation * degradation

    initial_condition <- abs(rnorm(1, 0,initial_condition_prior_sigma))

    if(integrate_ode45){
      params <- c(degradation = degradation, synthesis = synthesis);
      expression_true <-  ode( y = c(x = initial_condition), times = measurement_times, func = constant_synthesis_ODE, parms = params, method = "ode45")[,"x"];
    } else {
      expression_true <- numerical_integration(synthesis, degradation, initial_condition, 0, 0, 0, 0, max[measurement_times])[measurement_times]
    }

    if(all(expression_true > 0))
    {
      break;
    }
    n_rejections <- n_rejections + 1
  }
  expression_observed <- rtruncnorm(length(measurement_times), mean = expression_true, sd =  measurement_sigma_absolute + measurement_sigma_relative * expression_true, a = 0)

  cat(n_rejections, " rejections targets\n")

  return(list(
    true = list (
      initial_condition = initial_condition,
      synthesis_over_degradation = synthesis_over_degradation,
      degradation = degradation,
      expression = expression_true
    ), observed = list(
      num_measurements = length(measurement_times),
      measurement_times = measurement_times,
      expression = expression_observed,
      measurement_sigma_absolute = measurement_sigma_absolute,
      measurement_sigma_relative = measurement_sigma_relative,
      initial_condition_prior_sigma = initial_condition_prior_sigma,
      synthesis_over_degradation_prior_mean = synthesis_over_degradation_prior_mean,
      synthesis_over_degradation_prior_sigma = synthesis_over_degradation_prior_sigma,
      degradation_prior_mean = degradation_prior_mean,
      degradation_prior_sigma = degradation_prior_sigma
    )
  ))
}


numerical_integration <- function(basal_transcription, degradation, initial_condition, sensitivity, weight, bias, protein, num_time){

  numerical_result = numeric(num_time);
  numerical_result[1] = initial_condition
  integrationTimeStep = 1

  basal_over_degradation = basal_transcription / degradation;

  regulation_input = bias + weight * protein;
  synthesis = integrationTimeStep * ( 1 / (1 + exp(-regulation_input))) #Don't forget to change the line at the end

  residual = -0.5 * synthesis[1];
  degradationPerStep = exp(-degradation * integrationTimeStep)


  for(time in 2:num_time){
    integral = 0;
    integration_end <- time
    for(previousIndex in 1:integration_end)
    {
      if(previousIndex == 1 || previousIndex == integration_end)
      {
        h = 0.5;
      }
      else
      {
        h = 1;
      }
      degradationCoeff = exp(-degradation * (integration_end - previousIndex) * integrationTimeStep)
      integral = integral + h * synthesis[previousIndex] * degradationCoeff;
    }
    numerical_result[time] = basal_over_degradation + (initial_condition - basal_over_degradation) * exp(-degradation * (integration_end)) + sensitivity * integral
  }
  return(numerical_result)
}

