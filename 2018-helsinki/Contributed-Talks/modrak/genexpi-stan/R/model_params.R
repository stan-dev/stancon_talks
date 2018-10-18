
#IDEA: function print_assumptions that says what assumptions does the model make

measurement_sigma_given <- function(sigma_absolute, sigma_relative) {
  list(
    sigma_given = 1,
    sigma_absolute_data = array(sigma_absolute,1),
    sigma_relative_data = array(sigma_relative,1),
    sigma_relative_prior_mean = numeric(0),
    sigma_absolute_prior_mean = numeric(0),
    sigma_relative_prior_sigma = numeric(0),
    sigma_absolute_prior_sigma = numeric(0)
  )
}

measurement_sigma_prior <- function(sigma_absolute_prior_mean,
                                    sigma_absolute_prior_sigma,
                                    sigma_relative_prior_mean,
                                    sigma_relative_prior_sigma) {
  list(
    sigma_given = 0,
    sigma_absolute_data = numeric(0),
    sigma_relative_data = numeric(0),
    sigma_relative_prior_mean = array(sigma_relative_prior_mean,1),
    sigma_absolute_prior_mean = array(sigma_absolute_prior_mean,1),
    sigma_relative_prior_sigma = array(sigma_relative_prior_sigma,1),
    sigma_absolute_prior_sigma = array(sigma_absolute_prior_sigma,1)
  )
}

params_prior <- function(
  initial_condition_prior_sigma,
  asymptotic_normalized_state_prior_sigma,
  mean_regulatory_input_prior_sigma,
  sd_regulatory_input_prior_sigma,
  degradation_prior_mean,
  degradation_prior_sigma,
  intercept_prior_sigma = NULL
) {
  list(
    initial_condition_prior_sigma = initial_condition_prior_sigma,
    asymptotic_normalized_state_prior_sigma = asymptotic_normalized_state_prior_sigma,
    mean_regulatory_input_prior_sigma = mean_regulatory_input_prior_sigma,
    sd_regulatory_input_prior_sigma = sd_regulatory_input_prior_sigma,
    degradation_prior_mean = degradation_prior_mean,
    degradation_prior_sigma = degradation_prior_sigma,

    intercept_prior_sigma = if (is.null(intercept_prior_sigma)) { numeric(0) } else {array(intercept_prior_sigma, 1)}
  )
}

coeffs_prior_default <- function(num_spline_basis) {
  list(
    coeffs_prior_given = 0,
    coeffs_prior_mean = matrix(0,0,num_spline_basis),
    coeffs_prior_cov = array(0, c(0, num_spline_basis, num_spline_basis))
  )

}

coeffs_prior_given <- function(coeffs_prior_mean, coeffs_prior_cov) {
  list(
    coeffs_prior_given = 1,
    coeffs_prior_mean = coeffs_prior_mean,
    coeffs_prior_cov = coeffs_prior_cov
  )
}

coeffs_prior_from_fit <- function(fit, covariance_scale = 0.5) {
  samples_coeffs <- rstan::extract(fit,"coeffs")$coeffs
  num_regulators <- dim(samples_coeffs)[3]
  num_spline_basis <- dim(samples_coeffs)[2]

  coeffs_prior_mean <- array(-1, c(num_regulators,num_spline_basis))
  coeffs_prior_cov <- array(-1, c(num_regulators,num_spline_basis, num_spline_basis))

  for(r in 1:num_regulators) {
    coeffs_prior_cov[r,,] <- covariance_scale * cov(samples_coeffs[,,r])
    coeffs_prior_mean[r,] <- colMeans(samples_coeffs[,,r])
  }

  coeffs_prior_given(
    coeffs_prior_mean = coeffs_prior_mean,
    coeffs_prior_cov = coeffs_prior_cov
  )

}

spline_params <- function(
  spline_basis,
  scale
) {

  num_spline_basis <- dim(spline_basis)[2]

  list(
    num_spline_basis = num_spline_basis,
    spline_basis = spline_basis,
    scale = scale
  )
}

force_to_expression_matrix <- function(x) {
  if(is.matrix(x) || length(dim(x) == 2)) {
    x
  } else {
    as.matrix(x, 1, length(x))
  }
}

regulated_model_params <- function(
  measurement_times,
  measurement_sigma,
  spline_params,
  params_prior,
  regulation_signs = NULL,
  target_expression = NULL,
  regulator_expression = NULL,
  coeffs_prior = coeffs_prior_default(spline_params$num_spline_basis)
) {
  #TODO Checks:
  #expression has to be matrices (two dims)
  #measurement_times are integers
  #regulation signs is a matrix
  #coeffs prior, measurement_sigma, params_prior are of correct type

  if(is.null(regulator_expression) && is.null(target_expression)) {
    stop("Some expression data has to be given")
  }

  if(!is.null(regulator_expression)) {
    regulator_expression <- force_to_expression_matrix(regulator_expression)
  }

  if(!is.null(target_expression)) {
    target_expression <- force_to_expression_matrix(target_expression)
  }

  if(is.null(regulation_signs)) {
    if(!is.null(target_expression)) {
      stop("When regulation_signs are not given, target_expression must not be set")
    }
    num_targets <- 0
    num_regulators <- ncol(regulator_expression)
    regulation_signs <- matrix(0, num_regulators, 0)
  } else {
    num_targets <- ncol(regulation_signs)
    num_regulators <- nrow(regulation_signs)
  }

  if(is.null(regulator_expression) && length(params_prior$intercept_prior_sigma) > 0) {
    params_prior$intercept_prior_sigma <- numeric(0)
  }

  if(!is.null(regulator_expression) && length(params_prior$intercept_prior_sigma) == 0) {
    stop("When regulator_expression is given, the params_prior has to contain intercept_prior_sigma")
  }


  base_params <- list(
    num_time = max(measurement_times),
    num_measurements = length(measurement_times),
    num_regulators = num_regulators,
    num_targets = num_targets,
    regulators_measured = if (is.null(regulator_expression)) { 0 } else { 1 },


    measurement_times = measurement_times,

    regulator_expression = if (is.null(regulator_expression)) { matrix(0, 0, num_regulators) } else { regulator_expression },
    expression = if (is.null(target_expression)) { matrix(0, length(measurement_times), 0)} else { target_expression },
    regulation_signs = regulation_signs
  )

  unlist(list(base_params, spline_params, params_prior, measurement_sigma, coeffs_prior), recursive = FALSE)
}

csynth_model_params <- function(measurement_times, measurement_sigma, params_prior, target_expression){
  if(measurement_sigma$sigma_given != 1) {
    stop("Constant synthesis model requires sigma to be given")
  }
  base_params <- list(
    num_measurements = length(measurement_times),
    measurement_times = measurement_times,
    expression = target_expression,
    measurement_sigma_absolute = measurement_sigma$sigma_absolute_data[1],
    measurement_sigma_relative = measurement_sigma$sigma_relative_data[1]
  )

  unlist(list(base_params, params_prior), recursive = FALSE)
}
