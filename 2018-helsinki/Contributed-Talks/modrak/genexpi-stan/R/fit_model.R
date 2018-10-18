regulated_init_fun <- function(data){
  function(chain_id = 0) {
    initial_condition <- data$expression[1,]
    initial_condition[initial_condition < 1e-12] <- 1e-12
    list(
      initial_condition = array(initial_condition, data$num_targets)
    )
  }
}


regulated_control <- list(adapt_delta = 0.95)

fit_regulated <- function(data, model = stan_model(file = here::here('stan', 'regulated.stan')), ...) {
  rstan::sampling(model, data = data,
                  init = regulated_init_fun(data), control = regulated_control,
                  ...)
}

fit_regulated_multi <- function(data, output.dir,
                                       model = stan_model(file = here::here('stan', 'regulated.stan')),
                                       chains = 4, ids_to_compute = 1:length(data),
                                       ...) {
  #init_per_item <- lapply(data, function(x) { lapply(1:chains, regulated_init_fun(x)) })
  init_per_item <- lapply(data, regulated_init_fun)
  sampling_multi(model = model, data = data, output.dir = output.dir,
                 init_per_item = init_per_item, control = regulated_control,
                 chains = chains, ids_to_compute = ids_to_compute,
                 ...)
}


fit_csynth <- function(data, model = stan_model(file = here::here('stan', 'constant_synthesis.stan')), ...) {
  rstan::sampling(model, data = data,
                  ...)
}

fit_csynth_multi <- function(data, output.dir,
                                model = stan_model(file = here::here('stan', 'constant_synthesis.stan')),
                                chains = 4, ids_to_compute = 1:length(data),
                                ...) {
  sampling_multi(model = model, data = data, output.dir = output.dir,
                 chains = chains, ids_to_compute = ids_to_compute,
                 ...)
}


get_log_lik_genexpi <- function(fit, target = 1) {
  samples_log_lik <- extract_log_lik(fit, parameter_name = "log_likelihood", merge_chains = FALSE)
  if(length(dim(samples_log_lik)) == 3 && target == 1) {
    samples_log_lik
  } else {
    stop("Not implemented yet - how to extract log_lik for multiple targets")
    samples_log_lik[,,target]
  }
}

get_log_lik_csynth <- function(fit) {
  extract_log_lik(fit, parameter_name = "log_likelihood", merge_chains = FALSE)
}

get_loo_genexpi <- function(fit, target = 1) {
  log_lik <- get_log_lik_genexpi(fit, target)
  loo(log_lik, r_eff = relative_eff(exp(log_lik)))
}

get_loo_csynth <- function(fit) {
  log_lik <- get_log_lik_csynth(fit)
  loo(log_lik, r_eff = relative_eff(exp(log_lik)))
}

process_fit_multi <- function(results, fun, cores = parallel::detectCores()) {
  cl <- parallel::makeCluster(cores, useXDR = FALSE)
  on.exit(parallel::stopCluster(cl))


  process <- function(id) {
    fit <- sampling_multi_read_fit(results, id)
    fun(fit)
  }

  dependencies <- c("rstan","Rcpp","genexpiStan","loo")
  .paths <- unique(c(.libPaths(), sapply(dependencies, FUN = function(d) {
    dirname(system.file(package = d))
  })))
  .paths <- .paths[.paths != ""]
  parallel::clusterExport(cl, varlist = ".paths", envir = environment())
  parallel::clusterEvalQ(cl, expr = .libPaths(.paths))
  parallel::clusterEvalQ(cl, expr =
                           suppressPackageStartupMessages(require(rstan, quietly = TRUE)))
  parallel::clusterEvalQ(cl, expr =
                           suppressPackageStartupMessages(require(loo, quietly = TRUE)))
  parallel::clusterEvalQ(cl, expr =
                           suppressPackageStartupMessages(require(genexpiStan, quietly = TRUE)))

  parallel::clusterExport(cl, varlist = c("results"), envir = environment())

  parallel::parLapplyLB(cl, X = 1:length(results), fun = process)
}

waic_estimate <- function(log_lik) {
  waic(log_lik)$estimates["waic","Estimate"]
}

process_log_lik_regulated_multi <- function(results, target = 1, fun, cores = parallel::detectCores()) {

  process_fun <- function(fit) {
    fun(get_log_lik_regulated(fit, target))
  }

  process_fit_multi(results, process_fun, cores)
}

process_log_lik_csynth_multi <- function(results, fun = waic, cores = parallel::detectCores()) {

  process_fun <- function(fit) {
    fun(get_log_lik_csynth(fit))
  }

  process_fit_multi(results, process_fun, cores)
}
