sampling_multi_filename <- function(output.dir, id, chain) {
  file.path(output.dir, paste0("stan_out_",id,"_c",chain,".csv"))
}

sampling_multi_get_results <- function(count, output.dir, chains = 4) {
  results <- list()
  for(id in 1:count) {
    files <- character(chains)
    for(chain in 1:chains) {
      files[chain] <- sampling_multi_filename(output.dir, id, chain)
    }
    results[[id]] <- files
  }
  results
}

sampling_multi <- function(model, data, output.dir, chains = 4, cores = parallel::detectCores(), init_per_item = NULL, ids_to_compute = 1:length(data),  ...) {
  cl <- parallel::makeCluster(cores, useXDR = FALSE)
  on.exit(parallel::stopCluster(cl))

  if(!dir.exists(output.dir)) {
    warning(paste0("The directory ", output.dir, " does not exist, creating."))
    dir.create(output.dir)
    if(!dir.exists(output.dir)) {
      stop("Could not create output directory")
    }
  }

#  objects <- ls()
  fit_fun <- function(i) {
    id <- floor( (i - 1) / chains) + 1
    chain_id <- ((i - 1) %% chains) + 1
    .dotlist$chain_id <- chain_id

    out_file <- sampling_multi_filename(output.dir, id, chain_id)
    if(is.list(init_per_item)) {
      .dotlist$init <- init_per_item[[id]]
    }

    if(is.list(.dotlist$init)) .dotlist$init <- .dotlist$init[chain_id]

    .dotlist$sample_file <- out_file
    .dotlist$data <- data[[id]]
    out <- do.call(rstan::sampling, args = .dotlist)
    if (out@mode == 1L || out@mode == 2L) {
      #Sampling not conducted - test_grad or error
      return(NULL)
    } else {
      return(out_file)
    }
  }

#  .dotlist <- c(sapply(objects, simplify = FALSE, FUN = get,
#                       envir = environment()), list(...))
  .dotlist <- c(list(model, chains = 1L, cores = 0L), list(...))


  dependencies <- c("rstan", "Rcpp", "genexpiStan")
  .paths <- unique(c(.libPaths(), sapply(dependencies, FUN = function(d) {
    dirname(system.file(package = d))
  })))
  .paths <- .paths[.paths != ""]
  parallel::clusterExport(cl, varlist = ".paths", envir = environment())
  parallel::clusterEvalQ(cl, expr = .libPaths(.paths))
  parallel::clusterEvalQ(cl, expr =
                           suppressPackageStartupMessages(require(rstan, quietly = TRUE)))

  parallel::clusterExport(cl, varlist = ".dotlist", envir = environment())

  parallel::clusterExport(cl, varlist = c("data", "output.dir","init_per_item"), envir = environment())

  ids <- rep(ids_to_compute, each = chains)
  items <- ((rep(ids_to_compute, each = chains) - 1) * chains  + rep(1:chains, times = length(ids_to_compute)))

  results_flat <-  parallel::parSapplyLB(cl, X = items, FUN = fit_fun)

  results <- list()
  for(i in 1:length(data)) {
    results[[i]] <- results_flat[ids == i]
  }
  results
}


sampling_multi_read_fit <- function(results, index) {
  read_stan_csv(results[[index]])
}
