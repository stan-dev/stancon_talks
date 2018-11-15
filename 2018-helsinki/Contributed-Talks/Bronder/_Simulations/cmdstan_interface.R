library(rstan)
library(digest)
library(microbenchmark)

get_cmdstan_param_string <- function(num_samples = 100,
                                      num_warmup  = 100,
                                      seed = 0) {
  s <- paste0(
  " num_samples=", num_samples,
  " num_warmup=", num_warmup,
  " random seed=", seed
  )
  s
}

get_os <- function(){
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  tolower(os)
}

get_samples <- function(stan_file,
                        cmdstan_dir,
                        temp_dir,
                        stan_data,
                        experiment_hash,
                        param_string) {
  model_code <- readLines(stan_file)
  model_name <- paste0("model_", digest(list(model_code, cmdstan_dir)))
  experiment_name <- paste0("model_", experiment_hash) 
  data_file <- paste0(temp_dir, "/", experiment_name, "_data.R")
  samples_file <- paste0(temp_dir, "/", experiment_name, "_samples.csv")
  
  os <- get_os()
  
  if (os == "windows")
    exe_file <- paste0(temp_dir, "/", model_name, ".exe")
  else
    exe_file <- paste0(temp_dir, "/", model_name)
    
  model_file <- paste0(temp_dir, "/", model_name, ".stan")

  if (!file.exists(model_file)) writeLines(model_code, con = model_file)

  # setup dataset ----
  nms <- c()
  for (i in 1:length(stan_data)) {
    nm <- names(stan_data)[i]
    nms <- c(nms, nm)
    eval(parse(text = paste0(nm, " <- stan_data[[i]]")))
  }
  rstan::stan_rdump(nms, file = data_file)

  # compile Stan model ----
  if (!file.exists(exe_file)) system(paste0("make --directory=", cmdstan_dir, " ", exe_file))

  # sample from model ----
  x <- microbenchmark(
    system(paste0(
      exe_file,
      " method=sample",
      "", param_string,
      " data file=", data_file,
      " output file=", samples_file
    )),
    times = 1
  )

  list(
    samples = read.csv(samples_file, comment.char = "#"),
    sampling_time = x$time / 1000000000
  )
}
