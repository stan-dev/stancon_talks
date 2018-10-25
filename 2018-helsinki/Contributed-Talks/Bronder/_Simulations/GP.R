# directory information ----
setwd("D:/Projects/stancon2018")

cmdstan_dir <- "D:/Projects/cmdstan"

output_dir <- "./_Simulations/_Output"

# cmdstan interface ----
source("./_Simulations/cmdstan_interface.R")

# data generation script ----
source("./_Simulations/GP_data_generator.R")

# globals ----
notes_cpu <- "Intel i5-6600K"
notes_gpu <- "NVIDIA GTX 1080 Ti"

# flags for disabling/enabling simulations
run_cpu <- TRUE
run_gpu <- TRUE

# dataset sizes ----
N_range <- c(64, 128, 256, 512, 1024, 2048, 5096)

# warmup ----
N_w <- 300

# samples ----
N_s <- 200

# seed generation ----
set.seed(0)
seeds <- sample.int(.Machine$integer.max, length(N_range) * 2)
seed_index <- 1

for (N in N_range) {
  # seeds ----
  data_seed <- seeds[seed_index]
  stan_seed <- seeds[seed_index + 1]
  seed_index <- seed_index + 2
  
  # generate data ----
  set.seed(data_seed)
  stan_data <- generate_GP_data(N)
  
  # simulate cpu ----
  # param string ----
  param_string <- get_cmdstan_param_string(
    num_samples = N_s,
    num_warmup = N_w,
    seed = stan_seed
  )
  
  # Stan hash experiment settings ----
  experiment_hash_cpu <- digest(list(stan_data, param_string, notes_cpu))
  experiment_hash_gpu <- digest(list(stan_data, param_string, notes_gpu))
  
  # data files ----
  fn_stan_cpu <- paste0(output_dir, "/GP/GP_", experiment_hash_cpu, ".rds")
  fn_stan_gpu <- paste0(output_dir, "/GP/GP_", experiment_hash_gpu, ".rds")
  
  # run sampling cpu ----
  # don't run if n is large
  if (run_cpu)
  {
    print(paste0("---- Simulating on CPU, N = ", N, " ----"))
    
    res_stan <- get_samples(
      stan_file = normalizePath("./_Simulations/_Models/gp.stan", winslash = "/"),
      cmdstan_dir,
      normalizePath(paste0(output_dir, "/GP/_Temp"), winslash = "/"),
      stan_data,
      experiment_hash_cpu,
      param_string
    )
    
    # get parameters
    stan_fit <- data.frame(rho = mean(res_stan$samples$rho),
                           alpha = mean(res_stan$samples$alpha),
                           sigma = mean(res_stan$samples$sigma))
    
    saveRDS(
      list(
        time = res_stan$sampling_time,
        N = N,
        num_warmup = N_w,
        num_samples = N_s,
        GPU = FALSE,
        parameters_fit = stan_fit,
        notes = notes_cpu,
        y = stan_data$y,
        x = stan_data$x,
        x_predict = stan_data$x_predict,
        y_predict = colMeans(res_stan$samples[, -c(1:(10 + N))]),
        date = Sys.time()
      ),
      file = fn_stan_cpu
    )
  }
  
  # run sampling gpu ----
  if (run_gpu)
  {
    print(paste0("---- Simulating on GPU, N = ", N, " ----"))
    
    res_stan <- get_samples(
      stan_file = normalizePath("./_Simulations/_Models/gp_gpu.stan", winslash = "/"),
      cmdstan_dir,
      normalizePath(paste0(output_dir, "/GP/_Temp"), winslash = "/"),
      stan_data,
      experiment_hash_gpu,
      param_string
    )
    
    # get parameters
    stan_fit <- data.frame(rho = mean(res_stan$samples$rho),
                           alpha = mean(res_stan$samples$alpha),
                           sigma = mean(res_stan$samples$sigma))
    
    saveRDS(
      list(
        time = res_stan$sampling_time,
        N = N,
        num_warmup = N_w,
        num_samples = N_s,
        GPU = TRUE,
        parameters_fit = stan_fit,
        notes = notes_gpu,
        y = stan_data$y,
        x = stan_data$x,
        x_predict = stan_data$x_predict,
        y_predict = colMeans(res_stan$samples[, -c(1:(10 + N))]),
        date = Sys.time()
      ),
      file = fn_stan_gpu
    )
  }
}