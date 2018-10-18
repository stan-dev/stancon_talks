FitStanModel <- function(stan_model_file, fit_pars, model_data, output_file, fullrun = 0, n_samples = 100){
  if(fullrun) {
    sm <- stan_model(stan_model_file)
    stanfit <-  sampling(sm, 
                       chains = fit_pars$chains, 
                       init_r = fit_pars$init_r_val, 
                       warmup = n_samples/2,
                       iter = n_samples,
                       data = model_data,
                       save_warmup = FALSE,
                       control = list(max_treedepth = 15,
                                      adapt_delta = 0.99)
    )
    saveRDS(stanfit, output_file)
    return(stanfit)
  } else { stanfit <- readRDS(output_file)
            return(stanfit)
          }
}