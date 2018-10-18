# Fit model

fitOneStepAhead <- function(toggle_dynamic, 
                            sm_string = "models/epl_model.stan",
                            model_short_name = "default_short", 
                            model_long_name  = "default_long",
                            stanfit_include_pars = NA,
                            fitpars) {
  # unpack fit parameters
  nsamples <- fitpars$nsamples 
  nshifts_start <- fitpars$nshifts_start
  chains <- fitpars$chains 
  warmup <- fitpars$warmup
  nshifts_start <- fitpars$nshifts_start
  nshifts_end <- fitpars$nshifts_end
  init_r_val <- fitpars$init_r_val
  start_round <- fitpars$start_round
  end_round <- fitpars$end_round
  prev_perf_season <- fitpars$prev_perf_season
  target_folder <- fitpars$target_folder
  
  fullrun <- toggle_dynamic
  
  if(fullrun) {
    # compile stan model
    sm <- stan_model(sm_string)
    
    for(shift in 0:(nshifts_end - nshifts_start)){
      model_data <- Create_model_data_for_TS2(NL_ALL[round >= start_round & 
                                                       round < (end_round + nshifts_start + shift)], # in-sample
                                              NL_ALL[season == prev_perf_season], # prev_perf
                                              NL_ALL[round == (end_round + nshifts_start + shift)]) # prediction
      
      fit <- sampling(sm, chains = chains, 
                      init_r = init_r_val, 
                      warmup = warmup,
                      iter = nsamples,
                      save_warmup = FALSE,
                      include = T, # save ONLY specified pars in stanfit
                      pars = stanfit_include_pars,
                      control = list(max_treedepth = 15,
                                     adapt_delta = 0.99),
                      data = model_data)
      
      saveRDS(fit, paste(target_folder, model_short_name, "_", nshifts_start + shift, ".rds", sep=""))
    }
  } 
  model_res <- ReadFitsandCalculateRPS(n_shifts_start = nshifts_start,
                                       n_shifts_end = nshifts_end, 
                                       model_name = model_long_name,
                                       model_short = model_short_name,
                                       start_round = start_round,
                                       end_round = end_round,
                                       prev_perf_season = prev_perf_season,
                                       fitdir = target_folder
                       )
  model_res
}

