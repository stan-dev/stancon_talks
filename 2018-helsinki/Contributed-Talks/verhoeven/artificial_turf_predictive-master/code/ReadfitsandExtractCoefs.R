ReadFitsWrapper <- function(toggle_dynamic, 
                            sm_string = "models/epl_model.stan",
                            model_short_name = "default_short", 
                            model_long_name  = "default_long",
                            stanfit_extract_pars = NA,
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
  
  model_res <- ReadFitsandExtractCoefs(n_shifts_start = nshifts_start,
                                       n_shifts_end = nshifts_end, 
                                       model_name = model_long_name,
                                       model_short = model_short_name,
                                       start_round = start_round,
                                       end_round = end_round,
                                       prev_perf_season = prev_perf_season,
                                       fitdir = target_folder,
                                       stanfit_extract_pars = stanfit_extract_pars
  )
  model_res
}

ReadFitsandExtractCoefs <- function(n_shifts_start = NA,
                                    n_shifts_end = NA, 
                                    model_name, 
                                    model_short, 
                                    start_round = NA,
                                    end_round = NA,
                                    prev_perf_season = NA,
                                    fitdir = "c:/testversleutel/FITS/",
                                    stanfit_extract_pars = NA){
  
  model_fit <- c()
 
  for(shift in 0:(n_shifts_end - n_shifts_start)){
    #print(".")
    model_data <- Create_model_data_for_TS2(NL_ALL[round >= start_round & 
                                                     round < (end_round+n_shifts_start+shift)], 
                                            NL_ALL[season == prev_perf_season], 
                                            NL_ALL[round == (end_round+n_shifts_start+shift)]) 
    fit <- readRDS(paste(fitdir, model_short, "_", n_shifts_start + shift, ".rds", sep=""))
    sims <- extract(fit, pars = c(stanfit_extract_pars))
    sims <- data.table(melt(sims) )
    sims <- sims[, .(Q5 = quantile(value, 0.05),
                              Q50 = quantile(value, 0.5),
                              Q95 = quantile(value, 0.95)), .(L1)]
    # vector equal to n_samples
    sims <- sims[, shift_nr := shift]
    sims <- sims[, model_name := model_short]
    if(shift == 0){ model_fit <- sims} else { model_fit <- rbind(model_fit ,sims)}
    # extract posterior probs for parameters of interest contained in stanfit_extract_pars
    # add parameters of interest to model_list
  }
  model_fit
}

tidyCoefs <- function(coefs, my_model_list){
  tidy_res <- c()
  for (i in 1:nrow(my_model_list)){
    if(i == 1){tidy_res <- data.frame(model_nr = my_model_list[i,]$model_nr,
                                      coefs[[i]])
    } else {tidy_res <- rbind(tidy_res, data.frame(model_nr = my_model_list[i,]$model_nr,
                                                   coefs[[i]]))}
  }
  tidy_res <- data.table(tidy_res)
  tidy_res <- tidy_res[, model_name := NULL]
  setkey(tidy_res, model_nr)
  setkey(my_model_list, model_nr)
  tidy_res <- my_model_list[tidy_res]
  return(tidy_res)
}
  
ReadRhatNeffWrapper <- function(toggle_dynamic, 
                            sm_string = "models/epl_model.stan",
                            model_short_name = "default_short", 
                            model_long_name  = "default_long",
                            stanfit_extract_pars = NA,
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
  
  model_res <- ReadFitsandExtractRhatNeff(n_shifts_start = nshifts_start,
                                       n_shifts_end = nshifts_end, 
                                       model_name = model_long_name,
                                       model_short = model_short_name,
                                       start_round = start_round,
                                       end_round = end_round,
                                       prev_perf_season = prev_perf_season,
                                       fitdir = target_folder,
                                       stanfit_extract_pars = stanfit_extract_pars
  )
  model_res
}  

ReadFitsandExtractRhatNeff <- function(n_shifts_start = NA,
                                    n_shifts_end = NA, 
                                    model_name, 
                                    model_short, 
                                    start_round = NA,
                                    end_round = NA,
                                    prev_perf_season = NA,
                                    fitdir = "c:/testversleutel/FITS/",
                                    stanfit_extract_pars = NA){
  
  model_fit <- c()
  
  for(shift in 0:(n_shifts_end - n_shifts_start)){
    #print(".")
    model_data <- Create_model_data_for_TS2(NL_ALL[round >= start_round & 
                                                     round < (end_round+n_shifts_start+shift)], 
                                            NL_ALL[season == prev_perf_season], 
                                            NL_ALL[round == (end_round+n_shifts_start+shift)]) 
    fit <- readRDS(paste(fitdir, model_short, "_", n_shifts_start + shift, ".rds", sep=""))
    s <- summary(fit)
    Rhat <- s$summary[,'Rhat']
    n_eff <- s$summary[,'n_eff']
    varnames <- names(n_eff)
    sims <- data.table(varnames, n_eff, Rhat)
    sims <- sims[, shift_nr := shift]
    sims <- sims[, model_name := model_short]
    if(shift == 0){ model_fit <- sims} else { model_fit <- rbind(model_fit ,sims)}
    # extract posterior probs for parameters of interest contained in stanfit_extract_pars
    # add parameters of interest to model_list
  }
  model_fit
}

# tidyRhatNeff <- function(coefs, my_model_list){
#   tidy_res <- c()
#   for (i in 1:nrow(my_model_list)){
#     if(i == 1){tidy_res <- data.frame(model_nr = my_model_list[i,]$model_nr,
#                                       coefs[[i]])
#     } else {tidy_res <- rbind(tidy_res, data.frame(model_nr = my_model_list[i,]$model_nr,
#                                                    coefs[[i]]))}
#   }
#   tidy_res <- data.table(tidy_res)
#   tidy_res <- tidy_res[, model_name := NULL]
#   setkey(tidy_res, model_nr)
#   setkey(my_model_list, model_nr)
#   tidy_res <- my_model_list[tidy_res]
#   return(tidy_res)
# }