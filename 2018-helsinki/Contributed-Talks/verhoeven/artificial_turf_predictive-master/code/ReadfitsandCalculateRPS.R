# need RPS per match as well as ARPS for each model
# + Calculate bootstrapped intervals on the RPS
# store also the predictions and actual outcome of the predicted matches
ReadFitsandCalculateRPS <- function(n_shifts_start = NA,
                                    n_shifts_end = NA, 
                                    model_name, 
                                    model_short, 
                                    start_round = NA,
                                    end_round = NA,
                                    prev_perf_season = NA,
                                    fitdir = "c:/testversleutel/FITS/"){
  aprs_vec <- c()
  rps_vec <- c()
  model_fit <- c()
  pred_probz <- list()
  actual_scorez <- list()
  
  for(shift in 0:(n_shifts_end - n_shifts_start)){
    #print(".")
    model_data <- Create_model_data_for_TS2(NL_ALL[round >= start_round & 
                                                     round < (end_round+n_shifts_start+shift)], 
                                          NL_ALL[season == prev_perf_season], 
                                          NL_ALL[round == (end_round+n_shifts_start+shift)]) 
    fit <- readRDS(paste(fitdir, model_short, "_", n_shifts_start + shift, ".rds", sep=""))
    sims <- extract(fit)
    # vector equal to n_samples
    pred_score_diff_rep <- sims$goal_difference_pred_rep
    
    # not used anymore, use rps_vec
    #aprs_vec[shift+1] <- calculate_arps_wrap(pred_score_diff_rep, model_data$pred_goal_difference)
    
    #this calculates average probabilities for each match over all parameters (samples)
    pred_probz[shift+1] <- list(extract_pred_probs(pred_score_diff_rep))
    actual_scorez[shift+1] <- list(Convert_actual_to_win_draw_loss_vector(model_data$pred_goal_difference))
    # calculate rps per match
    if(shift == 0) {
      rps_vec <- calculate_rps_wrap(pred_score_diff_rep, model_data$pred_goal_difference)} else {
        # shift not equal to zero, add to previous
      rps_vec <- c(rps_vec, calculate_rps_wrap(pred_score_diff_rep, model_data$pred_goal_difference))
    }
  }
  print("*")
  # Bootstrap mean RPS
  require(boot)
  # function to obtain our statistic (the mean RPS aka ARPS) from the data
  arps <- function(data, indices) {
    d <- data[indices] # allows boot to select sample
    arps <- mean(d)
    return(arps)
  }
  # bootstrapping with 1000 replications
  results <- boot(data=rps_vec, statistic=arps,
                  R=1000)

  model_fit <- list(rps_vec = rps_vec,
                    aprs = mean(rps_vec, na.rm = T),
                    aprs_ci = boot.ci(results, type="bca"), # get 95% confidence interval
                    name = model_name,
                    short_name = model_short,
                    pred_probz = pred_probz,
                    actual_scorez = actual_scorez
  )
  model_fit
}

