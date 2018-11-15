CombineDataWithPredictions <- function(NL_ALL = NL_ALL, 
                                       fit_pars = fit_pars,
                                       osa_res = osa_res,
                                       model_nr, osa_i){

  start_round <- fit_pars$start_round
  end_round <- fit_pars$end_round
  nshifts_start <- fit_pars$nshifts_start
  nshifts_end <- fit_pars$nshifts_end
  prev_perf_season <- fit_pars$prev_perf_season
  rps_index <- 1
  for(shift in 0:(nshifts_end - nshifts_start)){
    model_data <- Create_model_data_for_TS2(NL_ALL[round >= start_round & 
                                                     round < (end_round + nshifts_start + shift)], # in-sample
                                            NL_ALL[season == prev_perf_season], # prev_perf
                                            NL_ALL[round == (end_round + nshifts_start + shift)]) 
    data <- NL_ALL[round == (end_round + nshifts_start + shift)]
    # prediction 
    predictions <- merge(osa_res[[osa_i]][[6]][[shift + 1]], 
                         osa_res[[osa_i]][[7]][[shift + 1]],
                         by = "game_id")
    predictions <- data.table(predictions, 
                              home_team_id = model_data$home_team_pred,
                              away_team_id = model_data$away_team_pred,
                              goal_difference = model_data$pred_goal_difference,
                              rps_vec = osa_res[[osa_i]][[1]][rps_index:(rps_index + 
                                                                              length(model_data$pred_goal_difference)-1)])
    rps_index <- rps_index + length(model_data$pred_goal_difference)
    
    setkey(model_data$id_lut, home_team_id)
    setkey(predictions, home_team_id)
    predictions <- model_data$id_lut[predictions]
    setnames(model_data$id_lut, "home_team_id", "away_team_id")
    setnames(model_data$id_lut, "HomeTeam", "AwayTeam")
    setkey(model_data$id_lut, away_team_id)
    setkey(predictions, away_team_id)
    predictions <- model_data$id_lut[predictions]
    # unieke key wordt HomeTeam_AwayTeam, pred_goal_difference 
    predictions <- predictions[, matchKey := paste(HomeTeam, AwayTeam, goal_difference, sep = '')]
    data <- data[, matchKey := paste(HomeTeam, AwayTeam, goal_difference, sep = '')]
    setkey(predictions, matchKey)
    setkey(data, matchKey)
    data <- predictions[,.(rps_vec, p_win,    p_draw,    p_loss, act_win, act_draw, act_loss, matchKey)][data]
    if(shift == 0){ data_w_predictions = data} else {data_w_predictions = rbind(data_w_predictions, data)}
  }
  data_w_predictions <- data_w_predictions[, model_nr := model_nr]
  return(data_w_predictions)
}