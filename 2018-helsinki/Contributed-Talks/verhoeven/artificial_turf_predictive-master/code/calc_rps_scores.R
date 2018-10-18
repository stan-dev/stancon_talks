
#sims <- extract(stanfit_epl_dpl)
#score_diff_rep <- sims$score_diff_rep

extract_pred_probs <- function(samples){
  # calculate probabilities from samples
  pred_probs <- data.table(melt(samples))
  
  setnames(pred_probs, "Var2", "game_id")
  
  # convert the predicted prob this to p_win, p_draw, p_loss
  # first on a per sample per match: convert prob to outcome
  pred_probs[, win := 0]
  pred_probs[value >= 0.5, win := 1]
  pred_probs[, loss := 0]
  pred_probs[value < -0.5, loss := 1]
  pred_probs[, draw := 0]
  pred_probs[value >= -0.5 & value < 0.5, draw := 1]
  
  # then take mean over predicted outcomes for all MCMC samples per match
  pred_probs <- pred_probs[, .(p_win = mean(win),
                                    p_draw = mean(draw),
                                    p_loss = mean(loss)), .(game_id)]
  
  return(pred_probs)
}

Convert_actual_to_win_draw_loss_vector <- function(actual) {
  actual_scores <- data.table(game_id = 1:length(actual), 
                            true_val = actual)

  # prep for loss function, actual result as (1,0,0), (0,1,0) or (0,0,1)
  actual_scores[, act_win := 0]
  actual_scores[true_val > 0, act_win := 1]
  actual_scores[, act_draw := 0]
  actual_scores[true_val == 0, act_draw := 1]
  actual_scores[, act_loss := 0]
  actual_scores[true_val < 0, act_loss := 1]
  
  return(actual_scores)
}

# quadratic loss for three categories, default scaling is `1` (i.e. Lit multiplies by 0.5, Hattum does not)
calculate_average_brier <- function(pred_probs, actual_scores, scale = 1){
  setkey(pred_probs, game_id)
  setkey(actual_scores, game_id)
  pred_probs <- pred_probs[actual_scores]
  
  rps <- pred_probs[, .(arps = mean((act_win - p_win)^2 +
                        (act_draw - p_draw)^2 +
                        (act_loss - p_loss)^2))]
  return(round(rps * scale, 4))
}  

calculate_brier_wrap <- function(samples, actual)  {
  pred_probz <- extract_pred_probs(samples)
  actual_scorez <- Convert_actual_to_win_draw_loss_vector(actual)
  return(calculate_average_brier(pred_probz, actual_scorez))
}

calculate_arps <- function(pred_probs, actual_scores){
  observed <- c()
  for(i in 1:nrow(actual_scores)){
    observed[i] <- which.max(actual_scores[i,])
  }
  rps <- rankProbScore(pred_probs, observed)
  return(round(mean(rps), 4))
}  

# no averaging, use to test function on examples in Paper
calculate_rps <- function(pred_probs, actual_scores){
  observed <- c()
  for(i in 1:nrow(actual_scores)){
    observed[i] <- which.max(actual_scores[i,])
  }
  rps <- rankProbScore(pred_probs, observed)
  return(round(rps, 6))
}  

calculate_arps_wrap <- function(samples, actual)  {
  pred_probz <- extract_pred_probs(samples)
  actual_scorez <- Convert_actual_to_win_draw_loss_vector(actual)
  return(calculate_arps(pred_probz[,.( p_win, p_draw, p_loss)], actual_scorez[,.(act_win, act_draw, act_loss)]))
}

calculate_rps_wrap <- function(samples, actual)  {
  pred_probz <- extract_pred_probs(samples)
  actual_scorez <- Convert_actual_to_win_draw_loss_vector(actual)
  return(calculate_rps(pred_probz[,.( p_win, p_draw, p_loss)], actual_scorez[,.(act_win, act_draw, act_loss)]))
}