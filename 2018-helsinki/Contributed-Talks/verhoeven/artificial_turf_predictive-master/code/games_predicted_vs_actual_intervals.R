games_predicted_vs_actual_intervals <- function(actual, samples){
  
  score_diff_rep <- data.table(melt(samples))
  
  setnames(score_diff_rep, "Var2", "game_id")
  
  # discretize continous models of goal difference
  score_diff_rep <- score_diff_rep[, value := round(value, 0)]
  

  res <- score_diff_rep[, .N, .(iterations, value)]
  
  # use cross join function from data.table to join to range of actual outcomes
  
  base <- CJ(value = -7:7, iterations = 1:nrow(samples))
  setkey(res, iterations, value)
  setkey(base, iterations, value)
  
  base <- res[base]
  base[is.na(N), N := 0]
  
  fres <- base[, .(Q5 = quantile(N, 0.05),
                   Q50 = quantile(N, 0.5),
                   Q95 = quantile(N, 0.95)), .(value)]
  
  actual_scores <- data.table(game_id = 1:length(actual), 
                              true_val = actual)
  
  actual_scores <- actual_scores[, .(N_ACT = .N), .(true_val)]
  setnames(actual_scores, "true_val", "value")
  setkey(actual_scores, value)
  
  setkey(fres, value)
  
  fres <- actual_scores[fres]
  fres <- fres[is.na(N_ACT), N_ACT := 0]
  
  gp <- ggplot(fres, aes(x= value, y = N_ACT)) + 
    geom_point(col = "red") + geom_line(col = "red", linetype = 2) + geom_point(aes(y = Q50), shape = 0) +
    geom_linerange(aes(ymin = Q5, ymax = Q95)) + xlab("score difference") + ylab("number of games")
  return(gp)
}

games_predicted_vs_actual_intervals_simple <- function(actual, samples, diff_limits = c(-7, 7)){
  
  score_diff_rep <- data.table(melt(samples))
  
  setnames(score_diff_rep, "Var2", "game_id")
  
  # discretize continous models of goal difference
  score_diff_rep <- score_diff_rep[, value := round(value, 0)]

  const <- nrow(score_diff_rep)/length(unique(score_diff_rep$game_id))
  
  res <- score_diff_rep[, .(N = .N/const), (value)]
  
  actual_scores <- data.table(game_id = 1:length(actual), 
                              true_val = actual)
  
  actual_scores <- actual_scores[, .(N_ACT = .N), .(true_val)]
  setnames(actual_scores, "true_val", "value")
  setkey(actual_scores, value)
  
  setkey(res, value)
  
  fres <- actual_scores[res]
  fres <- fres[is.na(N_ACT), N_ACT := 0]
  
  gp <- ggplot(fres[value >= diff_limits[1] & value <= diff_limits[2],], aes(x= value, y = N_ACT, colour = "Actual")) + 
    geom_point(size = 2) + geom_line(linetype = 2) + geom_point(aes(y = N, colour = "Model"), shape = 0, size = 2) +
    geom_line(aes(y = N, colour = "Model"), linetype = 1) +
     xlab("score difference") + ylab("number of games") + scale_colour_manual("", 
                                                                              breaks = c("Actual", "Model"),
                                                                              values = c("red", "black"))
  return(gp)
}