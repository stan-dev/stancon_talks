ppc_coverage_plot <- function(ppc_data){

  ggplot(ppc_data[order(true_val)], aes(x = 1:nrow(ppc), y=true_val))  +
  geom_ribbon(aes(ymin=Q5, ymax=Q95), fill = "gray70") + 
  geom_ribbon(aes(ymin=Q25, ymax=Q75), fill = "gray50") + 
  geom_line(aes(y=Q50), col = "red") + geom_point() + 
  ylab("score difference") + xlab("game")
}

ppc_prep_data <- function(score_diff_rep, true_score_diff) {
  score_diff_rep <- data.table(melt(score_diff_rep))
  
  setnames(score_diff_rep, "Var2", "game_id")
  
  ppc <- score_diff_rep[, .(Q5 = quantile(value, 0.05),
                            Q25 = quantile(value, 0.25),
                            Q50 = quantile(value, 0.5),
                            Q75 = quantile(value, 0.75),
                            Q95 = quantile(value, 0.95)), .(game_id)]
  
  actual_scores <- data.table(game_id = 1:length(true_score_diff), 
                              true_val = true_score_diff)
  
  setkey(actual_scores, game_id)
  setkey(ppc, game_id)
  
  ppc <- actual_scores[ppc]
  
  ppc <- ppc[, exp_val := round(Q50, 0)]
  ppc
}