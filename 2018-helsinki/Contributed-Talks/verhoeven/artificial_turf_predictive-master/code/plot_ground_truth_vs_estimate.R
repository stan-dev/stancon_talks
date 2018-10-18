plot_ground_truth_vs_estimate <- function(samples, true_par_vec, txt, art_turf_vec = NA, dotsize = 1){

  if(is.na(art_turf_vec)[1]) art_turf_vec <- rep(0, length(true_par_vec))
  samples <- data.table(melt(samples))
  setnames(samples, "Var2", "team_id")
  
  # 90% credible intervals
  res <- samples[, .(Q5 = quantile(value, probs = c(0.05)),
                     Q50 = quantile(value, probs = c(0.5)),
                     Q95 = quantile(value, probs = c(0.95))), .(team_id)]
  
  res <- res[, true_par_vec := true_par_vec]
  art_turf <- data.table(team_id = 1:length(art_turf_vec), art_turf = art_turf_vec)
  setkey(art_turf, team_id)
  setkey(res, team_id)
  res <- art_turf[res]
  
  gp <- ggplot(res, aes(x = true_par_vec, y = Q50)) + 
    geom_abline(slope = 1, intercept = 0) + 
    ggtitle(txt) + geom_linerange(aes(ymin = Q5, ymax = Q95)) + 
    geom_point(aes(col = factor(art_turf), shape = factor(art_turf)), size = dotsize) + 
    geom_text(aes(label = team_id, hjust = 1, vjust = 1), col = "red")
  gp
}

