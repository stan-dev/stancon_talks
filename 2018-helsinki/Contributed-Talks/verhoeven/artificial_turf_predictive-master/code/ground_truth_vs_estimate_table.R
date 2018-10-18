ground_truth_vs_estimate_table <- function(samples, true_par_vec, txt){
  
  
  samples <- data.table(melt(samples))
  setnames(samples, "Var2", "team_id")
  
  # 90% credible intervals
  res <- samples[, .(Q5 = quantile(value, probs = c(0.05)),
                     Q50 = quantile(value, probs = c(0.5)),
                     Q95 = quantile(value, probs = c(0.95))), .(team_id)]
  
  res <- res[, true_par_vec := true_par_vec]
  res
}