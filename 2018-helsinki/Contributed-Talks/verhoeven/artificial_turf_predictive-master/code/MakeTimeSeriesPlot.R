MakeTimeSeriesPlot <- function(ma_sims, title, id_lut, flip_sign = F, subset = NA) {
  sign <- 1
  if(flip_sign) sign <- -1
  
  res <- ma_sims[, .(Q5 = sign * quantile(value, 0.05),
                     Q50 = sign * quantile(value, 0.5),
                     Q95 = sign * quantile(value, 0.95)), .(nweek, nteam)]
  rm(ma_sims)
  setkey(res, nteam)
  setkey(id_lut, home_team_id)
  
  res <- id_lut[,.(home_team_id, HomeTeam)][res]
  
  if(!is.na(subset[1]))  {
    res <- res[HomeTeam %in% subset,]
  }
  gp <- ggplot(res, aes(x=nweek, y=Q50)) + #geom_point() + 
    geom_line() + facet_wrap(~ HomeTeam, ncol = 4) + 
    geom_ribbon(aes(ymin = Q5, ymax = Q95), alpha = 0.3) + ylab("Team ability") + 
    xlab("Week number") + ggtitle(title)
  gp
}

MakeTimeSeriesPlotSim <- function(a_sims, a_ground_truth, id_lut, title, flip_sign = F) {
  ma_sims <- data.table(melt(a_sims))
  setnames(ma_sims, "Var2", "nweek")
  setnames(ma_sims, "Var3", "nteam")
  
  sign <- 1
  if(flip_sign) sign <- -1
  
  res <- ma_sims[, .(Q5 = sign * quantile(value, 0.05),
                     Q50 = sign * quantile(value, 0.5),
                     Q95 = sign * quantile(value, 0.95)), .(nweek, nteam)]
  
  setkey(res, nteam)
  setkey(id_lut, home_team_id)
  
  res <- id_lut[,.(home_team_id, HomeTeam)][res]
  # prep ground truth
  #a_ground_truth <- model_data$a_defense
  ma_ground_truth <- data.table(melt(a_ground_truth))
  setnames(ma_ground_truth, "team_id", "nteam")
  setnames(ma_ground_truth, "variable", "nweek")
  setnames(ma_ground_truth, "value", "true_value")
  
  ma_ground_truth <- ma_ground_truth[, nweek := as.integer(substring(as.character(nweek), 2))]
  ma_ground_truth <- ma_ground_truth[, nteam := as.integer((as.character(nteam)))]
  
  setkey(ma_ground_truth, nteam)
  setkey(id_lut, home_team_id)
  
  ma_ground_truth <- id_lut[,.(home_team_id, HomeTeam)][ma_ground_truth]
  
  gp <- ggplot(res, aes(x=nweek, y=Q50)) + #geom_point() + 
    geom_line() + geom_line(data = ma_ground_truth, aes(y = true_value), col = "red") + facet_wrap(~ HomeTeam, ncol = 6) + 
    geom_ribbon(aes(ymin = Q5, ymax = Q95), alpha = 0.3) + ylab("Team ability") + 
    xlab("Week number") + ggtitle(title)
  gp
}