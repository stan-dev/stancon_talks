MakeEPLPlotAbility <- function(ability, my_data,  flip_sign = F){
  league_table_EPL <- create_league_table(my_data)
  
  if(flip_sign) sign <- -1
  else sign <- 1
  offense <- data.table(melt(ability))
  
  setnames(offense, "Var2", "team_id")
  setnames(offense, "value", "offense")
  
  setkey(league_table_EPL, home_team_id)
  setkey(offense, team_id)
  
  offense <- offense[league_table_EPL]
  
  # 90% credible intervals
  res <- offense[, .(Q5 = quantile(sign * offense, probs = c(0.05)),
                     Q50 = quantile(sign * offense, probs = c(0.5)),
                     Q95 = quantile(sign * offense, probs = c(0.95))), .(HomeTeam, total_points)]
  
  
  gp <- ggplot(res, aes(x = reorder(HomeTeam, total_points), y = Q50)) + geom_point() +
    geom_errorbar(aes(ymin=Q5, ymax=Q95), colour="black", width=.1) + xlab("Team") + ylab("parameter estimate") +
    coord_flip()
  return(gp)
}