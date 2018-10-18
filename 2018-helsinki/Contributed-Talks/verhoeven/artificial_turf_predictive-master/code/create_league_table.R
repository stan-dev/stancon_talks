create_league_table <- function(dt) {
  dt <- dt[, home_points := 0]
  dt <- dt[goal_difference == 0, home_points := 1]
  dt <- dt[goal_difference > 0, home_points := 3]
  
  dt <- dt[, away_points := 0]
  dt <- dt[goal_difference == 0, away_points := 1]
  dt <- dt[goal_difference < 0, away_points := 3]
  
  home_points <- dt[, .(home_points = sum(home_points)), .(HomeTeam, home_team_id)]
  away_points <- dt[, .(away_points = sum(away_points)), .(AwayTeam, away_team_id)]
  
  setkey(home_points, HomeTeam, home_team_id)
  setkey(away_points, AwayTeam, away_team_id)
  league_table <- home_points[away_points]
  
  league_table[, total_points := home_points + away_points]
}
