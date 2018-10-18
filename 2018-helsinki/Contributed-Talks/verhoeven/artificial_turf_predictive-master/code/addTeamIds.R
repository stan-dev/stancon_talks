addTeamIds <- function(dt) {
  # add team id's for Stan
  TeamList <- unique(c(dt$HomeTeam, as.character(dt$AwayTeam)))
  dt <- dt[, home_team_id := match(HomeTeam, TeamList)]
  dt <- dt[, away_team_id := match(AwayTeam, TeamList)]
}
