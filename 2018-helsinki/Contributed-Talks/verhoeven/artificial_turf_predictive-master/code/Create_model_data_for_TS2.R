# requires a dataset for starting values, 
# a dataset for parameter estimation and 
# a dataset for prediction
# NOTE always include at least two matches for prediction (stan model requires this in its current form)
Create_model_data_for_TS2 <- function(my_data, my_data_prev = NA, my_data_pred = NA){
  setkey(my_data, row_id)
  
  home_week <- c() # vector for each game in which week it was played for the home team
  away_week <- c() # vector for each game in which week it was played for the away team
  ngames = nrow(my_data)
  nteams = length(unique(c(my_data$HomeTeam, as.character(my_data$AwayTeam))))

  
  # make sure team_ids number continuously from 1 to nteam
  my_data <- addTeamIds(my_data)
  
  # adapted from Milad
  for(g in 1:ngames) {
    home_week[g] <- sum(my_data$home_team_id[1:g] == my_data$home_team_id[g]) + # vector of 1/0
      sum(my_data$away_team_id[1:g] == my_data$home_team_id[g])
    away_week[g] <- sum(my_data$away_team_id[1:g] == my_data$away_team_id[g]) + 
      sum(my_data$home_team_id[1:g] == my_data$away_team_id[g])
  }
  # prev_performance block
  if(!is.na(my_data_prev[1,1])){
    league_table_prev_perf <- create_league_table(my_data_prev)
    
    # make sure the ids are matched
    ids_my_data <- my_data[, .N, .(home_team_id, HomeTeam)]
    # hack for data set where a club is only present in AwayTeam
    if(length(unique(my_data$HomeTeam)) != length(unique(my_data$AwayTeam))){
      ids_my_data_away_only <- my_data[!(AwayTeam %in% ids_my_data$HomeTeam), 
                                  .N, .(away_team_id, AwayTeam)]
      setnames(ids_my_data_away_only, "AwayTeam", "HomeTeam")
      setnames(ids_my_data_away_only, "away_team_id", "home_team_id")
      ids_my_data <- rbind(ids_my_data, 
                           ids_my_data_away_only)
    }
    
    setkey(league_table_prev_perf, HomeTeam)
    setkey(ids_my_data, HomeTeam)
    
    ids_my_data <- league_table_prev_perf[ids_my_data]
    
    # map points to score between -1 and 1
    
    # for teams that did not play in prev season set starting performance at -1 
    # the logic here is that since earlier they fell out of the league, 
    # so they are presumably not that good
    
    ids_my_data <- ids_my_data[, prev_perf := -1]
    ids_my_data <- ids_my_data[!is.na(total_points), prev_perf := map_to_score(total_points)]
    
    prev_perf <- ids_my_data[order(i.home_team_id),]$prev_perf
  }
  id_lut <- ids_my_data[,.(i.home_team_id, HomeTeam)] # MakeTimeSeriesPlot() needs this
  setnames(id_lut, "i.home_team_id", "home_team_id" )
  
  # OUT OF SAMPLE PREDICTION PART
  # replace the id's in my_data_pred with the ids in my_data
  # drop matches with teams that are not included in my_data
  
  # make sure the ids are matched, for homeTeam as well as awayTeam
  ids_my_data_home <- unique(my_data[, .(home_team_id, HomeTeam)])
  ids_my_data_away <- unique(my_data[, .(away_team_id, AwayTeam)])
  
  
  my_data_pred <- my_data_pred[, home_team_id := NULL]
  my_data_pred <- my_data_pred[, away_team_id := NULL]
  
  setkey(ids_my_data_home, HomeTeam)
  setkey(my_data_pred, HomeTeam)
  my_data_pred <- ids_my_data_home[my_data_pred]
  
  setkey(ids_my_data_away, AwayTeam)
  setkey(my_data_pred, AwayTeam)
  my_data_pred <- ids_my_data_away[my_data_pred]
  
  # drop matches without team_ids
  my_data_pred <- na.omit(my_data_pred)
  # count matches to predict
  ngames_pred = nrow(my_data_pred)

  model_data <- list(
    n_teams = nteams,
    n_games = ngames,
    n_weeks = (ngames*2)/(18),
    home_week = home_week,
    away_week = away_week,
    home_team = my_data$home_team_id,
    away_team = my_data$away_team_id,
    goal_difference = my_data$goal_difference,
    prev_perf = prev_perf,
    n_games_pred = nrow(my_data_pred),
    home_team_pred = my_data_pred$home_team_id,
    away_team_pred = my_data_pred$away_team_id,
    pred_goal_difference = my_data_pred$goal_difference,
    art_turf = my_data$art_turf_advantage,
    id_lut = id_lut
  )
  return(model_data)
}

# map league points to score between -1 and 1
map_to_score <- function(x){
  x_max <- max(x);
  x_min <- min(x);
  return(2*x/(x_max - x_min) - (x_max + x_min) / (x_max - x_min))
}
