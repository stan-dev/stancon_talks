source("../code/rankProbScore.R")
source("../code/ppc_coverage_plot.R")
source("../code/MakeTimeSeriesPlot.R")
source("../code/Create_model_data_for_TS2.R")
source("../code/addTeamIds.R")
source("../code/create_league_table.R")
source("../code/MakeEPLPlotAbility.R")
source("../code/games_predicted_vs_actual_intervals.R")
source("../code/ppc_coverage_plot.R")
source("../code/calc_rps_scores.R")
source("../code/odds_to_probability.R")
source("../code/ReadfitsandCalculateRPS.R")
source("../code/array_to_table.R")
source("../code/FitStanModel.R")

NL_ALL <- readRDS("../data/NL Eredivisie 2000-2018.rds")

# set 2017/2018 season apart
NL_17 <- NL_ALL[Date > as.Date("2017-07-01")]

NL_ALL <- NL_ALL[Date < as.Date("2017-07-01")]
setkey(NL_ALL, Date)

# add round and season
nrounds <- nrow(NL_ALL)/9
nseasons <- nrow(NL_ALL)/(9*34)
NL_ALL <- NL_ALL[, round := rep(1:nrounds, each = 9)]
NL_ALL <- NL_ALL[, season := rep(1:nseasons, each = 34*9)]

# prep 2017/2018 separately
setkey(NL_17, Date)
nrounds <- ceiling(nrow(NL_17)/9)
start_nr_round <- max(NL_ALL$round)+1

round_vec <- rep(start_nr_round:(start_nr_round + nrounds), each = 9)
NL_17 <- NL_17[, round := round_vec[1:nrow(NL_17)]]
NL_17 <- NL_17[, season := 18]

# add to NL_ALL
NL_ALL <- rbind(NL_ALL, NL_17)

setkey(NL_ALL, Date)

NL_ALL <- NL_ALL[, row_id := 1:nrow(NL_ALL)]
#saveRDS(NL_ALL, "data/NL_ALL.rds")

# select only 2016/2017 season
NL <- NL_ALL[season == 17,]
setkey(NL, Date)
NL <- addTeamIds(NL)
# deze heeft al stan team ids 1-18

# select 2014/2015, 2015/2016 and 2016/2017 season
NLthree <- NL_ALL[season %in% 15:17]

# make sure team_ids number continuously from 1 to nteam
NLthree <- addTeamIds(NLthree)

NLsix <- NL_ALL[season %in% 12:17]

NLsix <- addTeamIds(NLsix)

# sort data again to be sure
setkey(NL, Date)
setkey(NLthree, Date)
setkey(NLsix, Date)
