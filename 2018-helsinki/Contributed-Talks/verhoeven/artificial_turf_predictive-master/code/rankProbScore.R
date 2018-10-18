# https://opisthokonta.net/?p=1333
# http://pena.lt/y/2013/03/21/how-accurate-are-the-ei-football-predictions/

# observed is the POSITION of the actual outcome (i.e. 1,2 or 3)
rankProbScore <- function(predictions, observed){
  ncat <- ncol(predictions)
  npred <- nrow(predictions)
  
  rps <- numeric(npred)
  
  for (rr in 1:npred){
    obsvec <- rep(0, ncat)
    obsvec[observed[rr]] <- 1
    cumulative <- 0
    for (i in 1:ncat){
      cumulative <- cumulative + (sum(predictions[rr,1:i]) - sum(obsvec[1:i]))^2
    }
    rps[rr] <- (1/(ncat-1))*cumulative
  }
  return(rps)
}

# predictions <- rbind(c(1, 0, 0),
#                         c(0.9, 0.1, 0),
#                         c(0.8, 0.1, 0.1),
#                         c(0.5, 0.25, 0.25))
# observed <- c(1, 1, 1, 1)
# 
# rankProbScore(predictions, observed)
