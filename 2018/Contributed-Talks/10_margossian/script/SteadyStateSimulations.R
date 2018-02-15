# rm(list = ls())
# gc()
set.seed(11091989) ## not required but assures repeatable results

library(ggplot2)
library(mrgsolve)  ## tested with version 0.8.6

## mrgsolve model
code <- '
$PARAM CL = 5, Q = 8, V2 = 20, V3 = 70, KA = 1.2

$CMT GUT CENT PERI

$GLOBAL
#define CP (CENT/V2)
 
$PKMODEL ncmt = 2, depot = TRUE

$SIGMA 0.001  // variance

$TABLE
capture DV = CP * exp(EPS(1));

$CAPTURE CP
'

mod <- mcode("accum", code) %>% Req(DV) %>% update(end=480, delta=0.1)

## Case 1: a patient starts a treatment and experiences changes
e1 <- ev(amt = 1000, ii = 8, addl = 9)
# mod %>% ev(e1) %>% mrgsim(end = 80) %>% plot
time <- seq(from = 0.01, to = 80, by = 0.01)
data <- mod %>% ev(e1) %>% 
        carry_out(cmt, ii, addl, amt, evid) %>%
        mrgsim(Req = "DV", end = -1, add = time, recsort = 3) %>%
        as.data.frame
data$Concentration <- data$DV
plot1 <- ggplot(data = data, aes(time, Concentration)) + geom_line(color = 'blue') + theme_bw()
plot_cObs1 <- plot1 + ylim(0, 45)

## Case 2: patient starts at a steady state
e2 <- ev(amt = 1000, ii = 8, addl = 9, ss = 1) # Create dosing events
# mod %>% ev(e2) %>% mrgsim(end = 80) %>% plot # plot data
data <- mod %>% ev(e2) %>% 
  carry_out(cmt, ii, addl, amt, evid) %>%
  mrgsim(Req = "DV", end = -1, add = time, recsort = 3) %>%
  as.data.frame
data$Concentration <- data$DV
plot2 <- ggplot(data = data, aes(time, Concentration)) + geom_line(color = 'blue') + theme_bw()
plot_cObs2 <- plot2 + ylim(0, 45)
