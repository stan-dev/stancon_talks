rm(list = ls())
gc()

modelName <- "SteadyState"

## Adjust directories to your settings.
scriptDir <- getwd()
projectDir <- dirname(scriptDir)
figDir <- file.path(projectDir, "deliv", "figure", modelName)
tabDir <- file.path(projectDir, "deliv", "table", modelName)
modelDir <- file.path(projectDir, "model")
outDir <- file.path(modelDir, modelName)
toolsDir <- file.path("tools")
stanDir <- file.path(projectDir, "cmdstan")
tempDir <- file.path(modelDir, modelName, "temp")

## set your lib path
# .libPaths(...)
library(rstan)
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(parallel)

## tools to use Stan and cmdStan from R,
source(file.path(toolsDir, "stanTools.R"))
# source(file.path(toolsDir, "cmdStanTools.R"))

## some helpful plots and generate data
source("SteadyStateSimulations.R")
data <- read_rdump("SteadyState.data.R")

rstan_options(auto_write = TRUE)
set.seed(11191951) ## not required but assures repeatable results

modelName <- "SteadyState"

# Specify the variables for which you want history and density plots
parametersToPlot <- c("CL", "Q", "VC", "VP", "ka", "sigma")

# Additional variables to monitor
otherRVs <- c("cObsPred")

parameters <- c(parametersToPlot, otherRVs)
parametersToPlot <- c("lp__", parametersToPlot)

# initial estimates
init <- function() {
  list(CL = exp(rnorm(1, log(10), 0.2)),
       Q = exp(rnorm(1, log(20), 0.2)),
       VC = exp(rnorm(1, log(70), 0.2)),
       VP = exp(rnorm(1, log(70), 0.2)),
       ka = exp(rnorm(1, log(1), 0.2)),
       sigma = runif(1, 0.5, 2))
}

nChains <- 4
nPost <- 1000 ## Number of post-warm-up samples per chain after thinning
nBurn <- 1000 ## Number of warm-up samples per chain after thinning
nThin <- 1
nIter <- (nBurn + nPost) * nThin
nBurnin <- nBurn * nThin

fit <- stan(file = file.path(modelDir, paste(modelName, ".stan", sep = "")),
            data = data,
            pars = parameters,
            iter = nIter,
            warmup = nBurnin,
            thin = nThin, 
            init = init,
            chains = nChains,
            cores = min(nChains, parallel::detectCores()))

dir.create(outDir)
save(fit, file = file.path(outDir, paste(modelName, "Fit.Rsave", sep = "")))

mcmcHistory(fit, parametersToPlot)
mcmcDensity(fit, parametersToPlot, byChain = TRUE)
pairs(fit, pars = parametersToPlot)

ptable <- parameterTable(fit, parametersToPlot)
ptable


# data <- read_rdump(file.path(modelDir, paste0(modelName,".data.R")))
data <- data.frame(data$cObs, data$time[data$iObs])
data <- plyr::rename(data, c("data.cObs" = "cObs", "data.time.data.iObs." = "time"))

pred <- as.data.frame(fit, pars = "cObsPred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.05),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.95)) %>%
  bind_cols(data)

p1 <- ggplot(pred, aes(x = time, y = cObs))
p1 <- p1 + geom_point() +
  labs(x = "time (h)", y = "plasma concentration (mg/L)") +
  theme(text = element_text(size = 12), axis.text = element_text(size = 12),
        legend.position = "none", strip.text = element_text(size = 8)) 
p1 + geom_line(aes(x = time, y = median)) +
  geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.25)




