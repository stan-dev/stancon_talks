# Set up -----------------------------------------------------------------------

library(rstan)
library(edstan)
library(loo)
library(reshape2)
library(doParallel)
options(mc.cores = 5)
options(loo.cores = 5)


# Functions --------------------------------------------------------------------

# Replacement for rstan::get_posterior_means() that returns object with same
# structure as rstan::extract()
# stan_fit: A fitted Stan model
better_posterior_means <- function(stan_fit) {
  draws <- extract(stan_fit, stan_fit@model_pars)
  f <- function(x) {
    dims <- dim(x)
    n_dims <- length(dims)
    if (n_dims == 1) {
      mean(x)
    } else {
      m <- apply(x, 2:n_dims, mean)
      array(m, dim = c(1, dims[-1]))
    }
  }
  lapply(draws, f)
}


# Function to obtain marginal likelihoods with parallel processing. 
# stan_fit: Fitted Stan model
# data_list: Data list used in fitting model
# MFUN: Function to calculate marginal likelihood for cluster at a node 
#   location. This is application specific.
# resid_name: Name of residual in Stan program to integrate out
# sd_name: Name of SD for residual in Stan program
# n_nodes: Number of adaptive quadrature nodes to use
# best_only: Whether to evaluate marginal likelihood only at posterior means
mll_parallel <- function(stan_fit, data_list, MFUN, resid_name, sd_name, n_nodes,
                         best_only = FALSE) {
  
  library(foreach)
  library(statmod)       # For gauss.quad.prob()
  library(matrixStats)   # For logSumExp()
  
  draws <- extract(stan_fit, stan_fit@model_pars)
  n_iter <- ifelse(best_only, 0, nrow(draws$lp__))
  post_means <- better_posterior_means(stan_fit)
  
  # Seperate out draws for residuals and their SD
  resid <- apply(draws[[resid_name]], 2, mean)
  stddev <- apply(draws[[resid_name]], 2, sd)
  
  # Get standard quadrature points
  std_quad <- gauss.quad.prob(n_nodes, "normal", mu = 0, sigma = 1)
  std_log_weights <- log(std_quad$weights)
  
  # Extra iteration is to evaluate marginal log-likelihood at parameter means.
  ll <- foreach(i = 1:(n_iter + 1), .combine = rbind,
                .packages = "matrixStats") %dopar% {
                  
    ll_j <- matrix(NA, nrow = 1, ncol = ncol(draws[[resid_name]]))
    
    for(j in 1:ncol(ll_j)) {
      
      # Set up adaptive quadrature using SD for residuals either from draws or
      # posterior mean (for best_ll).
      sd_i <- ifelse(i <= n_iter, draws[[sd_name]][i], post_means[[sd_name]])
      adapt_nodes <- resid[j] + stddev[j] * std_quad$nodes
      log_weights <- log(sqrt(2*pi)) + log(stddev[j]) + std_quad$nodes^2/2 +
        dnorm(adapt_nodes, sd = sd_i, log = TRUE) + std_log_weights
      
      # Evaluate mll with adaptive quadrature. If at n_iter + 1, evaluate
      # marginal likelihood at posterior means.
      if(i <= n_iter) {
        loglik_by_node <- sapply(adapt_nodes, FUN = MFUN, r = j, iter = i,
                                 data_list = data_list, draws = draws)
        weighted_loglik_by_node <- loglik_by_node + log_weights
        ll_j[1,j] <- logSumExp(weighted_loglik_by_node)
      } else {
        loglik_by_node <- sapply(adapt_nodes, FUN = MFUN, r = j, iter = 1,
                                 data_list = data_list, draws = post_means)
        weighted_loglik_by_node <- loglik_by_node + log_weights
        ll_j[1,j] <- logSumExp(weighted_loglik_by_node)
      }
      
    }
    
    ll_j
    
  }
  
  if(best_only) {
    return(ll[nrow(ll), ])
  } else {
    return(list(ll = ll[-nrow(ll), ], best_ll = ll[nrow(ll), ]))
  }
  
}


# Function to calculate likelihood for a cluster for an adaptive quad node
# specific to the IRT example. Similar functions would be written for other
# applications and passed to mll_parallel().
# node: node location
# r: index for cluster
# iter: mcmc iteration
# data_list: data used to fit Stan model
# draws: mcmc draws from fitted Stan model
f_marginal <- function(node, r, iter, data_list, draws) {
  y <- data_list$y[data_list$jj == r]
  theta_fix <- draws$theta_fix[iter, r]
  delta <- draws$delta[iter, data_list$ii[data_list$jj == r]]
  p <- boot::inv.logit(theta_fix + node - delta)
  sum(dbinom(y, 1, p, log = TRUE))
}


# Function to calculate DIC
# ll_obj: Object returned by mll_parallel()
dic <- function(ll_obj) {
  full_ll <- apply(ll_obj$ll, 1, sum)
  full_best <- sum(ll_obj$best_ll)
  mean_lpd <-  mean(full_ll)
  pdic <- 2 * (full_best - mean_lpd)
  elpd_dic <- full_best - pdic
  c(elpd_dic = elpd_dic, p_dic = pdic, dic = -2*elpd_dic,
    best_lpd = full_best, mean_lpd = mean_lpd)
}


# Example analysis -------------------------------------------------------------

# Assemble example dataset
dl <- irt_data(y = aggression$dich, jj = aggression$person,
               ii = aggression$item, covariates = aggression,
               formula = ~ 1 + male + anger)

# Fit model
fit <- stan("rasch_edstan_modified.stan", data = dl, iter = 500, chains = 5)

# Obtain marginal likelihoods
cl <- makeCluster(5)
registerDoParallel(cl)
ll_marg <- mll_parallel(fit, dl, f_marginal, "zeta", "sigma", 11)
stopCluster(cl)

# Obtain marginal information criteria
dic(ll_marg)
waic(ll_marg$ll)
loo(ll_marg$ll)


