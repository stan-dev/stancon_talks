samples_df_par2 <- function(sample_list, par) {
  
  res_df <- matrix(NA, 1, 4)
  
  colnames(res_df) <- c("length", "lower25", "mean", "upper75")
  
  
  for(i in sample_list) {
    
    s <- summary(i)$summary
    
    res <- s[grepl(paste0(par, "\\["), rownames(s)),c("25%", "50%", "75%")]
    
    if(class(res) == "numeric") {
      res_df <- res_df %>% rbind(c(1, res))
    } else {
      res_df <- res_df %>% rbind(cbind(nrow(res), res))
    }
    
  }
  
  res_df <- res_df[-1,]
  res_df <- res_df %>% as_tibble()
  
  res_df
  
}


plot_hier_posteriors <- function(res_df, sim_value) {
  p <- ggplot(res_df, aes(x=Observations, y=mean, group=Series, color=Series)) + 
    geom_errorbar(aes(ymin=lower25, ymax=upper75), width=1) +
    geom_line() +
    geom_point()  + geom_hline(yintercept = sim_value, linetype="dashed", color="black") + labs(y="Estimate", x="Length", title="Lambda estimates with 50% error bars vs. series length") + theme_bw() + scale_color_manual(values=c("#f1a340", "#998ec3")) + scale_y_continuous(limits = c(0, max(res_df$upper)))
  
  p
}

samples_df2 <- function(sample_list) {
  
  # make data frame for results
  res_df <- matrix(NA, 2*length(sample_list), 5) %>% set_colnames(c("Series", "Observations", "lower25", "mean", "upper75"))
  
  
  j <- 1
  for(i in sample_list) {
    res <- summary(i)$summary[grep("lambda\\[", rownames(summary(i)$summary)),c("25%", "50%", "75%")]
    
    res_df[c(j, j+1), c("lower25", "mean", "upper75")] <- res
    j <- j + 2
  }
  
  
  res_df[, "Observations"] <- rep(as.numeric(names(sample_list)), each=2)
  res_df <- res_df %>% as.data.frame()
  res_df[, "Series"] <- rep(c("A","B"), length(sample_list))
  
  res_df
} 

concatenate_series <- function(s_list) {
  
  len <- length(s_list)
  
  if(len == 1) {
    return(s_list[[1]])
  }
  
  obs <- c()
  t <- c()
  s <- c()
  for(i in 1:length(s_list)) {
    obs <- c(obs, s_list[[i]][["observations"]])
    t <- c(t, s_list[[i]][["time"]])
    s <- c(s, s_list[[i]][["samples_per_series"]])
  }
  
  kappa_log <- rep(s_list[[1]][["kappa_log"]], len)
  mu <- rep(s_list[[1]][["mu"]], len)
  
  T <- length(obs)
  
  
  list(observations=obs, time=t, n_series=len, samples_per_series=s, T=T, kappa_log=kappa_log, mu=mu)
  
}

samples_df <- function(sample_list) {
  res_df <- matrix(NA, length(sample_list), 4)
  
  colnames(res_df) <- c("Observations", "lower25", "mean", "upper75")
  
  j <- 1
  for(i in sample_list) {
    res <- summary(i)$summary[grep("lambda\\[", rownames(summary(i)$summary)),c("25%", "50%", "75%")]
    
    res_df[j, c("lower25", "mean", "upper75")] <- res
    j <- j + 1
  }
  
  res_df[, "Observations"] <- as.numeric(names(sample_list))
  res_df <- res_df %>% as.data.frame()
  
  res_df
}

plot_posteriors <- function(res_df, sim_value, par) {
  p <- ggplot(res_df, aes(x=Observations, y=mean)) + 
    geom_errorbar(aes(ymin=lower25, ymax=upper75), width=1) +
    geom_line() +
    geom_point() + geom_hline(yintercept = sim_value, linetype="dashed") + labs(y="Estimate", x="Length") + theme_bw() 
  # + scale_y_continuous(limits = c(0, 6.5))
  
  p
} 

generate_a_series <- function(kappa=NULL, sigma=NULL, lambda, mu, intervals, t.df = Inf, seed = 1){
  set.seed(seed)
  N <- length(intervals)
  lv.variates <- rnorm(N)
  dt <- outer(intervals,intervals,function(x,y) abs(x-y))
  
  if(!is.null(kappa)) {
    x <- kappa * exp(-lambda*dt)
  } else {
    x <- (sigma^2)/lambda * exp(-lambda*dt)
  }
  
  
  L <- chol(x)
  
  scale <- if(is.finite(t.df)) rep(sqrt(rgamma(1,t.df/2,(t.df-2)/2)),each=N) else 1
  out.data <- list()
  out.data$observations <- as.vector(t(L) %*% (rnorm(N) * scale)) + mu
  # out.data$observations <- rpois(N,exp(as.vector(t(L) %*% (rnorm(N) * scale))+mu))
  out.data$time <- intervals
  out.data$n_series <- 1L
  out.data$samples_per_series <- array(length(intervals))
  out.data$T <- length(intervals)
  out.data
}

generate_n_series <- function(n, kappa=NULL, sigma=NULL, lambda, mu, intervals, t.df = Inf, seed = 1, fix_mu=NULL, fix_kappa_log=NULL) {
  
  s_list <- list()
  
  for(i in 1:n) {
    s_list[[i]] <- generate_a_series(kappa=kappa, sigma=sigma, lambda=lambda, mu=mu, intervals=intervals, t.df = t.df, seed = seed + i)
    
    if(!is.null(fix_mu)) {
      s_list[[i]][["mu"]] <- as.array(fix_mu)
    }
    if(!is.null(fix_kappa_log)) {
      s_list[[i]][["kappa_log"]] <- as.array(fix_kappa_log)
    }
    
  }
  
  if(length(s_list) == 1) {
    return(s_list[[1]])
  }
  s_list
  
}
