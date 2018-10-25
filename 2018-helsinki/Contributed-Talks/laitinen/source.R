# Source for the objects needed in the Rmd.

#### Parameters ####
sigma <- 0.1
lambda <- 0.1
mu <- 5
t.df <- 7

chains = 2
iter = 4000

fix_kappa_log <- log(0.1)
fix_mu <- 5

#### Stan model ####
fixed_par_model <- stan_model("fixed_parameter_hierarchical_noncentered.stan")

#### Example plots ####

# lambda = 0.01, 0.1, 1
oup_example_plot <- sapply(c(0.01, 0.1, 1), function(x) generate_a_series(sigma = sigma, 
                                                                          lambda = x, mu = mu, intervals = 1:200, t.df = Inf, 
                                                                          seed = 1)[["observations"]]) %>% as_tibble() %>% 
  set_colnames(c("lambda = 0.01", "lambda = 0.1", "lambda = 1")) %>% 
  mutate(time = 1:200) %>% melt(id.vars = "time") %>% 
  ggplot(aes(x = time, y = value, color = variable)) + 
  geom_hline(yintercept = 5, linetype = "dashed") + 
  geom_line() + scale_y_continuous(limits = c(4, 7)) + 
  theme_bw() + labs(x = "Time", y = "Value", title = "mu = 5, sigma = 0.1") + 
  scale_color_manual(values = c("#4daf4a", "#377eb8", 
                                "#e41a1c")) + theme(legend.title = element_blank())

# sigma = 0.01, 0.1, 1
oup_example_plot2 <- sapply(c(0.05, 0.1, 0.2), function(x) generate_a_series(sigma = x, 
                                                                             lambda = lambda, mu = mu, intervals = 1:200, t.df = Inf, 
                                                                             seed = 1)[["observations"]]) %>% as_tibble() %>% 
  set_colnames(c("sigma = 0.05", "sigma = 0.1", "sigma = 0.2")) %>% 
  mutate(time = 1:200) %>% melt(id.vars = "time") %>% 
  ggplot(aes(x = time, y = value, color = variable)) + 
  geom_hline(yintercept = 5, linetype = "dashed") + 
  geom_line() + theme_bw() + labs(x = "Time", y = "Value", 
                                  title = "mu = 5, lambda = 0.1") + scale_color_manual(values = c("#4daf4a", 
                                                                                                  "#377eb8", "#e41a1c")) + theme(legend.title = element_blank())

#### Single series ####

## data: lambda = c(.1, .3, .5, .7)
single_series_set <- lapply(c(0.1, 0.3, 0.5, 0.7), function(l) lapply(seq(from = 5, 
                                                                          to = 100, by = 5), function(x) {
                                                                            s <- generate_n_series(n = 1, intervals = 1:x, mu = mu, 
                                                                                                   lambda = l, sigma = sigma, fix_mu = 5, fix_kappa_log = log(0.1))
                                                                            return(s)
                                                                          }) %>% set_names(seq(from = 5, to = 100, by = 5)))
names(single_series_set) <- as.character(c(0.1, 0.3, 
                                           0.5, 0.7))


## samples single_series_set_samples <- list() for(i
## in names(single_series_set)) {
## single_series_set_samples[[i]] <-
## lapply(single_series_set[[i]], function(x)
## sampling(fixed_par_model, x, chains=chains,
## iter=iter)) %>% set_names(as.character(seq(from=5,
## to=100, by = 5))) } save(single_series_set_samples,
## file='single_series_set_samples')

load(file = "single_series_set_samples")

## plots

single_series_plots <- lapply(single_series_set_samples, 
                              function(x) {
                                x %>% samples_df() %>% plot_posteriors(sim_value = 0.1, 
                                                                       par = "Lambda ")
                              }) %>% set_names(c(0.1, 0.3, 0.5, 0.7))



#### Two series ####

## Data
two_series_set <- lapply(seq(from = 5, to = 100, by = 5), 
                         function(x) {
                           generate_n_series(n = 2, intervals = 1:x, sigma = sigma, 
                                             lambda = lambda, mu = mu, t.df = t.df, fix_mu = mu, 
                                             fix_kappa_log = fix_kappa_log) %>% concatenate_series()
                         }) %>% set_names(as.character(seq(from = 5, to = 100, 
                                                           by = 5)))

## Samples two_series_set_samples <-
## lapply(two_series_set, function(x)
## sampling(fixed_par_model, x, chains=chains,
## iter=iter)) %>% set_names(as.character(seq(from=5,
## to=100, by = 5))) save(two_series_set_samples,
## file='two_series_set_samples')

load(file = "two_series_set_samples")


## Data frame for results
two_series_results <- two_series_set_samples %>% samples_df2()

## plot
two_series_plot <- two_series_results %>% plot_hier_posteriors(sim_value = lambda) + 
  ggtitle(" ")

## plot comparison
two_series_comparison_plot <- ggplot(rbind(single_series_set_samples[["0.1"]] %>% 
                                             samples_df() %>% mutate(Series = "single"), two_series_results), 
                                     aes(x = Observations, y = mean, group = Series, color = Series)) + 
  geom_errorbar(aes(ymin = lower25, ymax = upper75), 
                width = 1) + geom_line() + geom_point() + geom_hline(yintercept = lambda, 
                                                                     linetype = "dashed", color = "black") + labs(y = "Estimate", 
                                                                                                                  x = "Length") + theme_bw() + scale_color_manual(values = c("#f1a340", 
                                                                                                                                                                             "#998ec3", "black"))


#### Many short series ####

## Data number of observations and series
n_series <- c(2, 3, 5, 10, 20)
n_observations <- c(5, 10, 20, 40, 80)


short_series_set <- list()
for (i in n_observations) {
  
  series <- lapply(n_series, function(x) generate_n_series(n = x, 
                                                           sigma = sigma, mu = mu, lambda = 0.1, t.df = t.df, 
                                                           intervals = 1:i, seed = 1, fix_mu = fix_mu, fix_kappa_log = fix_kappa_log) %>% 
                     concatenate_series()) %>% set_names(as.character(n_series))
  
  series[[as.character(1)]] <- generate_n_series(n = 1, 
                                                 sigma = sigma, mu = mu, lambda = 0.1, t.df = t.df, 
                                                 intervals = 1:i, seed = 1, fix_mu = fix_mu, fix_kappa_log = fix_kappa_log)
  
  short_series_set[[as.character(i)]] <- series
}



## Samples

# short_series_samples <- list() for(i in
# n_observations) { samples <- list() for(j in c(1,
# n_series)) { samples[[as.character(j)]] <-
# sampling(fixed_par_model,
# short_series_set[[as.character(i)]][[as.character(j)]],
# chains=chains, iter=iter) }
# short_series_samples[[as.character(i)]] <- samples
# }

load(file = "short_series_samples")


## Results

short_series_results_list <- lapply(n_observations, function(x) samples_df_par2(short_series_samples[[as.character(x)]], 
                                                                                par = "lambda") %>% mutate(obs = x)) %>% set_names(as.character(c(n_observations)))

## Plot
short_series_plots <- do.call(rbind, short_series_results_list) %>% 
  ggplot(aes(x = length, y = mean)) + geom_point() + 
  geom_hline(yintercept = 0.1, linetype = "dashed") + 
  facet_grid(cols = vars(obs)) + theme_bw() + labs(y = "Posterior mean", 
                                                   x = "Number of series") + scale_x_continuous(breaks = c(1, 
                                                                                                           n_series))

#### Running times ####

# get times form fit list
running_times <- list()

for (i in names(short_series_samples)) {
  t_temp <- list()
  for (j in names(short_series_samples[[i]])) {
    t_temp[[j]] <- sapply(short_series_samples[[i]][[j]]@sim[[1]], 
                          function(x) attr(x, "elapsed_time")) %>% 
      sum
  }
  running_times[[i]] <- t_temp
}

# make data frame
for (i in names(running_times)) {
  running_times[[i]] <- unlist(running_times[[i]]) %>% 
    as_tibble() %>% mutate(n_series = names(running_times[[i]]), 
                           n_obs = i)
}
running_times <- do.call(rbind, running_times) %>% as_tibble

running_times$n_obs <- factor(running_times$n_obs, levels = as.character(c(5, 
                                                                           10, 20, 40, 80)))
running_times$n_series <- factor(running_times$n_series, 
                                 levels = unique(running_times$n_series))


# plot

running_times_plot <- running_times %>% ggplot(aes(x = as.numeric(as.character(n_series)), 
                                                   y = value/60, color = n_obs)) + geom_point() + labs(y = "Minutes", 
                                                                                                       x = "Number of series") + guides(color = guide_legend(title = "Observations")) + 
  theme_bw() + geom_line()

#### Diagnostics