default_expression_plot_num_samples <- 100
default_expression_plot_main_geom <- geom_line(alpha = 0.2)


get_total_samples <- function(fit) {
  sapply(fit@stan_args, FUN = function(x) {x$iter - x$warmup}) %>% sum()
}

fitted_regulator_plot <- function(fit, data, regulator = 1, name = NULL,
                                  num_samples = default_expression_plot_num_samples,
                                  main_geom = default_expression_plot_main_geom) {
  samples_to_show <-  sample(1:get_total_samples(fit), num_samples)
  samples_regulator <- rstan::extract(fit,"predicted_regulator_expression")$predicted_regulator_expression[samples_to_show,,regulator]

  if(!is.null(name)) {
    regulator_name <- name
  } else if(is.null(colnames(data$regulator_expression))) {
    regulator_name <- regulator
  } else {
    regulator_name <- colnames(data$regulator_expression)[regulator]
  }

  regulator_plot <- ggmatplot(1:data$num_time, t(samples_regulator), main_geom = main_geom,
                              x_title = "time", y_title = "expression") +
    ggtitle(paste0("Regulator - ", regulator_name))

  if(data$regulators_measured != 0) {
    regulator_plot <- regulator_plot +
      geom_point(data = data.frame(x = data$measurement_times,  y = data$regulator_expression[,regulator]),
                 aes(x=x, y=y), inherit.aes = FALSE, color = "#ba1b1d", size = 3)
  }

  regulator_plot
}

fitted_target_plot <- function(fit, data, target = 1, name = NULL,
                               num_samples = default_expression_plot_num_samples,
                               main_geom = default_expression_plot_main_geom) {

  samples_to_show <-  sample(1:get_total_samples(fit), num_samples)

  samples_expression <-
    rstan::extract(fit,"predicted_expression")$predicted_expression[samples_to_show,,target,drop=FALSE]

  if(!is.null(name)) {
    target_name <- name
  } else if(is.null(colnames(data$expression))) {
    target_name <- target
  } else {
    target_name <- colnames(data$expression)[target]
  }


  ggmatplot(1:data$num_time, t(samples_expression[,,1]),
            main_geom = main_geom, x_title = "time", y_title = "expression") +
          geom_point(
            data = data.frame(x = data$measurement_times,  y = data$expression[,target]),
            aes(x=x, y=y), inherit.aes = FALSE, color = "#ba1b1d", size = 3) +
    ggtitle(paste0("Expression - ", target_name))

}

fitted_target_observed_plot <- function(fit, data, target = 1, name = NULL,
                               num_samples = default_expression_plot_num_samples,
                               main_geom = default_expression_plot_main_geom) {

  samples_to_show <-  sample(1:get_total_samples(fit), num_samples)

  samples_expression <-
    rstan::extract(fit,"expression_replicates")$expression_replicates[samples_to_show,,target,drop=FALSE]

  if(!is.null(name)) {
    target_name <- name
  } else if(is.null(colnames(data$expression))) {
    target_name <- target
  } else {
    target_name <- colnames(data$expression)[target]
  }


  ggmatplot(data$measurement_times, t(samples_expression[,,1]),
            main_geom = main_geom, x_title = "time", y_title = "observed expression") +
    geom_point(
      data = data.frame(x = data$measurement_times,  y = data$expression[,target]),
      aes(x=x, y=y), inherit.aes = FALSE, color = "#ba1b1d", size = 3) +
    ggtitle(paste0("Expression obs. - ", target_name))

}

fitted_csynth_plot <- function(fit, data, name = NULL,
                               num_samples = default_expression_plot_num_samples,
                               main_geom = default_expression_plot_main_geom) {

  samples_to_show <-  sample(1:get_total_samples(fit), num_samples)

  samples_expression <-
    rstan::extract(fit,"predicted_expression")$predicted_expression[samples_to_show,,drop=FALSE]

  if(!is.null(name)) {
    target_name <- name
  } else {
    target_name <- "N/A"
  }


  ggmatplot(data$measurement_times, t(samples_expression),
            main_geom = main_geom, x_title = "time", y_title = "expression") +
    geom_point(
      data = data.frame(x = data$measurement_times,  y = data$expression),
      aes(x=x, y=y), inherit.aes = FALSE, color = "#ba1b1d", size = 3) +
    ggtitle(paste0("Constant synthesis - ", target_name))

}

fitted_csynth_observed_plot <- function(fit, data, target = 1, name = NULL,
                                        num_samples = default_expression_plot_num_samples,
                                        main_geom = default_expression_plot_main_geom) {

  samples_to_show <-  sample(1:get_total_samples(fit), num_samples)

  samples_expression <-
    rstan::extract(fit,"expression_replicates")$expression_replicates[samples_to_show,,drop=FALSE]

  if(!is.null(name)) {
    target_name <- name
  } else {
    target_name <- "N/A"
  }


  ggmatplot(data$measurement_times, t(samples_expression),
            main_geom = main_geom, x_title = "time", y_title = "observed expression") +
    geom_point(
      data = data.frame(x = data$measurement_times,  y = data$expression),
      aes(x=x, y=y), inherit.aes = FALSE, color = "#ba1b1d", size = 3) +
    ggtitle(paste0("Constant synthesis obs. - ", target_name))

}
