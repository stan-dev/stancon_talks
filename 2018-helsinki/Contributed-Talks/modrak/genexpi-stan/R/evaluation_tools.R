missing_arg <- function() quote(expr=)

evaluate_single_param_indices <- function(samples, param_name, indices, true_value) {
  if(is.null(indices)) {
    param_samples = samples[[param_name]];
  }
  else {
    #A magic form to run samples[[param_name]][,indices[1], ... , indices[N]] based on the length of indices
    param_samples = do.call(`[`,append(list(samples[[param_name]],missing_arg()),indices))
  }

  if(length(indices) > 0) {
    indices_str = do.call(paste, append(indices, list(sep = ",")))
    fullName =   paste0(param_name,"[", indices_str, "]")
  } else {
    fullName = param_name
  }

  mad_val = mad(param_samples, center = true_value)
  rmse_val = sqrt(mean((param_samples - true_value) ^ 2))
  return(data.frame(
    param_name = fullName,
    true_value = true_value,
    median = median(param_samples),
    IQR = IQR(param_samples),
    quantile = ecdf(param_samples)(true_value),
    mad = mad_val,
    relative_mad = mad_val / true_value,
    relative_rmse = rmse_val / true_value
  ))
}

evaluate_single_param <- function(samples, param_name, param_values)
{
  result = list();
  dimensions <- dim(samples[[param_name]])[-1] #The first dimension is the number of samples
  num_dimensions <- length(dimensions)
  next_element = 1
  if(num_dimensions == 0) {
    result[[next_element]] = evaluate_single_param_indices(samples, param_name, NULL, param_values)
    next_element = next_element + 1
  } else if (num_dimensions == 1) {
    for(i in 1:dimensions[1]) {
      result[[next_element]] = evaluate_single_param_indices(samples, param_name, list(i), param_values[i])
      next_element = next_element + 1
    }
  }
  else if(num_dimensions == 2) {
      for(i in 1:dimensions[1]) {
        for(j in 1:dimensions[2]) {
          result[[next_element]] = evaluate_single_param_indices(samples, param_name, list(i,j), param_values[i,j])
          next_element = next_element + 1
        }
      }
  } else {
    stop("3+ dimensional parameters not supported yet");
  }
  return(do.call(rbind.data.frame, result))
}

evaluate_all_params <- function(samples, true_params) {
  result = list();
  next_element = 1;
  for(param_name in names(true_params)) {
    if(!param_name %in% names(samples)) {
      next;
    }
    param_values = get(param_name, true_params);
    result[[next_element]] = evaluate_single_param(samples, param_name, param_values)
    next_element = next_element + 1
  }
  return(do.call(rbind.data.frame, result));
}

evaluation_summary <- function(samples, true_params, printParamsResults = TRUE) {
  eval_result = evaluate_all_params(samples, true_params);
  if(printParamsResults) {
    print(eval_result);
  }
  quantiles = eval_result$quantile;
  within25 = mean(quantiles >= 0.375 & quantiles <= 0.625);
  within50 = mean(quantiles >= 0.25 & quantiles <= 0.75);
  within95 = mean(quantiles >= 0.025 & quantiles <= 0.975);
  cat("\nWithin 25% interval:", within25,"\nWithin 50% interval:", within50, "\nWithin 95% interval:",within95,"\n")
}

is_vector_like <- function(xs) {
  is.null(dim(xs)) ||
    length(dim(xs)) == 1 ||
    (sum(dim(xs) == 1) >= length(dim(xs)) - 1) #multiple dimensions, but all except one are of size 1
}

ggmatplot <- function(xs, ys, x_title = "x", y_title = "value", main_geom = geom_line()) {
  x_is_vector = is_vector_like(xs)
  y_is_vector = is_vector_like(ys)

  if(y_is_vector) {
    if(x_is_vector) {
      num_rows = length(xs)
      num_cols = 1
      if(num_rows !=  length(ys)) {
        stop("Incompatible vector lengths")
      }
      values = ys
    } else {
      num_rows = dim(xs)[1]
      num_cols = dim(xs)[2]
      if(num_rows != length(ys)) {
        stop("Incompatible y vector with x matrix")
      }
      values = rep(ys, times = num_cols)
    }
  }
  else {
    num_rows = dim(ys)[1]
    num_cols = dim(ys)[2]
    if(x_is_vector) {
      if(num_rows != length(xs)){
        stop("Incompatible y matrix with x vector")
      }
    } else {
      if(num_rows != dim(xs)[1] || num_cols != dim(ys)[2]){
        stop("Incompatible y matrix with x matrix")
      }
    }
    values = as.numeric(ys)
  }

  num_elements = num_rows * num_cols
  data <- data.frame(
    x = array(xs, num_elements),
    column = as.factor(rep(seq(num_cols), each=num_rows)),
    value = values,
    stringsAsFactors=FALSE
  )

  if(num_cols > 10) {
    color_scale = scale_color_discrete()
  } else {
    color_scale = scale_color_brewer(type="qual", palette="Paired")
  }

  ggplot(data, aes(x = x, y = value, color = column)) + main_geom +
    color_scale +
    scale_x_continuous(name = x_title) +
    scale_y_continuous(name = y_title) +
    guides(color = FALSE)
}

averageSamplingTime <- function(fits)
{
  timeList = lapply(fits, get_elapsed_time)
  allTimes = Reduce(rbind,timeList, array(0,c(0,2)))
  warmupTimes = allTimes[,"warmup"]
  sampleTimes = allTimes[,"sample"]
  return(list(total = mean(warmupTimes + sampleTimes), sample = mean(sampleTimes)))
}


launch_shinystan_nonblocking <- function(fit) {
  library(future)
  plan(multisession)
  future(launch_shinystan(fit))
}
