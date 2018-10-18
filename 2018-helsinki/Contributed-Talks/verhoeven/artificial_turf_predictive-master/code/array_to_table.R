# Author: Mark Klik

array_to_data_table <- function(multi_array) {
  
  # get dimensions
  dims <- attr(multi_array, "dim")
  
  # convert to vector (strip attribute)
  attr(multi_array, "dim") <- NULL
  
  # put data in first column
  dt <- data.table(value = multi_array)
  
  # safe up to 26 dimensions
  col_names <- LETTERS[1:length(dims)]
  
  # total number of rows
  nr_of_rows <- 1
  sapply(dims, function(x) {
    nr_of_rows <<- nr_of_rows * x
  })
  
  count_vec   <- 1:nr_of_rows
  multiplyer  <- 1
  multiplyers <- dims
  
  for (dim in 1:length(dims)) {
    multiplyers[dim] <- multiplyer
    
    dt[, c(col_names[dim]) := 1 + (as.integer((count_vec - 1) / multiplyer)) %% dims[dim]]
    
    multiplyer <- multiplyer * dims[dim]
  }
  
  dt
}

