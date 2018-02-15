# Return cached RDS object or otherwise run an R expression
#
# Return a stored RDS object from a file cache, or if no RDS object is found,
# then evaluate an expression an store the result.
#
# @param expr The expression to evaluate.
# @param filename The name of the object to save or find.
# @param cache_dir The directory file path for the cache.
# @param parse Whether to parse the expression as text before evaluating
#   (\code{TRUE}) or to evaluate the expression as is (\code{FALSE}). If
#   \code{NULL} then the expression will be parsed only if it is a character
#   string.
#
# @author J Buros?
#
with_filecache <- function(expr, filename, cache_dir = "Rcache", parse = NULL) {
  ## prepare expr for evaluation
  expr <- substitute(expr)
  if (is.null(parse)) {
    if ("character" %in% class(expr)) {
      parse=TRUE;
      warning("Detected a character expr; consider wrapping in {} instead of quotes.");
    } else {
      parse=FALSE;
    }
  }
  
  cache_file <- file.path(cache_dir, filename)
  if (!dir.exists(cache_dir)) {
    warning('Cache directory does not exist. Creating one.')
    dir.create(cache_dir, recursive=TRUE)
  }
  if (file.exists(cache_file)) {
    try({obj <- readRDS(cache_file)})
    if (!inherits(obj, 'try-error')) {
      return(obj)
    } else {
      warning('Error reading RDS file -- re-executing expr')
    }
  }
  
  ## evaluate expr
  if (parse) {
    obj = eval(parse(text=expr));
  } else {
    obj = eval(expr)
  }
  
  saveRDS(obj, file = cache_file, compress = TRUE)
  return(obj)
}
