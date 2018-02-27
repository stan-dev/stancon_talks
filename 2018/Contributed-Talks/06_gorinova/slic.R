library(stringr)

tostan <- function(model_code, file){
  
  if(missing(model_code))
  {
    stan_name <- gsub(".slic", ".stan", basename(file))
    chars <- system2('slicstan/slicstan.exe', args=c("--fromfile", file), stdout = TRUE, stderr = 'stderr.txt')
  
    return(paste(chars, sep="", collapse="\n"))
  }
  else
  {
    chars <- system2('slicstan/slicstan.exe', args=c(paste(c("\"",model_code,"\""))), 
                     stdout = TRUE, stderr = 'stderr.txt')
    
    return(paste(chars, sep="", collapse="\n"))
  }
}


slic <- function(slicstan_code, file, data){
  
  stan_code = tostan(slicstan_code, file)
  if(missing(data))
  {
    return(stan(model_code=stan_code))
  }
  else
  {
    return(stan(model_code=stan_code, data=data))
  }
}