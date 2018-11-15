#
# Functions for constructing
# Mixed Effects Bayesian Network
#
# Jari Turkia
#

mebn.load_datadesc <- function(DataDescriptionFile = "Data description.xlsx")
{
  library(xlsx)    
    
  # Read data description
  datadesc <- read.xlsx(DataDescriptionFile,sheetIndex=1,header=TRUE,rowIndex=c(1:50),colIndex=c(1:9))
  
  return(datadesc)
}

##################################################

mebn.new_graph_with_randomvariables <- function(datadesc)
{
  library(igraph)  

  # Initialize graph structure
  reaction_graph <- make_empty_graph() 
  
  # TODO: Get nodes styles from description, and not like this
  min_order <- min(datadesc$Order)
  predictor_columns <- datadesc[datadesc$Order==min_order,]
  
  next_order <- max(datadesc$Order)
  assumed_targets <- datadesc[datadesc$Order==next_order,]
  
  # TODO: Get shape?
  
  # Add nodes to reaction graph for all the random variables  
  reaction_graph <- reaction_graph + vertices(as.vector(assumed_targets$Name), 
                                              label=as.vector(assumed_targets$Description), 
                                              type="rv",
                                              color = "#74aaf2", 
                                              size = 1,
                                              shape = "circle")
                                              
                                              # shape = "distbox",
                                              # mean=50,
                                              # l95CI=30,
                                              # u95CI=70,
                                              # scalemin=0,
                                              # scalemax=100)
  
  reaction_graph <- reaction_graph + vertices(as.vector(predictor_columns$Name), 
                                              label=as.vector(predictor_columns$Description), 
                                              type="rv",
                                              color = "#3cd164", 
                                              size = 1, 
                                              shape = "circle")
  
                                              # shape = "distbox",
                                              # mean=50,
                                              # l95CI=30,
                                              # u95CI=70,
                                              # scalemin=0,
                                              # scalemax=100)
  

  return(reaction_graph)
}

##################################################

mebn.add_priornodes <- function(datadesc, reaction_graph)
{
  require(igraph)  
  
  # Add prior information from datadesc for random variables
  # - filter prior information for the nutrition predictors
  
  predictor_columns <- datadesc[datadesc$Order==100,]
  
  predictor_names <- as.vector(predictor_columns$Name)
  predprior <- datadesc[(!is.na(datadesc$Lowerbound) | !is.na(datadesc$Upperbound)) & datadesc$Order == 100,]
  predprior$PriorName <- paste0(predprior$Name, "_prior")
  
  # - add prior nodes
  reaction_graph <- reaction_graph + vertices(as.vector(predprior$PriorName), 
                                              label=as.vector(predprior$PriorName), 
                                              color = "#92d3a5", 
                                              size = 1, 
                                              shape = "distbox",
                                              mean=as.vector(predprior$Lowerbound+(predprior$Upperbound-predprior$Lowerbound)/2),
                                              l95CI=as.vector(predprior$Lowerbound),
                                              u95CI=as.vector(predprior$Upperbound),
                                              scalemin=as.vector(predprior$ScaleMin),
                                              scalemax=as.vector(predprior$ScaleMax))
  
  # Connect priors to predictors
  for (p in 1:nrow(predprior))
  {
    print(paste0(predprior[p,]$PriorName, " -> ", predprior[p,]$Name))
    reaction_graph <- reaction_graph + edges(as.vector(predprior[p,]$PriorName), as.vector(predprior[p,]$Name), shape = "arrow", weight = 1)
  }
  
  return(reaction_graph)
}

##################################################

mebn.scale_gaussians <- function(r, data, datadesc) 
{ 
  # data parameters contains whole dataset 
  # - subset only field(s) in datadesc so that index r matches
  
  data <- data[as.vector(datadesc$Name)]
  s <- data[,r] # data from column r
  
  if (datadesc[r,]$Distribution == "Gaussian")
  {
    # TODO: Scale also returns the scaling factor. Store it to restore the original scale.
    s <- scale(s, center = FALSE, scale = TRUE)
  }
  
  return (s) 
}

##################################################

mebn.set_model_parameters <- function(predictor_columns, target_column, group_column, inputdata, normalize_values, reg_params = NULL)
{
  normalize <- function(x) { return ((x - min(x)) / (max(x) - min(x))) }
  
  ident <- function(x) { return (x) }
  predictors <- inputdata[as.vector(predictor_columns$Name)]
  
  target_name <- as.vector(target_column$Name)
  
  prior_sigma <- rep(-1, nrow(predictor_columns))
  dim(prior_sigma) <- nrow(predictor_columns)
  prior_mean <- rep(0, nrow(predictor_columns))
  dim(prior_mean) <- nrow(predictor_columns)
  
  # Set informative priors 
  if (!is.na(predictor_columns$Lowerbound) && !is.na(predictor_columns$Upperbound))
  {  
    prior_sigma <- c(predictor_columns$Upperbound - predictor_columns$Lowerbound)
    prior_mean <- c(predictor_columns$Lowerbound + prior_sigma / 2)
    
    dim(prior_sigma) <- length(prior_sigma)
    dim(prior_mean) <- length(prior_mean)
  }
  
  # Scale if the predictor is Gaussian
  N <- nrow(inputdata)
  
  if (normalize_values == TRUE)
  {
    X <- sapply(1:nrow(assumedpredictors), mebn.scale_gaussians, data = sysdimet, datadesc = assumedpredictors)
    
    # append intercept 
    X <- cbind(rep(1,N), X)
    
    # TODO: condscale Y, too
    Y <- scale(inputdata[target_name][,], center = FALSE, scale = TRUE)[,1]
  }
  else
  {
    X <- cbind(rep(1,N), apply(predictors, 2, ident))
    Y <- inputdata[target_name][,]
  }

  params <- within(list(),
                   {
                     N <- N
                     X <- X
                     p <- k <- ncol(X)               # all predictors may have random effects
                     # Mean and variance of Gaussian prior predictors
                     X_prior_sigma <- prior_sigma    # no prior for the intercept 
                     X_prior_mean <- prior_mean      # sigma < 0 means noninformative prior
                     Y <- Y
                     Z <- X                          # Is this ok? 
                     J <- length(levels(inputdata[[group_column]]))
                     group <- as.integer(inputdata[[group_column]])
                   })
  
  params <- c(params, reg_params)
  
  return(params)
}

##################################################

mebn.localsummary <- function(fit)
{
  draws <- extract(fit)
  
  fixef_CI <- apply(draws$beta, 2, quantile, probs = c(0.05, 0.95)) # 2 = by column
  sigma_e <- mean(draws$sigma_e)

  ModelSummary <- within(list(),
                         {
                           intmean     <- mean(draws$beta_Intercept)
                           fixef       <- colMeans(draws$beta)   
                           fixef_l95CI <- fixef_CI[1,]
                           fixef_u95CI <- fixef_CI[2,]
                           ranef_sd    <- colMeans(draws$sigma_b)
                           std_error   <- sigma_e
                           Sigma_b     <- colMeans(draws$Sigma_b)
                           C           <- colMeans(draws$C)
                         })
  
  # Create matrix D
  sdM <- diag(ModelSummary$ranef_sd)
  ModelSummary$D <- sdM %*% ModelSummary$C %*% t(sdM)
  
  return(ModelSummary)
}

##################################################

mebn.predict_nodevalue <- function(localsummary, random_effects, parent_values)
{
  # TODO: This should return a distribution

  X0 <- as.matrix(parent_values)  
  Z0 <- X0  
  
  # fixed effects
  beta.hat <- as.matrix(localsummary$fixef)
  int.hat <- as.matrix(localsummary$intmean) 
  
  # random effects
  b_int.hat <- random_effects[1]
  b.hat <- random_effects[2:length(random_effects)]
  
  # Prediction
  y.hat <- int.hat + b_int.hat + X0 %*% beta.hat + Z0 %*% b.hat
  
  return(y.hat)
}

##################################################

mebn.get_prediction_CI95 <- function(parent_values, error_variances)
{
  X0 <- as.matrix(parent_values)
  X0 <- cbind(1, X0)
  
  # assume that all predictors may have random effects
  Z0 <- X0

  C1 <- error_variances$C1
  C2 <- error_variances$C2
  C12 <- error_variances$C12
  
  err_sigma2 <- X0 %*% C1 %*% t(X0) + Z0 %*% C2 %*% t(Z0) + 2*X0 %*% C12 %*% t(Z0)
  err_sigma <- sqrt(err_sigma2)
  
  n <- nrow(X0)
  
  # 95% confidence interval from Normal distribution
  CI95 <- qnorm(0.975)*(err_sigma/sqrt(n))
  
  return(CI95)
}

##################################################

mebn.get_personal_arcweight <- function(target, graph, parent_values)
{
  # fixed effects
  beta.hat <- as.matrix(localsummary$fixef)

  # random effects
  b.hat <- ranefs[2:length(ranefs)]
  
  weight <- beta.hat + b.hat
  
  return(weight)
}

##################################################

mebn.predict_ranefs <- function(localsummary, new_data, new_response)
{
  normalize_values <- TRUE

  n<-nrow(new_data)
  k<-length(localsummary$ranef_sd)  # nb of ranefs
  
  if (normalize_values == TRUE)
  {
    # TODO: myös tänne condscale!
    
    X <- as.matrix(cbind(rep(1,n), apply(new_data, 2, scale, center = TRUE, scale = TRUE)))
    Y <- as.vector(scale(new_response, center = TRUE, scale = TRUE)[,1])
    
    # dirty fix - works only for sysdimet data
    X[,2] <- as.vector(new_data[,1])
    X[,3] <- as.vector(new_data[,2])
  }
  else
  {
    X<-as.matrix(cbind(1,new_data))
    Y<-as.vector(new_response)
  }
  
  D<-as.matrix(localsummary$D)
  Z<-X

  sigma<-localsummary$std_error
  
  # TODO: Update matching autocorrelation structure
  R<-diag(sigma^2,nrow=n,ncol=n)
  
  beta<-as.vector(c(localsummary$intmean, localsummary$fixef))

  ranef <- D %*% t(Z) %*% solve(Z %*% D %*% t(Z) + R) %*% as.matrix(Y-X %*% beta)
  return(ranef[,1])
}

##################################################

mebn.ranef_designmatrix <- function(fixef_matrix, groups, observations)
{
  X <- fixef_matrix
  n <- nrow(X)  
  p <- ncol(X)
  g <- groups
  
  Z <- matrix(0, nrow = n, ncol = p*g)
  
  group_observations <- observations
  group_begins <- seq(1,n,group_observations)
  
  k <- 1
  for (r in group_begins)
  {
    Z[r:(r+group_observations-1),k:(k+p-1)] <- X[r:(r+group_observations-1),1:p]
    k <- k + p
  }
  
  return(Z)
}

##################################################

mebn.ranef_BLUP_vars <- function(localparams, localsummary)
{
  # Number of fixed and random effects 
  p <- localparams$p
  k <- localparams$k
  
  n <- localparams$N  # used to calculate R
  groups <- length(unique(localparams$group))
  group_obs <- n/groups
  
  # Here X and Z are design matrices for whole data, not just one group
  X <- localparams$X     # fixef design matrix 
  Z <- mebn.ranef_designmatrix(X,groups,group_obs)
  
  # We also need D for whole data
  D <- as.matrix(localsummary$D)
  D <- diag(1,groups) %x% D # block diagonal matrix for every group with Kronecker product

  sigma<-localsummary$std_error
  
  R <- diag(sigma^2,nrow=n,ncol=n) # residual var-cov matrix 

  # Henderson's Mixed-Effect Equations
  
  H11 <- t(X)%*%solve(R)%*%X
  H12 <- t(X)%*%solve(R)%*%Z
  H22 <- t(Z)%*%solve(R)%*%Z+solve(D)
  
  H <- cbind(rbind(H11,t(H12)),rbind(H12,H22))
  
  # Estimates of beta and b could be calculated from MEE 
  #y <- localparams$Y
  #r <- rbind(t(X)%*%solve(R)%*%y,t(Z)%*%solve(R)%*%y)
  #print(solve(H)%*%r)

  # The inverse of H provides variance-covariance matrix of the estimation errors
  Hinv <- solve(H)
  
  # Submatrices of Hinv
  C1 <- Hinv[1:p,1:p]                   # var(beta.hat) 
  C2 <- Hinv[(p+1):(p+k),(p+1):(p+k)]   # var(b.tilde-b), var(BLUP)
  C12 <- Hinv[(p+1):(p+k),1:p]          # cov(beta.hat, b.tilde-b) 
  
  # VarCov of errors for all the observations
  #err_varcov <- X %*% C1 %*% t(X) + Z %*% C2 %*% t(Z) + 2*X %*% C12 %*% t(Z)
  
  var_BLUP <- diag(C2)

  return(list(var_BLUP=var_BLUP, Hinv = Hinv, C1 = as.matrix(C1), C2 = as.matrix(C2), C12 = as.matrix(C12)))
}

##################################################

mebn.get_localfit <- function(target_name, local_model_cache = "models")
{
  modelcache <- paste0(local_model_cache, "\\", target_name, "_blmm", ".rds")
  print(paste0("Loading ", modelcache))
  
  localfit <- NULL
  
  if (file.exists(modelcache))
  {
    localfit <- readRDS(modelcache)
  }
  
  return(localfit)
}

##################################################

mebn.get_parents <- function(g, nodename)
{
  vindex <- as.numeric(V(g)[nodename])
  parents <- neighbors(g, vindex, mode = c("in"))
  
  return(parents)
}

##################################################

mebn.get_parents_with_type <- function(g, nodename, type)
{
  vindex <- as.numeric(V(g)[nodename])
  parents <- neighbors(g, vindex, mode = c("in"))
  
  # filter by type
  parents <- parents[parents$type==type]
  
  return(parents)
}

##################################################

mebn.sampling <- function(inputdata, predictor_columns, target_column, group_column, local_model_cache = "models", stan_mode_file = "BLMM.stan", normalize_values = TRUE, reg_params = NULL)
{
  require(rstan)
  
  # Run Stan parallel on multiple cores
  rstan_options(auto_write=TRUE)
  options(mc.cores=parallel::detectCores()) 
  
  target_name <- as.vector(target_column$Name)

  localfit <- mebn.get_localfit(target_name, local_model_cache)
  
  # Use cached model if it exists
  if (is.null(localfit))
  {
    stanDat <- mebn.set_model_parameters(predictor_columns, target_column, group_column, inputdata, normalize_values, reg_params)
    
    localfit <- stan(file=stan_mode_file, data=stanDat, warmup = 1000, iter=2500, chains=4, control = list(adapt_delta = 0.95, max_treedepth = 15))
    
    modelcache <- paste0(local_model_cache, "\\", target_name, "_blmm", ".rds")
    saveRDS(localfit, file=modelcache)
  }
  
  return(localfit)
}  

##################################################

mebn.variational_inference <- function(inputdata, predictor_columns, target_column, group_column, local_model_cache = "models", stan_mode_file = "BLMM.stan", normalize_values = TRUE, reg_params = NULL)
{
  require(rstan)
  
  # Run Stan parallel on multiple cores
  rstan_options (auto_write=TRUE)
  options (mc.cores=parallel::detectCores ()) 
  
  target_name <- as.vector(target_column$Name)
  
  localfit <- mebn.get_localfit(target_name, local_model_cache)
  
  # Use cached model if it exists
  if (is.null(localfit))
  {
    stanDat <- mebn.set_model_parameters(predictor_columns, target_column, group_column, inputdata, normalize_values, reg_params)

    print(stanDat)
        
    localmodel <- stan_model(file = stan_mode_file)
    localfit <- vb(localmodel, data=stanDat, output_samples=2500, iter=10000, seed=123)
    
    modelcache <- paste0(local_model_cache, "\\", target_name, "_blmm", ".rds")
    saveRDS(localfit, file=modelcache)
  }
  
  return(localfit)
}  



##################################################

mebn.iterate <- function(reaction_graph, inputdata, predictor_columns, assumed_targets, group_column, local_model_cache, stan_model_file, local_estimation, edge_significance_test, normalize_values = TRUE, use_regularization = FALSE)
{
  for (c in 1:dim(assumed_targets)[1])
  {
    target_column <- assumed_targets[c,]
    target_name <- as.vector(target_column$Name)
    
    localfit <- local_estimation(inputdata, predictor_columns, target_column, group_column, local_model_cache, stan_model_file, normalize_values, use_regularization)
    
    # Extract model summary
    localsummary <- mebn.localsummary(localfit)
    # localsummary
    
    if (FALSE) {
      
      # TODO: Missä tapauksessa nodella on nämä tiedot?
    
      vindex <- as.numeric(V(reaction_graph)[target_name])
      reaction_graph <- set_vertex_attr(reaction_graph, "mean", index = vindex, value = localsummary$target_mean)
      reaction_graph <- set_vertex_attr(reaction_graph, "l95CI", index = vindex, value = localsummary$target_l95CI)
      reaction_graph <- set_vertex_attr(reaction_graph, "u95CI", index = vindex, value = localsummary$target_u95CI)
      reaction_graph <- set_vertex_attr(reaction_graph, "scalemin", index = vindex, value = target_column$ScaleMin)
      reaction_graph <- set_vertex_attr(reaction_graph, "scalemax", index = vindex, value = target_column$ScaleMax)
  
    }    

    # - Loop through betas for current target
    predictor_names <- as.vector(predictor_columns$Name)
    
    for (p in 1:length(predictor_names))
    {
      predictor_name <- predictor_names[p]
      
      if (FALSE) {
        
        # TODO: Missä tapauksessa nodella on nämä tiedot?

        vindex <- as.numeric(V(reaction_graph)[predictor_name])     
        reaction_graph <- set_vertex_attr(reaction_graph, "mean", index = vindex, value = localsummary$parent_mean[p])
        reaction_graph <- set_vertex_attr(reaction_graph, "l95CI", index = vindex, value = localsummary$parent_l95CI[p])
        reaction_graph <- set_vertex_attr(reaction_graph, "u95CI", index = vindex, value = localsummary$parent_u95CI[p])
        reaction_graph <- set_vertex_attr(reaction_graph, "scalemin", index = vindex, value = predictor_columns[p,]$ScaleMin)
        reaction_graph <- set_vertex_attr(reaction_graph, "scalemax", index = vindex, value = predictor_columns[p,]$ScaleMax)
        
        #reaction_graph <- set_vertex_attr(reaction_graph, "scalemin", index = vindex, value = 1500.0)
        #reaction_graph <- set_vertex_attr(reaction_graph, "scalemax", index = vindex, value = 3000.0)
      }
      
      # Add significant edges between variables
      if (edge_significance_test(localsummary, p) == TRUE)
      {
        # Muutettu confband-kaari näyttämään betan vahvuutta ja CI esittää henk.koht. vaihtelua
        
        reaction_graph <- reaction_graph + edge(c(predictor_name, target_name), 
                                                weight = localsummary$fixef[p], 
                                                shape   = "confband", 
                                                mean   = localsummary$fixef[p],
                                                l95CI  = localsummary$fixef[p] - localsummary$ranef_sd[p+1]/2,
                                                u95CI  = localsummary$fixef[p] + localsummary$ranef_sd[p+1]/2)
        
        # reaction_graph <- reaction_graph + edge(c(predictor_name, target_name), 
        #                                         weight = localsummary$fixef[p], 
        #                                         shape   = "confband", 
        #                                         mean   = localsummary$fixef[p],
        #                                         l95CI  = localsummary$fixef_l95CI[p],
        #                                         u95CI  = localsummary$fixef_u95CI[p])
        
        # Add random-effect for significant predictors
        #reaction_graph <- reaction_graph + vertex(paste0("b_", predictor_name, "_", target_name), label=paste0("b_", predictor_name), color="#AAAAAA", size = 0.5, shape = "disc")
        #reaction_graph <- reaction_graph + vertex(paste0("b_sigma_", predictor_name, "_", target_name), label="b_sigma", color="#AAAAAA", size = localsummary$ranef_sd[p], shape = "disc")
        
        #reaction_graph <- reaction_graph + edge(paste0("b_sigma_", predictor_name, "_", target_name), paste0("b_", predictor_name, "_", target_name), shape = "arrow", weight = localsummary$ranef_sd[p]) 
        #reaction_graph <- reaction_graph + edge(paste0("b_", predictor_name, "_", target_name), target_name, shape = "arrow", weight=1)
      }
      else
      {
        #print(paste0("pruning out edge ", predictor_name, "->", target_name))
      }
    }
    
  } # loop targets
  
  return(reaction_graph)
}

##################################################

mebn.get_rootnodes <- function(g)
{
  which(sapply(sapply(V(g), function(x) neighbors(g,x, mode="in")), length) == 0)
}

##################################################

mebn.predict_graph <- function(reaction_graph, personal_data, alldata, predictor_columns, assumed_targets, local_model_cache = "models", edge_significance_test)
{
  predictor_names <- as.vector(predictor_columns$Name)
  new_input <- subset(personal_data, select = predictor_names)
  
  for (c in 1:dim(assumed_targets)[1])
  {
    target_column <- assumed_targets[c,]
    target_name <- as.vector(target_column$Name)
    
    localfit <- mebn.get_localfit(target_name, local_model_cache)
    
    if (is.null(localfit))
    {
      stop(paste0("Local model for ", target_name, " not found in '", local_model_cache, "'"))
    }
      
    localsummary <- mebn.localsummary(localfit)
    new_response <- subset(personal_data, select = target_name)
  
    # Predict random-effects with BLUP
    ranefs <- mebn.predict_ranefs(localsummary, new_input, new_response)
  
    # Get variances of BLUP, i.e. prediction errors
    localparams <- mebn.set_model_parameters(assumedpredictors, assumedtargets, "notinuse", alldata, TRUE, TRUE)
    err <- mebn.ranef_BLUP_vars(localparams, localsummary)
    ranef_BLUP_vars <- err$var_BLUP
    
    # - Loop through betas for current target
    for (p in 1:length(predictor_names))
    {
      predictor_name <- predictor_names[p]
      
      predicted_coef <- localsummary$fixef[p] + ranefs[p+1]
      prediction_error <- ranef_BLUP_vars[p]
      
      if (edge_significance_test(predicted_coef) == TRUE)
      {
        # Muutettu confband-kaari näyttämään betan vahvuutta ja CI esittää henk.koht. vaihtelua
        
        reaction_graph <- reaction_graph + edge(c(predictor_name, target_name), 
                                                weight = predicted_coef, 
                                                shape   = "confband", 
                                                mean   = predicted_coef,
                                                l95CI  = prediction_error/2,
                                                u95CI  = prediction_error/2)
        
        # Add random-effect for significant predictors
        #reaction_graph <- reaction_graph + vertex(paste0("b_", predictor_name, "_", target_name), label=paste0("b_", predictor_name), color="#AAAAAA", size = 0.5, shape = "disc")
        #reaction_graph <- reaction_graph + vertex(paste0("b_sigma_", predictor_name, "_", target_name), label="b_sigma", color="#AAAAAA", size = localsummary$ranef_sd[p], shape = "disc")
        
        #reaction_graph <- reaction_graph + edge(paste0("b_sigma_", predictor_name, "_", target_name), paste0("b_", predictor_name, "_", target_name), shape = "arrow", weight = localsummary$ranef_sd[p]) 
        #reaction_graph <- reaction_graph + edge(paste0("b_", predictor_name, "_", target_name), target_name, shape = "arrow", weight=1)
      }
      else
      {
        #print(paste0("pruning out edge ", predictor_name, "->", target_name))
      }
    }
    
  } # loop targets
  
  return(reaction_graph)
}

##################################################

mebn.write_gexf <- function(reaction_graph, gexf_path = "C:\\projects\\Responses\\reactions.gexf")
{
  require(rgexf)
  
  # TODO: https://github.com/jacomyal/sigma.js/wiki/Settings
  
  # - hover descriptions for nodes (current vs bounds) 
  # - beta-weights to edges and labels?
  
  MakeRGBA <- function(RGBstring, alpha)
  {
    strtodec <- function(rgb, b, e) { strtoi(paste0("0x", substr(rgb, b, e))) }
    
    RGBA <- data.frame(strtodec(RGBstring, 2, 3), strtodec(RGBstring, 4, 5), strtodec(RGBstring, 6, 7), alpha)  
    colnames(RGBA) <- c("r", "g", "b", "alpha")  
    
    return(RGBA)
  }
  
  graphdata <- get.data.frame(reaction_graph)
  
  nodeviz <- list(color = MakeRGBA(V(reaction_graph)$color, 1.0), size = V(reaction_graph)$size, shape = V(reaction_graph)$shape)
  edgeviz <- list(shape = E(reaction_graph)$shape)
  
  edgesatt <- data.frame(E(reaction_graph)$mean, E(reaction_graph)$l95CI, E(reaction_graph)$u95CI, E(reaction_graph)$b_sigma)
  
  if (length(edgesatt) == 4)
  {
    colnames(edgesatt) <- c("mean", "l95CI", "u95CI", "b_sigma")
  }
  
  nodesatt <- data.frame(V(reaction_graph)$mean, V(reaction_graph)$l95CI, V(reaction_graph)$u95CI, V(reaction_graph)$scalemin, V(reaction_graph)$scalemax)
  if (length(nodesatt) == 5)  
  {
    colnames(nodesatt) <- c("mean", "l95CI", "u95CI", "scalemin", "scalemax")
  }
  
  edgelabels <- data.frame(paste0("beta = ", round(E(reaction_graph)$mean, 3), ", b_sigma = ", round(2*(E(reaction_graph)$u95CI - E(reaction_graph)$mean), 3)))
  
  write.gexf(
    defaultedgetype = "directed",
    nodes = data.frame(V(reaction_graph)$name, V(reaction_graph)$label),
    edges = get.edgelist(reaction_graph),
    edgesWeight = graphdata[,3],
    edgesLabel = edgelabels,
    nodesVizAtt = nodeviz,
    edgesVizAtt = edgeviz,
    edgesAtt = edgesatt,
    nodesAtt = nodesatt,
    output = gexf_path
  )
}

###################################

mebn.BetaLevelTest <- function(LocalModelSummary, PredictorId)
{
  abs(LocalModelSummary$fixef[PredictorId]) > 0.001
}

###################################

mebn.RanefTest <- function(localsummary, PredictorId)
{
  abs(localsummary$fixef[PredictorId]) > 0.001 ||
    abs(localsummary$ranef_sd[PredictorId]) > 0.001
}

###################################

mebn.PersonalSignificanceTest <- function(personal_coef)
{
  abs(personal_coef) > 0.001
}

##################################################

mebn.typical_graph <- function(reaction_graph, inputdata, predictor_columns, assumed_targets, group_column, local_model_cache, stan_model_file, local_estimation, edge_significance_test, normalize_values = TRUE, reg_params = NULL)
{
  for (c in 1:dim(assumed_targets)[1])
  {
    target_column <- assumed_targets[c,]
    target_name <- as.vector(target_column$Name)
    
    localfit <- local_estimation(inputdata, predictor_columns, target_column, group_column, local_model_cache, stan_model_file, normalize_values, reg_params)
    
    # Extract model summary
    localsummary <- mebn.localsummary(localfit)

    # - Loop through betas for current target
    predictor_names <- as.vector(predictor_columns$Name)
    
    for (p in 1:length(predictor_names))
    {
      predictor_name <- predictor_names[p]
      
      # Add significant edges between variables
      if (edge_significance_test(localsummary, p) == TRUE)
      {
        # All the knowledge is stored to the graph
        
        # Attach the random variable
        reaction_graph <- reaction_graph + edge(c(predictor_name, target_name), 
                                                weight = localsummary$fixef[p], 
                                                b_sigma = localsummary$ranef_sd[p+1],
                                                shape   = "confband", 
                                                mean   = localsummary$fixef[p],
                                                l95CI  = localsummary$fixef[p] - localsummary$ranef_sd[p+1]/2,
                                                u95CI  = localsummary$fixef[p] + localsummary$ranef_sd[p+1]/2)
        
        # Fixed-effect
        reaction_graph <- reaction_graph + vertex(paste0("beta_", predictor_name, "_", target_name), 
                                                  label=paste0("beta_", predictor_name), 
                                                  type="beta", color="#AAAAAA", 
                                                  value = localsummary$fixef[p], 
                                                  value_l95CI = localsummary$fixef_l95CI[p],
                                                  value_u95CI = localsummary$fixef_u95CI[p],
                                                  shape = "circle")
        
        reaction_graph <- reaction_graph + edge(paste0("beta_", predictor_name, "_", target_name), paste0("beta_", predictor_name, "_", target_name), shape = "arrow", weight = 1, type = "beta") 
        
        # Add random-effect for significant predictors
        reaction_graph <- reaction_graph + vertex(paste0("b_", predictor_name, "_", target_name), 
                                                  label=paste0("b_", predictor_name), 
                                                  type="b", 
                                                  color="#AAAAAA", 
                                                  size = 0.5, 
                                                  shape = "circle")
        
        reaction_graph <- reaction_graph + vertex(paste0("b_sigma_", predictor_name, "_", target_name), 
                                                  label="b_sigma", 
                                                  type="b_sigma", 
                                                  color="#AAAAAA", 
                                                  value = localsummary$ranef_sd[p],
                                                  size = localsummary$ranef_sd[p], 
                                                  shape = "circle")
        
        reaction_graph <- reaction_graph + edge(paste0("b_sigma_", predictor_name, "_", target_name), paste0("b_", predictor_name, "_", target_name), shape = "arrow", weight = localsummary$ranef_sd[p], type = "b_sigma") 
        reaction_graph <- reaction_graph + edge(paste0("b_", predictor_name, "_", target_name), target_name, shape = "arrow", weight=1, type = "b")
      }
      else
      {
        #print(paste0("pruning out edge ", predictor_name, "->", target_name))
      }
    }
    
  } # loop targets
  
  return(reaction_graph)
}

##################################################

mebn.specific_graph <- function(reaction_graph, personal_data, alldata, predictor_columns, assumed_targets, local_model_cache = "models", edge_significance_test)
{
  predictor_names <- as.vector(predictor_columns$Name)
  new_input <- subset(personal_data, select = predictor_names)
  
  for (c in 1:dim(assumed_targets)[1])
  {
    target_column <- assumed_targets[c,]
    target_name <- as.vector(target_column$Name)
    
    localfit <- mebn.get_localfit(target_name, local_model_cache)
    
    if (is.null(localfit))
    {
      stop(paste0("Local model for ", target_name, " not found in '", local_model_cache, "'"))
    }
    
    localsummary <- mebn.localsummary(localfit)
    new_response <- subset(personal_data, select = target_name)
    
    # Predict random-effects with BLUP
    ranefs <- mebn.predict_ranefs(localsummary, new_input, new_response)
    
    # Get variances of BLUP, i.e. prediction errors
    localparams <- mebn.set_model_parameters(assumedpredictors, assumedtargets, "notinuse", alldata, TRUE, TRUE)
    err <- mebn.ranef_BLUP_vars(localparams, localsummary)
    ranef_BLUP_vars <- err$var_BLUP
    
    if (!is.null(prediction_data))
    {
      # Set predicted value and prediction error for target node
      prediction_input <- subset(prediction_data, select = predictor_names)
      targetvalue <- mebn.predict_nodevalue(localsummary, ranefs, prediction_input)
      target_95CI <- mebn.get_prediction_CI95(personal_data, err)
      
      vindex <- as.numeric(V(reaction_graph)[target_name])
      reaction_graph <- set_vertex_attr(reaction_graph, "mean", index = vindex, value = targetvalue)
      reaction_graph <- set_vertex_attr(reaction_graph, "l95CI", index = vindex, value = targetvalue-target_95CI)
      reaction_graph <- set_vertex_attr(reaction_graph, "u95CI", index = vindex, value = targetvalue+target_95CI)
      reaction_graph <- set_vertex_attr(reaction_graph, "scalemin", index = vindex, value = target_column$ScaleMin)
      reaction_graph <- set_vertex_attr(reaction_graph, "scalemax", index = vindex, value = target_column$ScaleMax)
    }
    
    # - Loop through betas for current target
    for (p in 1:length(predictor_names))
    {
      predictor_name <- predictor_names[p]
      
      predicted_coef <- localsummary$fixef[p] + ranefs[p+1]
      prediction_error <- ranef_BLUP_vars[p]
      
      if (edge_significance_test(predicted_coef) == TRUE)
      {
        # Muutettu confband-kaari näyttämään betan vahvuutta ja CI esittää henk.koht. vaihtelua
        
        reaction_graph <- reaction_graph + edge(c(predictor_name, target_name), 
                                                weight = predicted_coef, 
                                                shape   = "confband", 
                                                mean   = predicted_coef,
                                                l95CI  = prediction_error/2,
                                                u95CI  = prediction_error/2)
        
        # TODO: Should prior distributions of predictor be placed here? 
        if (FALSE) {
          vindex <- as.numeric(V(reaction_graph)[predictor_name])     
          reaction_graph <- set_vertex_attr(reaction_graph, "mean", index = vindex, value = localsummary$parent_mean[p])
          reaction_graph <- set_vertex_attr(reaction_graph, "l95CI", index = vindex, value = localsummary$parent_l95CI[p])
          reaction_graph <- set_vertex_attr(reaction_graph, "u95CI", index = vindex, value = localsummary$parent_u95CI[p])
          reaction_graph <- set_vertex_attr(reaction_graph, "scalemin", index = vindex, value = predictor_columns[p,]$ScaleMin)
          reaction_graph <- set_vertex_attr(reaction_graph, "scalemax", index = vindex, value = predictor_columns[p,]$ScaleMax)
          
          #reaction_graph <- set_vertex_attr(reaction_graph, "scalemin", index = vindex, value = 1500.0)
          #reaction_graph <- set_vertex_attr(reaction_graph, "scalemax", index = vindex, value = 3000.0)
        }
        
        # Add random-effect for significant predictors
        #reaction_graph <- reaction_graph + vertex(paste0("b_", predictor_name, "_", target_name), label=paste0("b_", predictor_name), color="#AAAAAA", size = 0.5, shape = "disc")
        #reaction_graph <- reaction_graph + vertex(paste0("b_sigma_", predictor_name, "_", target_name), label="b_sigma", color="#AAAAAA", size = localsummary$ranef_sd[p], shape = "disc")
        
        #reaction_graph <- reaction_graph + edge(paste0("b_sigma_", predictor_name, "_", target_name), paste0("b_", predictor_name, "_", target_name), shape = "arrow", weight = localsummary$ranef_sd[p]) 
        #reaction_graph <- reaction_graph + edge(paste0("b_", predictor_name, "_", target_name), target_name, shape = "arrow", weight=1)
      }
      else
      {
        print(paste0("pruning out edge ", predictor_name, "->", target_name))
      }
    }
    
  } # loop targets
  
  return(reaction_graph)
}

##################################################

mebn.visualization_graph <- function(mebn_graph)
{
  # All this information should be in attributes for visualization
  visual_graph <- mebn_graph
  
  # Remove edges to latent variables
  visual_graph <- delete.edges(visual_graph, which(E(visual_graph)$type=="beta"))
  visual_graph <- delete.edges(visual_graph, which(E(visual_graph)$type=="b_sigma"))
  visual_graph <- delete.edges(visual_graph, which(E(visual_graph)$type=="b"))
  
  # Remove nodes of latent variable
  visual_graph <- delete.vertices(visual_graph, which(V(visual_graph)$type=="beta"))
  visual_graph <- delete.vertices(visual_graph, which(V(visual_graph)$type=="b_sigma"))
  visual_graph <- delete.vertices(visual_graph, which(V(visual_graph)$type=="b"))
  
  return(visual_graph)
}

##################################################

mebn.set_evidence <- function(reaction_graph, evidence)
{
  # Evidence is a data frame with observed values of the predictors
  
}

##################################################

# Loop through the graph and evaluate

mebn.evaluate <- function(reaction_graph, targets)
{
  for (c in 1:dim(targets)[1])
  {
    target_column <- targets[c,]
    target_name <- as.vector(target_column$Name)
    
    localfit <- mebn.get_localfit(target_name, local_model_cache)
    
    if (is.null(localfit))
    {
      stop(paste0("Local model for ", target_name, " not found in '", local_model_cache, "'"))
    }
    
    localsummary <- mebn.localsummary(localfit)
    new_response <- subset(personal_data, select = target_name)
    
    # Predict random-effects with BLUP
    ranefs <- mebn.predict_ranefs(localsummary, new_input, new_response)
    
    # Get variances of BLUP, i.e. prediction errors
    localparams <- mebn.set_model_parameters(assumedpredictors, assumedtargets, "notinuse", alldata, TRUE, TRUE)
    err <- mebn.ranef_BLUP_vars(localparams, localsummary)
    ranef_BLUP_vars <- err$var_BLUP
    
    # TODO: Get evidence value
    
    if (!is.null(prediction_data))
    {
      # Set predicted value and prediction error for target node
      prediction_input <- subset(prediction_data, select = predictor_names)
      targetvalue <- mebn.predict_nodevalue(localsummary, ranefs, prediction_input)
      target_95CI <- mebn.get_prediction_CI95(personal_data, err)
      
      vindex <- as.numeric(V(reaction_graph)[target_name])
      reaction_graph <- set_vertex_attr(reaction_graph, "mean", index = vindex, value = targetvalue)
      reaction_graph <- set_vertex_attr(reaction_graph, "l95CI", index = vindex, value = targetvalue-target_95CI)
      reaction_graph <- set_vertex_attr(reaction_graph, "u95CI", index = vindex, value = targetvalue+target_95CI)
      reaction_graph <- set_vertex_attr(reaction_graph, "scalemin", index = vindex, value = target_column$ScaleMin)
      reaction_graph <- set_vertex_attr(reaction_graph, "scalemax", index = vindex, value = target_column$ScaleMax)
    }
    
    # - Loop through betas for current target
    for (p in 1:length(predictor_names))
    {
      predictor_name <- predictor_names[p]
      
      predicted_coef <- localsummary$fixef[p] + ranefs[p+1]
      prediction_error <- ranef_BLUP_vars[p]
      
      if (edge_significance_test(predicted_coef) == TRUE)
      {
        # Muutettu confband-kaari näyttämään betan vahvuutta ja CI esittää henk.koht. vaihtelua
        
        reaction_graph <- reaction_graph + edge(c(predictor_name, target_name), 
                                                weight = predicted_coef, 
                                                shape   = "confband", 
                                                mean   = predicted_coef,
                                                l95CI  = prediction_error/2,
                                                u95CI  = prediction_error/2)
        
        # TODO: Should prior distributions of predictor be placed here? 
        if (FALSE) {
          vindex <- as.numeric(V(reaction_graph)[predictor_name])     
          reaction_graph <- set_vertex_attr(reaction_graph, "mean", index = vindex, value = localsummary$parent_mean[p])
          reaction_graph <- set_vertex_attr(reaction_graph, "l95CI", index = vindex, value = localsummary$parent_l95CI[p])
          reaction_graph <- set_vertex_attr(reaction_graph, "u95CI", index = vindex, value = localsummary$parent_u95CI[p])
          reaction_graph <- set_vertex_attr(reaction_graph, "scalemin", index = vindex, value = predictor_columns[p,]$ScaleMin)
          reaction_graph <- set_vertex_attr(reaction_graph, "scalemax", index = vindex, value = predictor_columns[p,]$ScaleMax)
          
          #reaction_graph <- set_vertex_attr(reaction_graph, "scalemin", index = vindex, value = 1500.0)
          #reaction_graph <- set_vertex_attr(reaction_graph, "scalemax", index = vindex, value = 3000.0)
        }
        
        # Add random-effect for significant predictors
        #reaction_graph <- reaction_graph + vertex(paste0("b_", predictor_name, "_", target_name), label=paste0("b_", predictor_name), color="#AAAAAA", size = 0.5, shape = "disc")
        #reaction_graph <- reaction_graph + vertex(paste0("b_sigma_", predictor_name, "_", target_name), label="b_sigma", color="#AAAAAA", size = localsummary$ranef_sd[p], shape = "disc")
        
        #reaction_graph <- reaction_graph + edge(paste0("b_sigma_", predictor_name, "_", target_name), paste0("b_", predictor_name, "_", target_name), shape = "arrow", weight = localsummary$ranef_sd[p]) 
        #reaction_graph <- reaction_graph + edge(paste0("b_", predictor_name, "_", target_name), target_name, shape = "arrow", weight=1)
      }
      else
      {
        #print(paste0("pruning out edge ", predictor_name, "->", target_name))
      }
    }
    
  } # loop targets
  
  
  return(reaction_graph)
}

##################################################
