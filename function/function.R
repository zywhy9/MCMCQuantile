library(zeallot)
library(progress)

#' Build Dataset
#'
#' @param data An MCMC array from jagsfit$BUGSoutput$sims.array.
#' @param niter Number of iterations of interest.
#' @param set The set number of the dataset. For now only support one dataset a time.
#' @param location The location of the dataset. Default dataset name is 'res'+index+'.rds'.
#' @param param.name Name(s) of the parameter(s) of interest. if it's not null, data will only contain these parameters.
#' @return A list includes MCMC sample matrix, number of iterations, number of chains, number of parameters.

mcmc_data <- function(data=NULL, niter=NULL, set=1, location="res", param.name=NULL){
  
  ## If no dataset assigned, then find the dataset from its location
  if(is.null(data)){
    if(is.null(param.name)){
      data <- readRDS(paste0(location, "/res", set,".rds"))
      param.name <- dimnames(data)[[3]]
    }else{
      data <- readRDS(paste0(location, "/res", set,".rds"))[,,param.name]
    }
  }
  
  ## Find the number of iterations, number of chains and number of parameters
  if(is.null(niter)){
    c(niter, nchain, npar) %<-% dim(data)
  }else{
    c(nchain, npar) %<-% dim(data)[2:3]
  }
  
  mcmc.samples <- matrix(0, nrow = niter, ncol = (npar * nchain)) ## Initialization
  for(i in 1:npar){
    mcmc.samples[,(1:nchain)+((i-1)*nchain)] <- data[1:niter,,i]  ## Read Data
  }
  
  if(is.null(niter)){
    res <- list(data = mcmc.samples, niter = niter, nchain = nchain, npar = npar, param.name = param.name)
  }else{
    res <- list(data = mcmc.samples, nchain = nchain, npar = npar, param.name = param.name)
  }
  
  return(res)
}


#' Analysis with Method of Momemts
#'
#' @param data A matrix of MCMC samples.
#' @param niter Number of iterations.
#' @param npar Number of parameters.
#' @param nchain Number of MCMC chains.
#' @param transf A vector indicates transformation if needed.
#' @return A list including mean, SD, lower and upper bound of 95% CI, computing time.
mm_analysis <- function(data, niter=NULL, npar=NULL, nchain=NULL, transf=NULL){
  
  pb <- progress_bar$new(format = "[:bar] :current/:total (:percent)", total = niter) ## Set up prograss bar
  
  mean_mm <- sd_mm <- low_mm <- upp_mm <- matrix(0, niter, npar) ## Initialization
  
  starttime_mm <- Sys.time()  ## Set up start time
  
  ## Transformation based on type of data
  for(i in 1:npar){
    if(is.null(transf[i])){
      next
    }else{
      if(transf[i]=="log"){
        data[,(1:nchain)+((i-1)*nchain)] <- log(data[,(1:nchain)+((i-1)*nchain)])
      }else if(transf[i]=="logit"){
        data[,(1:nchain)+((i-1)*nchain)] <- boot::logit(data[,(1:nchain)+((i-1)*nchain)])
      }
    }
  }
  
  for(i in 1:niter){
    tempres <- data[1:i,]  ## Current data which only include first i iterations
    meanvec <- sdvec <- lowvec <- uppvec <- rep(0, npar) ## Initialization
    
    if(i==1){ ## If only one iteration, then tempres is a vector
      for(j in 1:npar){
        meanvec[j] <- mean(tempres[(1:nchain)+((j-1)*nchain)])
        sdvec[j] <- sd(tempres[(1:nchain)+((j-1)*nchain)])
        lowvec[j] <- meanvec[j] - 1.96 * sdvec[j]
        uppvec[j] <- meanvec[j] + 1.96 * sdvec[j]
      }
    }else{  ## If more than one iteration, then tempres is a matrix
      for(j in 1:npar){
        meanvec[j] <- mean(tempres[,(1:nchain)+((j-1)*nchain)])
        sdvec[j] <- sd(tempres[,(1:nchain)+((j-1)*nchain)])
        lowvec[j] <- meanvec[j] - 1.96 * sdvec[j]
        uppvec[j] <- meanvec[j] + 1.96 * sdvec[j]
      }
    }
    
    ## Transform back to original scale
    for(j in 1:npar){
      if(transf[j]=="log"){
        mean_mm[i,j] <- exp(meanvec[j])
        sd_mm[i,j] <- sdvec[j] * mean_mm[i,j]
        low_mm[i,j] <- exp(lowvec[j])
        upp_mm[i,j] <- exp(uppvec[j])
      }else if(transf[j]=="logit"){
        mean_mm[i,j] <- boot::inv.logit(meanvec[j])
        sd_mm[i,j] <- sdvec[j] * mean_mm[i,j] * (1 - mean_mm[i,j])
        low_mm[i,j] <- boot::inv.logit(lowvec[j])
        upp_mm[i,j] <- boot::inv.logit(uppvec[j])
      }
    }
    
    pb$tick(1) ## Progress bar
  }
  endtime_mm <- Sys.time()  ## Record end time
  time_mm <- endtime_mm - starttime_mm ## Calculate computing time
  
  res <- list(mean_mm, sd_mm, low_mm, upp_mm, time_mm)
}

#' Analysis with Empirical Quantiles
#'
#' @param data A matrix of MCMC samples.
#' @param niter Number of iterations.
#' @param npar Number of parameters.
#' @param nchain Number of MCMC chains.
#' @return A list including 2.5% quantile and 97.5% quantile, computing time.
eq_analysis <- function(data, niter=NULL, npar=NULL, nchain=NULL){
  
  pb <- progress_bar$new(format = "[:bar] :current/:total (:percent)", total = niter) ## Set up progress bar
  
  low_eq <- upp_eq <- matrix(0, niter, npar) ## Initialization
  
  starttime_eq <- Sys.time() ## Record start time
  
  for(i in 1:niter){
    tempres <- data[1:i,] ## Current data which only include first i iterations
    lowvec <- uppvec <- rep(0, npar) ## Initialization
    
    if(i==1){ ## If only one iteration, then tempres is a vector
      for(j in 1:npar){
        lowvec[j] <- quantile(tempres[(1:nchain)+((j-1)*nchain)], 0.025)
        uppvec[j] <- quantile(tempres[(1:nchain)+((j-1)*nchain)], 0.975)
      }
    }else{  ## If more than one iteration, then tempres is a matrix
      for(j in 1:npar){
        lowvec[j] <- quantile(tempres[,(1:nchain)+((j-1)*nchain)], 0.025)
        uppvec[j] <- quantile(tempres[,(1:nchain)+((j-1)*nchain)], 0.975)
      }
    }
    
    low_eq[i,] <- lowvec ## Save values
    upp_eq[i,] <- uppvec ## Save values
    
    pb$tick(1) ## Progress bar
  }
  endtime_eq <- Sys.time() ## Record end time
  time_eq <- endtime_eq - starttime_eq ## Calculate computing time
  
  res <- list(low_eq, upp_eq, time_eq)
}


#' Analysis with Maximum Likelihood Estimation
#'
#' @param data A matrix of MCMC samples.
#' @param niter Number of iterations.
#' @param npar Number of parameters.
#' @param nchain Number of MCMC chains.
#' @param transf A vector indicates transformation if needed.
#' @param initial Optional; Initial values, if NULL then all 0.
#' @param extra If TRUE, then also return mu, psi, beta, gamma, delta.
#' @param nburnin Optional; number of burn-in.
#' @return A list includes mean vector and vairance-covariance matrix.
mle_analysis <- function(data, niter=NULL, npar=NULL, nchain=NULL, transf=NULL, initial=NULL, extra=FALSE, nburnin=NULL){
  
  mean_mle <- sd_mle <- low_mle <- upp_mle <- matrix(0, niter, npar) ## Initialization
  
  pb <- progress_bar$new(format = "[:bar] :current/:total (:percent)", total = niter-(npar-1)) ## Set progress bar
  starttime_mle <- Sys.time()  ## Set start time
  
  ## Transformation based on type of data
  for(i in 1:npar){
    if(is.null(transf[i])){
      next ## No transformation
    }else{
      if(transf[i]=="log"){
        data[,(1:nchain)+((i-1)*nchain)] <- log(data[,(1:nchain)+((i-1)*nchain)]) ## Log transformation
      }else if(transf[i]=="logit"){
        data[,(1:nchain)+((i-1)*nchain)] <- boot::logit(data[,(1:nchain)+((i-1)*nchain)])  ## Logit transformation
      }
    }
  }
  
  if(is.null(initial)){
    temp <- rbind(rep(0, npar * nchain), data) ## If no initial values set, use 0
  }else{
    temp <- rbind(initial, data)  ## If initial values exist, use them
  }
  
  ## Set initial values
  ninit <- ifelse(is.null(nburnin), npar, nburnin)
  
  ## Start iteration
  for(iter in (ninit+1):niter){
    
    temp00 <- c() ## Initialization
    temp01 <- c() ## Initialization
    for(i in 1:npar){
      temp00 <- cbind(temp00, as.vector(temp[1:(iter-1),(1:nchain)+((i-1)*nchain)]))
      temp01 <- cbind(temp01, as.vector(temp[2:iter,(1:nchain)+((i-1)*nchain)]))
    }
    X0bar <- colMeans(temp00)
    X1bar <- colMeans(temp01)
    
    temp1 <- temp2 <- temp3 <- matrix(0, npar, npar) ## Initialization
    for(t in 1:iter){
      for(j in 1:nchain){
        temp1 <- temp1 + (temp[(t+1),(1:npar)*nchain - (nchain-j)] - X1bar) %*% t(temp[t,(1:npar)*nchain - (nchain-j)] - X0bar)
        temp2 <- temp2 + (temp[t,(1:npar)*nchain - (nchain-j)] - X0bar) %*% t(temp[t,(1:npar)*nchain - (nchain-j)] - X0bar)
      }
    }
    
    beta <- temp1 %*% solve(temp2)
    gamma <- X1bar - beta %*% X0bar
    
    for(t in 1:iter){
      for(j in 1:nchain){
        temp3 <- temp3 + 
          (temp[(t+1),1:npar*nchain - (nchain-j)] - beta %*% temp[t,1:npar*nchain - (nchain-j)] - gamma) %*% t(temp[(t+1),1:npar*nchain - (nchain-j)] - beta %*% temp[t,1:npar*nchain - (nchain-j)] - gamma)
      }
    }
    delta <- temp3/(iter*nchain)
    
    mu_mle <- solve(a = (diag(npar)-beta), b = gamma)
    
    A <- solve(beta)
    B <- -t(beta)
    C <- solve(beta) %*% delta
    psi_mle <- maotai::sylvester(A, B, C)
    
    ## If variance is negative, print the corresponding iteration number
    if(sum(diag(psi_mle)<0)>0){
      print(iter)
      print(psi_mle)
      break
      return()
    }
    
    sigma_mle <- sqrt(diag(psi_mle))
    
    ## Transform back to original scale
    for(j in 1:npar){
      if(transf[j]=="log"){
        mean_mle[iter,j] <- exp(mu_mle[j])
        sd_mle[iter,j] <- sigma_mle[j] * mean_mle[iter,j]
        low_mle[iter,j] <- exp(mu_mle[j] - 1.96 * sigma_mle[j])
        upp_mle[iter,j] <- exp(mu_mle[j] + 1.96 * sigma_mle[j])
      }else if(transf[j]=="logit"){
        mean_mle[iter,j] <- boot::inv.logit(mu_mle[j])
        sd_mle[iter,j] <- sigma_mle[j] * mean_mle[iter,j] * (1 - mean_mle[iter,j])
        low_mle[iter,j] <- boot::inv.logit(mu_mle[j] - 1.96 * sigma_mle[j])
        upp_mle[iter,j] <- boot::inv.logit(mu_mle[j] + 1.96 * sigma_mle[j])
      }else{
        mean_mle[iter,j] <- mu_mle[j]
        sd_mle[iter,j] <- sigma_mle[j]
        low_mle[iter,j] <- mu_mle[j] - 1.96 * sigma_mle[j]
        upp_mle[iter,j] <- mu_mle[j] + 1.96 * sigma_mle[j]
      }
    }
    pb$tick(1) ## Progress bar
  }
  
  endtime_mle <- Sys.time() ## Record end time
  time_mle <- endtime_mle - starttime_mle ## Calculate computing time

  if(extra){
    res <- list(mean_mle, sd_mle, low_mle, upp_mle, time_mle, mu_mle, psi_mle, beta, gamma, delta)
  }else{
    res <- list(mean_mle, sd_mle, low_mle, upp_mle, time_mle)
  }
  
  return(res)
}


#' Markov-Normal Simulation
#'
#' @param niter Number of iterations.
#' @param npar Number of parameters.
#' @param nchain Number of MCMC chains.
#' @param transf A vector indicates transformation if needed.
#' @param init Initial values, if NULL then all 0.
#' @param beta The last Beta from the MLE analysis.
#' @param gamma The last Gamma from the MLE analysis.
#' @param delta The last Delta from the MLE analysis.
#' @param seed Seed for random sampler.
#' @return A matrix with simulated data.
ar_simulation <- function(niter, npar, nchain, transf, init, beta, gamma, delta, seed=1234) {
  
  set.seed(seed)
  
  pb <- progress_bar$new(format = "[:bar] :current/:total (:percent)", total = niter) ## Set progress bar
  samples <- matrix(0, nrow = (niter + 1), ncol = (npar * nchain)) ## Initialize a matrix
  
  samples[1,] <- init ## The first row is initial values
  
  for(i in 2:(niter + 1)){
    for(j in 1:nchain){
      parid <- c(1:npar) * nchain - (nchain - j)  ## The parameter location for each chain
      samples[i,parid] <- mvrnorm(1, as.vector(beta%*%samples[(i-1), parid] + gamma), delta)  ## Follow the transition distribution in the paper
    }
    pb$tick(1) ## Add 1 to the progress bar
  }
  samples <- samples[-1,] ## Remove the initial values
  
  ## Transform back to the original scale
  for(i in 1:npar){
    if(is.null(transf[i])){
      next
    }else{
      if(transf[i]=="log"){
        samples[,(1:nchain)+((i-1)*nchain)] <- exp(samples[,(1:nchain)+((i-1)*nchain)])
      }else if(transf[i]=="logit"){
        samples[,(1:nchain)+((i-1)*nchain)] <- boot::inv.logit(samples[,(1:nchain)+((i-1)*nchain)])
      }
    }
  }
  
  return(samples)
}
