library(zeallot)
library(progress)
library(coda)

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
  if(length(dim(data))==3){
    if(is.null(niter)){
      c(niter, nchain, npar) %<-% dim(data)
    }else{
      c(nchain, npar) %<-% dim(data)[2:3]
    }
  }else{
    if(is.null(niter)){
      c(niter, nchain) %<-% dim(data)
      npar <- 1
    }else{
      nchain <- ncol(data)
      npar <- 1
    }
  }
  
  
  mcmc.samples <- matrix(0, nrow = niter, ncol = (npar * nchain)) ## Initialization
  if(npar>1){
    for(i in 1:npar){
      mcmc.samples[,(1:nchain)+((i-1)*nchain)] <- data[1:niter,,i]  ## Read Data
    }
  }else{
    mcmc.samples[,(1:nchain)] <- data[1:niter,]
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
#' @param short An indicator specifying if need shorter sequence under log10 scale.
#' @return A list including mean, SD, lower and upper bound of 95% CI, computing time.
mm_analysis <- function(data, niter=NULL, npar=NULL, nchain=NULL, transf=NULL, short=FALSE){
  
  if(short){
    actual_iter <- unique(round(10^seq(0,log10(niter),length.out=100)))
  }else{
    actual_iter <- 1:niter
  }
  actual_niter <- length(actual_iter)
  
  pb <- progress_bar$new(format = "[:bar] :current/:total (:percent)", total = actual_niter) ## Set up prograss bar
  
  mean_mm <- sd_mm <- low_mm <- upp_mm <- matrix(0, actual_niter, npar) ## Initialization
  
  starttime_mm <- Sys.time()  ## Set up start time
  
  ## Transformation based on type of data
  for(i in 1:npar){
    if(transf[i]=="log"){
      data[actual_iter,(1:nchain)+((i-1)*nchain)] <- 
        log(data[actual_iter,(1:nchain)+((i-1)*nchain)]) ## Log transformation
    }else if(transf[i]=="logit"){
      data[actual_iter,(1:nchain)+((i-1)*nchain)] <- 
        boot::logit(data[actual_iter,(1:nchain)+((i-1)*nchain)])  ## Logit transformation
    }else{
      next ## No transformation
    }
  }
  
  for(i in 1:actual_niter){
    iter <- actual_iter[i]
    tempres <- as.matrix(data[1:iter,])  ## Current data which only include first i iterations
    meanvec <- sdvec <- lowvec <- uppvec <- rep(0, npar) ## Initialization
    
    if(iter==1){ ## If only one iteration, then tempres is a vector
      for(j in 1:npar){
        meanvec[j] <- mean(tempres[(1:nchain)+((j-1)*nchain)])
        sdvec[j] <- sd(tempres[(1:nchain)+((j-1)*nchain)])
        lowvec[j] <- meanvec[j] - qnorm(0.975) * sdvec[j]
        uppvec[j] <- meanvec[j] + qnorm(0.975) * sdvec[j]
      }
    }else{  ## If more than one iteration, then tempres is a matrix
      for(j in 1:npar){
        meanvec[j] <- mean(tempres[,(1:nchain)+((j-1)*nchain)])
        sdvec[j] <- sd(tempres[,(1:nchain)+((j-1)*nchain)])
        lowvec[j] <- meanvec[j] - qnorm(0.975) * sdvec[j]
        uppvec[j] <- meanvec[j] + qnorm(0.975) * sdvec[j]
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
      }else{
        mean_mm[i,j] <- meanvec[j]
        sd_mm[i,j] <- sdvec[j]
        low_mm[i,j] <- lowvec[j]
        upp_mm[i,j] <- uppvec[j]
      }
    }
    
    pb$tick(1) ## Progress bar
  }
  endtime_mm <- Sys.time()  ## Record end time
  time_mm <- endtime_mm - starttime_mm ## Calculate computing time
  
  res <- list(mean_mm, sd_mm, low_mm, upp_mm, time_mm)
  if(short){
    res$actual_iter <- actual_iter
  }
  
  return(res)
}

#' Analysis with Empirical Quantiles
#'
#' @param data A matrix of MCMC samples.
#' @param niter Number of iterations.
#' @param npar Number of parameters.
#' @param nchain Number of MCMC chains.
#' @param short An indicator specifying if need shorter sequence under log10 scale.
#' @return A list including 2.5% quantile and 97.5% quantile, computing time.
eq_analysis <- function(data, niter=NULL, npar=NULL, nchain=NULL, short=FALSE){
  
  if(short){
    actual_iter <- unique(round(10^seq(0,log10(niter),length.out=100)))
  }else{
    actual_iter <- 1:niter
  }
  actual_niter <- length(actual_iter)
  
  pb <- progress_bar$new(format = "[:bar] :current/:total (:percent)", total = actual_niter) ## Set up progress bar
  
  low_eq <- upp_eq <- matrix(0, actual_niter, npar) ## Initialization
  
  starttime_eq <- Sys.time() ## Record start time
  
  for(i in 1:actual_niter){
    iter <- actual_iter[i]
    tempres <- as.matrix(data[1:iter,]) ## Current data which only include first i iterations
    lowvec <- uppvec <- rep(0, npar) ## Initialization
    
    if(iter==1){ ## If only one iteration, then tempres is a vector
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
  if(short){
    res$actual_iter <- actual_iter
  }
  
  return(res)
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
#' @param short An indicator specifying if need shorter sequence under log10 scale.
#' @return A list includes mean vector and vairance-covariance matrix.
mle_analysis <- function(data, niter=NULL, npar=NULL, nchain=NULL, transf=NULL, initial=NULL, extra=FALSE, nburnin=NULL, short=FALSE){
  
  ## Set initial values
  ninit <- ifelse(is.null(nburnin), npar, nburnin)
  
  if(short){
    actual_iter <- unique(round(10^seq(0,log10(niter),length.out=200)))
  }else{
    actual_iter <- 1:niter
  }
  actual_iter <- actual_iter[which(actual_iter>ninit)]
  actual_niter <- length(actual_iter)
  
  mean_mle <- sd_mle <- low_mle <- upp_mle <- matrix(0, actual_niter, npar) ## Initialization
  
  
  pb <- progress_bar$new(format = "[:bar] :current/:total (:percent)", total = actual_niter) ## Set progress bar
  starttime_mle <- Sys.time()  ## Set start time
  
  ## Transformation based on type of data
  for(i in 1:npar){
    if(transf[i]=="log"){
      data[actual_iter,(1:nchain)+((i-1)*nchain)] <- log(data[actual_iter,(1:nchain)+((i-1)*nchain)]) ## Log transformation
    }else if(transf[i]=="logit"){
      data[actual_iter,(1:nchain)+((i-1)*nchain)] <- boot::logit(data[actual_iter,(1:nchain)+((i-1)*nchain)])  ## Logit transformation
    }else{
      next ## No transformation
    }
  }
  
  if(is.null(initial)){
    temp <- rbind(rep(0, npar * nchain), data) ## If no initial values set, use 0
  }else{
    temp <- rbind(initial, data)  ## If initial values exist, use them
  }
  
  
  ## Start iteration
  for(i in 1:actual_niter){
    iter <- actual_iter[i]
    
    temp00 <- c() ## Initialization
    temp01 <- c() ## Initialization
    for(j in 1:npar){
      temp00 <- cbind(temp00, as.vector(temp[1:(iter-1),(1:nchain)+((j-1)*nchain)]))
      temp01 <- cbind(temp01, as.vector(temp[2:iter,(1:nchain)+((j-1)*nchain)]))
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
      for(j in 1:npar){
        mean_mle[i,j] <- NA
        sd_mle[i,j] <- NA
        low_mle[i,j] <- NA
        upp_mle[i,j] <- NA
      }
      pb$tick(1) ## Progress bar
      next
    }
    
    sigma_mle <- sqrt(diag(psi_mle))
    
    ## Transform back to original scale
    for(j in 1:npar){
      if(transf[j]=="log"){
        mean_mle[i,j] <- exp(mu_mle[j])
        sd_mle[i,j] <- sigma_mle[j] * mean_mle[i,j]
        low_mle[i,j] <- exp(mu_mle[j] - qnorm(0.975) * sigma_mle[j])
        upp_mle[i,j] <- exp(mu_mle[j] + qnorm(0.975) * sigma_mle[j])
      }else if(transf[j]=="logit"){
        mean_mle[i,j] <- boot::inv.logit(mu_mle[j])
        sd_mle[i,j] <- sigma_mle[j] * mean_mle[i,j] * (1 - mean_mle[i,j])
        low_mle[i,j] <- boot::inv.logit(mu_mle[j] - qnorm(0.975) * sigma_mle[j])
        upp_mle[i,j] <- boot::inv.logit(mu_mle[j] + qnorm(0.975) * sigma_mle[j])
      }else{
        mean_mle[i,j] <- mu_mle[j]
        sd_mle[i,j] <- sigma_mle[j]
        low_mle[i,j] <- mu_mle[j] - qnorm(0.975) * sigma_mle[j]
        upp_mle[i,j] <- mu_mle[j] + qnorm(0.975) * sigma_mle[j]
      }
    }
    pb$tick(1) ## Progress bar
  }
  
  if(!short){
    mean_mle <- rbind(matrix(rep(0, ninit),ncol=1), mean_mle)
    sd_mle <- rbind(matrix(rep(0, ninit),ncol=1), sd_mle)
    low_mle <- rbind(matrix(rep(0, ninit),ncol=1), low_mle)
    upp_mle <- rbind(matrix(rep(0, ninit),ncol=1), upp_mle)
  }
  
  endtime_mle <- Sys.time() ## Record end time
  time_mle <- endtime_mle - starttime_mle ## Calculate computing time
  
  res <- list(mean_mle, sd_mle, low_mle, upp_mle, time_mle)
  
  if(extra){
    res$mu_mle <- mu_mle
    res$psi_mle <- psi_mle
    res$beta <- beta
    res$gamma <- gamma 
    res$delta <- delta
  }
  if(short){
    res$actual_iter <- actual_iter
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


#' Analysis with High Posterior Density Interval
#'
#' @param data A matrix of MCMC samples.
#' @param niter Number of iterations.
#' @param npar Number of parameters.
#' @param nchain Number of MCMC chains.
#' @param transf A vector indicates transformation if needed.
#' @return A list including mean, SD, lower and upper bound of 95% CI, computing time.
hpd_analysis <- function(data, niter=NULL, npar=NULL, nchain=NULL, transf=NULL){
  
  pb <- progress_bar$new(format = "[:bar] :current/:total (:percent)", total = niter - 1) ## Set up prograss bar
  
  low_hpd <- upp_hpd <- matrix(0, niter, npar) ## Initialization
  
  starttime_hpd <- Sys.time()  ## Set up start time
  
  ## Transformation based on type of data
  for(i in 1:npar){
    if(transf[i]=="log"){
      data[,(1:nchain)+((i-1)*nchain)] <- log(data[,(1:nchain)+((i-1)*nchain)]) ## Log transformation
    }else if(transf[i]=="logit"){
      data[,(1:nchain)+((i-1)*nchain)] <- boot::logit(data[,(1:nchain)+((i-1)*nchain)])  ## Logit transformation
    }else{
      next ## No transformation
    }
  }
  
  for(i in 2:niter){
    tempres <- as.matrix(data[1:i,])  ## Current data which only include first i iterations
    lowvec <- uppvec <- rep(0, npar) ## Initialization
    
    for(j in 1:npar){
      mcmc.obj <- as.mcmc(as.vector(tempres[,(1:nchain)+((j-1)*nchain)]))
      hpd <- HPDinterval(mcmc.obj)
      lowvec[j] <- hpd[1,1]
      uppvec[j] <- hpd[1,2]
    }
    
    ## Transform back to original scale
    for(j in 1:npar){
      if(transf[j]=="log"){
        low_hpd[i,j] <- exp(lowvec[j])
        upp_hpd[i,j] <- exp(uppvec[j])
      }else if(transf[j]=="logit"){
        low_hpd[i,j] <- boot::inv.logit(lowvec[j])
        upp_hpd[i,j] <- boot::inv.logit(uppvec[j])
      }else{
        low_hpd[i,j] <- lowvec[j]
        upp_hpd[i,j] <- uppvec[j]
      }
    }
    
    pb$tick(1) ## Progress bar
  }
  endtime_hpd <- Sys.time()  ## Record end time
  time_hpd <- endtime_hpd - starttime_hpd ## Calculate computing time
  
  res <- list(low_hpd, upp_hpd, time_hpd)
}

#' AR(1) Simulation
#'
#' @param niter Number of iterations.
#' @param nchain Number of MCMC chains.
#' @param mu Mean of Xt.
#' @param sd Standard deviation of Xt.
#' @param phi AR(1) model coefficient.
#' @param initial Initial value of AR(1) process.
#' @param seed Seed for random sampler.
#' @return A matrix with simulated data.
ts_sampler <- function(niter, nchain, mu, sd, phi, initial=NULL, seed = 1234) { ## x_t = delta + phi * x_{t-1} + eps_t
  set.seed(seed)
  
  samples <- matrix(0, nrow = (niter+1), ncol = nchain)
  ## If no initial value set, the default initial value is 0
  if(!is.null(initial)){
    samples[1,] <- initial
  }
  
  ## Simulation
  delta <- mu * (1 - phi)
  for(i in 2:(niter+1)){
    for(j in 1:nchain){
      samples[i,j] <- delta + phi * samples[(i-1),j] + rnorm(1, 0, sd)
    }
  }
  samples <- samples[-1,]
  return(samples)
}
