
#' Log likelihood of univariate AR(1) model
#'
#' @param data A vector or a matrix of AR(1) samples.
#' @param param A vector of parameters including phi, delta, and sigma.
#' @return The log-likelihood of the AR(1) model.
loglikar <- function(data, param){
  
  ## Parameters of interest
  phi <- param[1]
  delta <- param[2]
  sigma <- exp(param[3])
  
  
  if(is.null(dim(data))){  ## One chain
    t <- length(data)
    temp <- 0
    for(i in 2:t){
      temp <- temp + (data[i] - delta - phi * data[i-1])^2
    }
    res <- -t * log(2 * pi)/2 - 1/2 * log(sigma^2 / (1 - phi^2)) - 
      (1 - phi^2)/(2 * sigma^2) * (data[1] - delta/(1 - phi))^2 - (t-1) * log(sigma^2)/2 - temp/(2 * sigma^2)
  }else{  ## Multiple chains
    c(t, nchain) %<-% dim(data)
    res <- 0
    
    for(j in 1:nchain){
      temp <- 0
      for(i in 2:t){
        temp <- temp + (data[i,j] - delta - phi * data[i-1,j])^2
      }
      res <- res - t * log(2 * pi)/2 - 1/2 * log(sigma^2 / (1 - phi^2)) - 
        (1 - phi^2)/(2 * sigma^2) * (data[1,j] - delta/(1 - phi))^2 - 
        (t-1)/2 * log(sigma^2) - temp/(2 * sigma^2)
    }
  }
  return(res)
}

#' Using optimization to find the MLEs of univariate AR(1) model
#'
#' @param data A vector or a matrix of AR(1) samples.
#' @param inits A vector initial values for phi, delta, and sigma.
#' @return A vector containing MLE of phi, delta and sigma.
mle_univ <- function(data, inits = c(0.5, 0.5, log(0.5))){
  res <- optim(par=inits, fn=loglikar,
               data=data,
               control = list(fnscale=-1))
  out <- c(res$par[1:2], exp(res$par[3]))
  return(out)
}


#' Analysis with Maximum Likelihood Estimation using Univariate Optimization
#'
#' @param data A matrix of MCMC samples.
#' @param niter Number of iterations.
#' @param npar Number of parameters.
#' @param nchain Number of MCMC chains.
#' @param transf A vector indicates transformation if needed.
#' @param initial Optional; Initial values, if NULL then all 0.
#' @param extra If TRUE, then also return phi, delta, sigma
#' @param nburnin Optional; number of burn-in.
#' @return A list includes mean vector and vairance-covariance matrix.
mle_analysis_uni <- function(data, niter=NULL, npar=NULL, nchain=NULL, transf=NULL, initial=NULL, extra=FALSE, nburnin=NULL){
  
  mean_mle <- sd_mle <- low_mle <- upp_mle <- matrix(0, niter, npar) ## Initialization
  
  ## Set initial values
  ninit <- ifelse(is.null(nburnin), npar, nburnin)
  
  pb <- progress_bar$new(format = "[:bar] :current/:total (:percent)", total = niter-ninit) ## Set progress bar
  starttime_mle <- Sys.time()  ## Set start time
  
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
  
  if(is.null(initial)){
    temp <- rbind(rep(0, npar * nchain), data) ## If no initial values set, use 0
  }else{
    temp <- rbind(initial, data)  ## If initial values exist, use them
  }
  
  
  ## Start iteration
  for(iter in (ninit+1):niter){
    
    c(phi, delta, sigma) %<-% mle_univ(data = temp)
    
    mu_mle <- delta / (1 - phi)
    sigma_mle <- sigma
    
    ## If variance is negative, print the corresponding iteration number
    if(sigma_mle<0){
      print(iter)
      print(sigma_mle)
      break
      return()
    }
    
    ## Transform back to original scale
    for(j in 1:npar){
      if(transf[j]=="log"){
        mean_mle[iter,j] <- exp(mu_mle[j])
        sd_mle[iter,j] <- sigma_mle[j] * mean_mle[iter,j]
        low_mle[iter,j] <- exp(mu_mle[j] - qnorm(0.975) * sigma_mle[j])
        upp_mle[iter,j] <- exp(mu_mle[j] + qnorm(0.975) * sigma_mle[j])
      }else if(transf[j]=="logit"){
        mean_mle[iter,j] <- boot::inv.logit(mu_mle[j])
        sd_mle[iter,j] <- sigma_mle[j] * mean_mle[iter,j] * (1 - mean_mle[iter,j])
        low_mle[iter,j] <- boot::inv.logit(mu_mle[j] - qnorm(0.975) * sigma_mle[j])
        upp_mle[iter,j] <- boot::inv.logit(mu_mle[j] + qnorm(0.975) * sigma_mle[j])
      }else{
        mean_mle[iter,j] <- mu_mle[j]
        sd_mle[iter,j] <- sigma_mle[j]
        low_mle[iter,j] <- mu_mle[j] - qnorm(0.975) * sigma_mle[j]
        upp_mle[iter,j] <- mu_mle[j] + qnorm(0.975) * sigma_mle[j]
      }
    }
    pb$tick(1) ## Progress bar
  }
  
  endtime_mle <- Sys.time() ## Record end time
  time_mle <- endtime_mle - starttime_mle ## Calculate computing time
  
  if(extra){
    res <- list(mean_mle, sd_mle, low_mle, upp_mle, time_mle, phi, delta, sigma)
  }else{
    res <- list(mean_mle, sd_mle, low_mle, upp_mle, time_mle)
  }
  
  return(res)
}
