source("../function/function.R")
library(ggplot2)
library(mcmcse)
library(timsac)


true_sd <- sqrt(0.5^2 / (1 - 0.9^2)) ## SD of stationary distribution
true_4rm <- 3 * true_sd^4 ## Fourth moment of stationary distribution
varx <- true_sd^2 ## Variance of stationary distribution 
varx2 <- true_4rm - varx

nsim <- 10000
nburnin <- 10000
h <- c(1, qnorm(0.975)/(2*true_sd))

par(mfrow=c(3,1))
phichoice <- c(0.4,0.6,0.9)

for(phi in phichoice){
  Sigma11 <- varx * (1+phi)/(1-phi)
  Sigma22 <- varx2 * (1+phi^2)/(1-phi^2)
  Final_var <- Sigma11 + qnorm(0.975)^2/(4 * varx) * Sigma22
  
  sd <- sqrt(0.5^2 / (1 - 0.9^2) * (1 - phi^2))
  data <- ts_sampler(niter=nsim, nchain=1, mu=0, sd=sd, phi=phi, seed=1234)
  data <- rbind(data,data^2)[,101:10000]
  
  res <- rep(0, 9900)
  for(i in 1:9900){
    temp <- autoarmafit(data[1,])
    rho <- ARMAacf(ar=temp$best.model$arcoef, ma=temp$best.model$macoef, lag.max = ncol(data))
    res[i] <- (sum(rho[1:(i+1)])*2-1)*varx
  }
  
  plot(101:10000, res, type="l", ylab="Estimated Variance", xlab="Iteration", main=paste0("Estimated Variance by Batch Means with phi=",phi," and half burn-in"))
  abline(h=Sigma11, col="red")
}
