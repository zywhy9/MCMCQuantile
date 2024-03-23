source("../function/function.R")
library(ggplot2)
library(mcmcse)


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
  start.time <- Sys.time()
  Sigma11 <- varx * (1+phi)/(1-phi)
  Sigma22 <- varx2 * (1+phi^2)/(1-phi^2)
  Final_var <- Sigma11 + qnorm(0.975)^2/(4 * varx) * Sigma22
  
  sd <- sqrt(0.5^2 / (1 - 0.9^2) * (1 - phi^2))
  data <- ts_sampler(niter=nsim+nburnin, nchain=1, mu=0, sd=sd, phi=phi, seed=1234)
  data <- rbind(data,data^2)
  
  # out <- t(data[,(nburnin+1):(nsim+nburnin)])
  # est <- mcse.multi(x = out, method = "obm")
  # h%*%est$cov%*%h
  
  
  # res <- rep(0, length(102:(nsim+nburnin)))
  # for(i in 102:(nsim+nburnin)){
  #   out <- t(data[,100:i])
  #   est <- mcse.multi(x = out, method = "obm")
  #   res[i] <- h%*%est$cov%*%h
  # }
  
  res <- rep(0, length(100:nburnin))
  for(i in 100:nburnin){
    out <- t(data[,(i+1):(2*i)])
    est <- mcse.multi(x = out, method = "obm")
    res[i-99] <- h%*%est$cov%*%h
  }
  end.time <- Sys.time()
  plot(seq(200,20000,2),res, type="l", ylab="Estimated Variance", xlab="Iteration", main=paste0("Estimated Variance by Batch Means with phi=",phi, "and half burn-in"))
  abline(h=Final_var, col="red")
}