source("../function/function.R")
library(ggplot2)
library(scales)
library(Matrix)
library(rstan)
library(ggpubr)


## autocorrelation
phi <- 0.9
sd <- 0.5 
true_sd <- sqrt(0.5^2 / (1 - 0.9^2)) ## SD of stationary distribution
true_4rm <- 3 * true_sd^4 ## Fourth moment of stationary distribution
varx <- true_sd^2 ## Variance of stationary distribution 
varx2 <- true_4rm - varx
Sigma11 <- varx * (1+phi)/(1-phi)
Sigma22 <- varx2 * (1+phi^2)/(1-phi^2)
Final_var <- Sigma11 + qnorm(0.975)^2/(4 * varx) * Sigma22


## Simulation
nsim <- 5000 ## Number of iterations
nset <- 100 ## Number of simulated sets
simest <- matrix(0, ncol = nset, nrow = nsim - 1)
for(set in 1:nset){
  temp <- readRDS(paste0("CLT/res",set,".rds")) ## Read Data
  simest[,set] <- as.vector(temp$upp_mm[2:nsim,1])
}
simvar <- apply(simest, 1, var)

# simvar1 <- rep(0, nsim - 1)
# simvar1[1] <- var(simest[1,])
# pb <- txtProgressBar(min = 2,      # Minimum value of the progress bar
#                      max = nsim-1, # Maximum value of the progress bar
#                      style = 3,    # Progress bar style (also available style = 1 and style = 2)
#                      width = 50,   # Progress bar width. Defaults to getOption("width")
#                      char = "=")   # Character used to create the bar
# for(i in 2:(nsim - 1)){
#   mean1 <- apply(simest[1:i,], 2, mean)
#   B <- i * var(mean1)
#   sm <- apply(simest[1:i,], 2, var)
#   W <- mean(sm)
#   simvar1[i] <- (i - 1) / i * W + B / i
#   setTxtProgressBar(pb, i)
# }
# close(pb) # Close the connection

## Plot
hex <- hue_pal()(2)

df <- data.frame(id = rep(2:nsim, 2), 
                 var = c(Final_var / c(2:nsim), simvar),
                 Type = as.factor(rep(1:2, each = (nsim - 1))))

plot <- ggplot(df, aes(x = id, group = Type)) +  # Plot for lower bound
  geom_line(aes(y = var, col = Type, linetype = Type), lwd=0.75) +
  ylim(c(0, 1)) +
  labs(x = "Iterations", y = "Variance of method of moments estimator") +
  scale_linetype_manual(values = rep(1, 2), labels = c("Theoretical", "Simulated")) + 
  scale_color_manual(values = hex, labels = c("Theoretical", "Simulated")) +
  ggtitle(paste0("Iterations vs Variance of method of moments estimator")) 
plot

df.var <- data.frame(id = rep(2:nsim, 2), 
                 var = c(rep(Final_var, (nsim - 1)), simvar * c(2:nsim)),
                 Type = as.factor(rep(1:2, each = (nsim - 1))))

plot <- ggplot(df.var, aes(x = id, group = Type)) +  # Plot for lower bound
  geom_line(aes(y = var, col = Type, linetype = Type), lwd=0.75) +
  labs(x = "Iterations", y = "Variance of method of moments estimator") +
  scale_linetype_manual(values = rep(1, 2), labels = c("Theoretical", "Simulated")) + 
  scale_color_manual(values = hex, labels = c("Theoretical", "Simulated")) +
  ggtitle(paste0("Iterations vs Variance of method of moments estimator")) 
plot

## Batch Means Estimator
phi <- 0.9
sd <- sqrt(0.5^2 / (1 - 0.9^2) * (1 - phi^2))
nsim <- 5000
nburnin <- 5000
data <- ts_sampler(niter=nsim+nburnin, nchain=1, mu=0, sd=sd, phi=phi, seed=1234)
data <- rbind(data,data^2)[,(nburnin+1):(nsim+nburnin)]

b <- 20
m <- nsim/b

muk <- matrix(0, nrow=2, ncol=m)
for(i in 1:m){
  for(j in 1:2){
    muk[j,i] <- mean(data[j,c(1:b)+(i-1)*b])
  }
}

mun <- apply(data, 1, mean)

temp <- 0
for(i in 1:m){
  temp <- temp + (muk[,i]-mun) %*% t(muk[,i]-mun)
}
Sigma <- b / (m-1) * temp
Sigma
