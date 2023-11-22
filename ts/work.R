source("function/function.R")
library(ggplot2)
library(Matrix)

## AR(1) Simulation

niter <- 10000 ## Number of iterations
nchain <- 3
mu <- 10 
sd <- 5 
phi <- 0.9
initial <- c(0,10,20)

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

data <- ts_sampler(niter, nchain, mu, sd, phi, initial)

df <- data.frame(id = rep(1:niter, nchain), x = as.vector(data), chain = as.factor(rep(1:3, each = niter)))
ggplot(data = df, aes(x = id, y = x, group = chain)) +
  geom_line(aes(color=chain)) +
  xlab("Iteration") +
  ylab("Sample") +
  scale_color_discrete(name = "Chain")

## Analysis
transf <- "no" ## Indicators of transformation

## moment estimate
c(mean_mm, sd_mm, low_mm, upp_mm, time_mm) %<-% mm_analysis(data=data, niter=niter, npar=1, nchain=nchain, transf=transf)

## Empirical Quantile
c(low_eq, upp_eq, time_eq) %<-% eq_analysis(data=data, niter=niter, npar=1, nchain=nchain)

## MLE
c(mean_mle, sd_mle, low_mle, upp_mle, time_mle, mu_mle, psi_mle, beta, gamma, delta) %<-% mle_analysis(data=data, 
                                                                                                       niter=niter, 
                                                                                                       npar=1, 
                                                                                                       nchain=nchain, 
                                                                                                       transf=transf, 
                                                                                                       initial = initial,
                                                                                                       extra = T,
                                                                                                       nburnin = 15)
