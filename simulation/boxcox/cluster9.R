setwd("~/scratch/quantile/boxcox/")

source("~/scratch/quantile/function.R")
library(MASS)

niter <- 5000 ## Number of iterations
nchain <- 1
npar <- 1
mu <- 0
phi <- 0.9
sd <- sqrt(1 - phi ^ 2)
true_sd <- 1
initial <- 0

transf <- "no" ## Indicators of transformation

for(seed in 1:100) {
  data <- ts_sampler(niter, nchain, mu, sd, phi, initial, seed)
  data <- qgamma(pnorm(data, mean = 0, sd = true_sd), 2, 1)
  if(is.null(dim(data))) {data <- as.matrix(data)}
  
  ## Analysis
  ## Empirical Quantile
  c(low_eq, upp_eq, time_eq) %<-% eq_analysis(data=data, niter=niter, npar=npar, nchain=nchain)
  
  ## moment estimate
  b <- boxcox(lm(data ~ 1), plotit = F)
  lambda <- b$x[which.max(b$y)]
  data <- (data ^ lambda - 1) / lambda
  c(mean_mm, sd_mm, low_mm, upp_mm, time_mm) %<-% mm_analysis(data=data, niter=niter, npar=npar, nchain=nchain, transf=transf)
  low_mm <- (low_mm * lambda + 1) ^ (1 / lambda)
  upp_mm <- (upp_mm * lambda + 1) ^ (1 / lambda)
  
  savelist <- list(low_mm = low_mm, upp_mm = upp_mm, time_mm = time_mm,
                   low_eq = low_eq, upp_eq = upp_eq, time_eq = time_eq)
  
  saveRDS(savelist, paste0("res9/res",seed,".rds"))
  
}


