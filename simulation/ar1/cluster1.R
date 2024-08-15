setwd("~/scratch/quantile/ar1/")

source("~/scratch/quantile/function.R")
library(doParallel)

niter <- 5000 ## Number of iterations
nchain <- 1
npar <- 1
mu <- 0
phi <- 0.1
sd <- sqrt(1 - phi ^ 2)
true_sd <- 1
initial <- 0

transf <- "no" ## Indicators of transformation

ncores = Sys.getenv("SLURM_CPUS_PER_TASK") 
registerDoParallel(cores=ncores)# Shows the number of Parallel Workers to be used
print(ncores) # this how many cores are available, and how many you have requested.
getDoParWorkers()# you can compare with the number of actual workers

foreach(seed=1:100) %dopar% {
  data <- ts_sampler(niter, nchain, mu, sd, phi, initial, seed)
  if(is.null(dim(data))) {data <- as.matrix(data)}
  
  ## Analysis
  
  ## moment estimate
  c(mean_mm, sd_mm, low_mm, upp_mm, time_mm) %<-% mm_analysis(data=data, niter=niter, npar=npar, nchain=nchain, transf=transf)
  
  ## Empirical Quantile
  c(low_eq, upp_eq, time_eq) %<-% eq_analysis(data=data, niter=niter, npar=npar, nchain=nchain)
  
  
  savelist <- list(mean_mm = mean_mm, sd_mm = sd_mm, low_mm = low_mm, upp_mm = upp_mm, time_mm = time_mm,
                   low_eq = low_eq, upp_eq = upp_eq, time_eq = time_eq)
  
  saveRDS(savelist, paste0("res1/res",seed,".rds"))
  
}


