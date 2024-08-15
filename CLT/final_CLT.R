source("../function/function.R")
library(progress)
library(mcmcse)
library(ggplot2)
library(scales)
library(mvtnorm)
library(coda)
library(gridExtra)
library(ggpubr)

#### Simulation
nadapt <- 9  ## Avoid first NA in EQ
nburnin <- 100 ## Number of burn-in
nset <- 1000 ## Number of datasets
niter <- 30000 ## Number of iterations
nchain <- 1 ## Number of chains
npar <- 1 ## Number of parameters
mu <- 0 ## AR(1) mean
phi <- 0.9 ## AR(1) coefficient
sd <- sqrt(1 - phi ^ 2) ## AR(1) SD, which guarantees a N(0,1) stationary distribution
initial <- 0 ## Initial value for AR(1)

transf <- "no" ## Indicators of transformation

# if(!dir.exists("simdata/")){
#   dir.create("simdata/")
# }
# if(!dir.exists("CLT/")){
#   dir.create("CLT/")
# }

pb <- progress_bar$new(format = "[:bar] :current/:total (:percent)", total = nset) ## Set up prograss bar

for(seed in 1:nset){
  data <- ts_sampler(niter, nchain, mu, sd, phi, initial, seed) ## Simulate an AR(1) process
  if(is.null(dim(data))) {data <- as.matrix(data)}

  saveRDS(data, paste0("simdata/res",seed,".rds"))
  c(mean_mm, sd_mm, low_mm, upp_mm, time_mm) %<-% mm_analysis(data=data, niter=niter, npar=npar, nchain=nchain, transf=transf) ## Calculate MM estimator
  
  c(low_eq, upp_eq, time_eq) %<-% eq_analysis(data=data, niter=niter, npar=npar, nchain=nchain) ## Calculate EQ estimator
  
  savelist <- list(upp_mm=upp_mm, upp_eq=upp_eq)
  
  saveRDS(savelist, paste0("CLT/res",seed,".rds"))
  pb$tick(1) ## Progress bar
}


#### MM CLT Variance
true_sd <- 1 ## SD of stationary distribution
true_quantile <-
  qnorm(0.975, mean = mu, sd = true_sd) ## True posterior quantile
varx <- true_sd ^ 2 ## Variance of stationary distribution sigma^2
varx2 <- varx ^ 2 ## sigma^4
Sigma11 <- varx * (1 + phi) / (1 - phi)
Sigma22 <- 2 * varx2 * (1 + phi ^ 2) / (1 - phi ^ 2)
Final_var_MM <- Sigma11 + qnorm(0.975) ^ 2 / (4 * varx) * Sigma22
MMclt <- sqrt(Final_var_MM / c((nburnin + 1):niter))


#### EQ CLT Variance (Doss 2014)
covk <- function(k, phi, true_sd, true_quantile) {
  true_var <- true_sd ^ 2
  sigma <- matrix(c(true_var, phi ^ k * true_var, phi ^ k * true_var, true_var), 2, 2)
  pmvnorm(
    upper = rep(true_quantile, 2),
    mean = rep(0, 2),
    sigma = sigma,
    keepAttr = F
  ) - 0.975 ^ 2
}

temp <- 0
for (i in 1:100) {
  temp <- temp + covk(i, phi, true_sd, true_quantile)
}
Final_sigma2_EQ <- 0.975 * (1 - 0.975) + 2 * temp
fv <- dnorm(true_quantile, mean = 0, sd = true_sd)
Final_var_EQ <- Final_sigma2_EQ / (fv ^ 2)
EQclt <- sqrt(Final_var_EQ / c((nburnin + 1):niter))

#### True MCSE

MM <- matrix(NA, nrow = nset, ncol = niter)
EQ <- matrix(NA, nrow = nset, ncol = niter)

for (set in 1:nset) {
  temp <- readRDS(paste0("CLT/res", set, ".rds"))
  MM[set,] <- as.vector(temp$upp_mm)
  EQ[set,] <- as.vector(temp$upp_eq)
}
rm(temp)

MMtrue <- rep(NA, niter - nburnin)
EQtrue <- rep(NA, niter - nburnin)
for (iter in 1:(niter - nburnin)) {
  MMtrue[iter] <- sum((MM[, (nburnin + iter)] - true_quantile) ^  2) / (nset - 1)
  EQtrue[iter] <- sum((EQ[, (nburnin + iter)] - true_quantile) ^  2) / (nset - 1)
}
rm(MM)
rm(EQ)

#### Estimated CLT
# q <- qnorm(0.975)
# 
# res <- foreach(set = 1:nset, .combine = "rbind", .packages = "mcmcse") %dopar% {
#   temp_mm <- rep(NA, niter)
#   temp <- readRDS(paste0("simdata/res", set, ".rds"))
#   temp2 <- temp^2
#   mmtemp <- cbind(temp, temp2)
#   for (iter in 3:niter) {
  #   res <- mcse.multi(mmtemp[1:iter, ], method = "obm")
  #   Sigma <- res$cov
  #   mux <- res$est[1]
  #   sigmax <- sd(temp[1:iter])
  #   temp_mm[iter] <- (1 - q * mux / sigmax)^2 * Sigma[1, 1] + 
  #     (q / (2 * sigmax)) ^ 2 * Sigma[2, 2] +
  #     (q / sigmax) * (1 - q * mux / sigmax) * Sigma[1, 2]
  # }
#   temp_mm
# }
# res <- foreach(set = 1:nset, .combine = "rbind", .packages = "mcmcse") %dopar% {
#   temp_eq <- rep(NA, niter - 1)
#   temp <- readRDS(paste0("simdata/res", set, ".rds"))
#   
#   for (iter in 10:niter) {
#     temp_eq[(iter - 1)] <- mcse.q(temp[1:iter,1], q=0.975)$se
#   }
#   temp_eq
# }


load("mcseeq.RData")

eq_se <- res[, nburnin: (niter - 1)]
rm(res)

mm_var <- matrix(NA, nrow = niter - nburnin, ncol = nset)
for (set in 1:nset) {
  mm_var[, set] <- readRDS(paste0("MMCLT/res", set, ".rds"))[(nburnin + 1):niter]
}
mm_se <- t(apply(mm_var, 2, function(x){sqrt(x/c((nburnin+1): niter))}))
rm(mm_var)

#### Plots
hex <- hue_pal()(2)
id <- (nburnin + 1):niter
# mean_eq <- apply(eq_se, 2, mean)
# sigma_eq <- apply(eq_se, 2, sd)
# lb_eq <- mean_eq - sigma_eq
# ub_eq <- mean_eq + sigma_eq
hpd_eq <- apply(eq_se,2,function(x){HPDinterval(as.mcmc(x))})

# mean_mm <- apply(mm_se, 2, mean)
# sigma_mm <- apply(mm_se, 2, sd)
# lb_mm <- mean_mm - sigma_mm
# ub_mm <- mean_mm + sigma_mm
hpd_mm <- apply(mm_se,2,function(x){HPDinterval(as.mcmc(x))})

df <- data.frame(id = id,
                 # mu_mm = mean_mm,
                 # mu_eq = mean_eq,
                 # lb_mm = lb_mm,
                 # lb_eq = lb_eq,
                 # ub_mm = ub_mm,
                 # ub_eq = ub_eq,
                 true_mm = sqrt(MMtrue),
                 true_eq = sqrt(EQtrue),
                 clt_mm = MMclt,
                 clt_eq = EQclt)

thick <- 0.5

mmplot <- ggplot(df, aes(x = id)) +
  # geom_line(aes(y = mu_mm, linetype = "CLTest"), lwd = thick) +
  geom_line(aes(y = true_mm, color = "True", linetype = "True"), lwd = thick) +
  geom_line(aes(y = clt_mm, color = "CLT", linetype = "CLT"), lwd = thick) +
  geom_ribbon(aes(ymin=hpd_mm[1,], ymax=hpd_mm[2,]), alpha = 0.2) +
  labs(x = "Iterations", y = "Monte Carlo Standard Error", linetype = "Type", color = "Type") +
  theme_minimal() +
  theme(legend.position = "bottom")

eqplot <- ggplot(df, aes(x = id)) +
  # geom_line(aes(y = mu_eq, linetype = "CLTest"), lwd = thick) +
  geom_line(aes(y = true_eq, color = "True", linetype = "True"), lwd = thick) +
  geom_line(aes(y = clt_eq, color = "CLT", linetype = "CLT"), lwd = thick) +
  geom_ribbon(aes(ymin=hpd_eq[1,], ymax=hpd_eq[2,]), alpha = 0.2) +
  labs(x = "Iterations", y = "Monte Carlo Standard Error", linetype = "Type", color = "Type") +
  theme_minimal() +
  theme(legend.position = "bottom")

combplots <- ggarrange(mmplot, eqplot, ncol = 2, nrow = 1, common.legend = T, legend = "bottom")
combplots

ggplot(df, aes(x = id)) +
  geom_line(aes(y = true_eq, color = "EQ", linetype = "True"), lwd = thick) +
  geom_line(aes(y = clt_eq, color = "EQ", linetype = "CLT"), lwd = thick) +
  geom_line(aes(y = true_mm, color = "MM", linetype = "True"), lwd = thick) +
  geom_line(aes(y = clt_mm, color = "MM", linetype = "CLT"), lwd = thick) +
  geom_ribbon(aes(ymin=hpd_eq[1,], ymax=hpd_eq[2,]), fill = hex[1], alpha = 0.2) +
  geom_ribbon(aes(ymin=hpd_mm[1,], ymax=hpd_mm[2,]), fill = hex[2], alpha = 0.2) +
  labs(x = "Iterations", y = "Monte Carlo Standard Error", linetype = "Type", color = "Estimator") +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14))

