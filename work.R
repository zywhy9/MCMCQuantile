source("../function/function.R")
library(ggplot2)

niter <- 30000 ## Number of iterations
transf <- c("log", "logit", "logit") ## Indicators of transformation
c(mcmc.samples, nchain, npar, param.name) %<-% mcmc_data(niter=niter) ## Read Data

## moment estimate
c(mean_mm, sd_mm, low_mm, upp_mm, time_mm) %<-% mm_analysis(data=mcmc.samples, niter=niter, npar=npar, nchain=nchain, transf=transf)

## Empirical Quantile
c(low_eq, upp_eq, time_eq) %<-% eq_analysis(data=mcmc.samples, niter=niter, npar=npar, nchain=nchain)

## MLE
simdata <- readRDS("sim/sim1.rds")
H.inits <- (simdata$census.unkn - round(simdata$M.maybe*0.7))/0.7

initial <- c(log(rep(H.inits, nchain)), boot::logit(c(0.5, 0.3, 0.1, 0.5, 0.6, 0.4)))
rm(simdata, H.inits)

# initial <- rep(c(log(1000), boot::logit(0.5), boot::logit(0.5)), each = nchain)
c(mean_mle, sd_mle, low_mle, upp_mle, time_mle, mu_mle, psi_mle, beta, gamma, delta) %<-% mle_analysis(data=mcmc.samples, 
                                                                                          niter=niter, 
                                                                                          npar=npar, 
                                                                                          nchain=nchain, 
                                                                                          transf=transf, 
                                                                                          initial = initial,
                                                                                          extra = T)
# save.image("full.RData")

#### Plot
full.samples <- readRDS("res/res1.rds")
full.samples.H <- full.samples[,,"H"]
log.full.samples.H <- log(full.samples.H)

true.quantile <- quantile(full.samples.H, c(0.025, 0.975)) ## True empirical quantile

true.mean <- mean(log.full.samples.H)
true.sd <- sd(log.full.samples.H)
na.quantile <- exp(c(true.mean - 1.96 * true.sd, true.mean + 1.96 * true.sd)) ## True normal-approximated quantile

n.burnin <- 20 ## Number of burn-ins
par <- 1 ## ID for the parameter of interest
df <- data.frame(id = rep((n.burnin+1):niter,5), 
                 low = c(low_mle[(n.burnin+1):niter,par], ## MLE
                         low_mm[(n.burnin+1):niter,par], ## MM
                         low_eq[(n.burnin+1):niter,par], ## Empirical
                         rep(true.quantile[1], (niter - n.burnin)), ## True empirical
                         rep(na.quantile[1], (niter - n.burnin))),  ## True NA
                 upp = c(upp_mle[(n.burnin+1):niter,par], ## MLE
                         upp_mm[(n.burnin+1):niter,par], ## MM
                         upp_eq[(n.burnin+1):niter,par], ## Empirical
                         rep(true.quantile[2], (niter - n.burnin)), ## True empirical
                         rep(na.quantile[2], (niter - n.burnin))),  ## True NA
                 Type = as.factor(rep(0:4, each = (niter - n.burnin))))

ggplot(df, aes(x = id, group = Type)) +
  geom_line(aes(y = low, col = Type, linetype = Type), lwd=0.75) +
  geom_line(aes(y = upp, col = Type, linetype = Type), lwd=0.75) +
  labs(x = "Iterations", y = "H") +
  ylim(1200, 2000) + 
  scale_linetype_manual(values=c(1,1,1,2,3), labels = c("MLE", "MM", "Empirical Quantile", "True Interval", "NA Interval")) + 
  scale_color_manual(values=c("#F8766D","#00BA38","#619CFF","black","purple"), labels = c("MLE", "MM", "Empirical Quantile", "True Interval", "NA Interval"))

#### Simulation using Markov-normal Model

library(MASS)
set.seed(1234)

## Simulate Data
res <- ar_simulation(niter=niter, npar=npar, nchain=nchain, transf=transf, init=initial, beta=beta, gamma=gamma, delta=delta)

# save(res, niter, nchain, npar, transf, initial, file = "simdata.RData")

## moment estimate
c(mean_mm, sd_mm, low_mm, upp_mm, time_mm) %<-% mm_analysis(data=res, niter=niter, npar=npar, nchain=nchain, transf=transf)

## Empirical Quantile
c(low_eq, upp_eq, time_eq) %<-% eq_analysis(data=res, niter=niter, npar=npar, nchain=nchain)

## MLE
c(mean_mle, sd_mle, low_mle, upp_mle, time_mle, mu_mle, psi_mle, beta, gamma, delta) %<-% mle_analysis(data=res,
                                                                                                       niter=1000, 
                                                                                                       npar=npar,
                                                                                                       nchain=nchain,
                                                                                                       transf=transf,
                                                                                                       initial = initial,
                                                                                                       extra = T,
                                                                                                       nburnin=15)

## Plot
full.samples.H <- res[,1:nchain]
log.full.samples.H <- log(full.samples.H)

true.quantile <- quantile(full.samples.H, c(0.025, 0.975)) ## True empirical quantile

true.mean <- mean(log.full.samples.H)
true.sd <- sd(log.full.samples.H)
na.quantile <- exp(c(true.mean - 1.96 * true.sd, true.mean + 1.96 * true.sd)) ## True normal-approximated quantile

n.burnin <- 20
par <- 1
df <- data.frame(id = rep((n.burnin+1):niter,5), 
                 low = c(low_mle[(n.burnin+1):niter,par], 
                         low_mm[(n.burnin+1):niter,par],
                         low_eq[(n.burnin+1):niter,par],
                         rep(true.quantile[1], (niter - n.burnin)),
                         rep(na.quantile[1], (niter - n.burnin))), 
                 upp = c(upp_mle[(n.burnin+1):niter,par], 
                         upp_mm[(n.burnin+1):niter,par], 
                         upp_eq[(n.burnin+1):niter,par],
                         rep(true.quantile[2], (niter - n.burnin)),
                         rep(na.quantile[2], (niter - n.burnin))), 
                 Type = as.factor(rep(0:4, each = (niter - n.burnin))))

ggplot(df, aes(x = id, group = Type)) +
  geom_line(aes(y = low, col = Type, linetype = Type), lwd=0.75) +
  geom_line(aes(y = upp, col = Type, linetype = Type), lwd=0.75) +
  labs(x = "Iterations", y = "H") +
  ylim(750, 3000) + 
  scale_linetype_manual(values=c(1,1,1,2,3), labels = c("MLE", "MM", "Empirical Quantile", "True Interval", "NA Interval")) + 
  scale_color_manual(values=c("#F8766D","#00BA38","#619CFF","black","purple"), labels = c("MLE", "MM", "Empirical Quantile", "True Interval", "NA Interval"))




