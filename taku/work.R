source("../function/function.R")
library(ggplot2)

niter <- 30000 ## Number of iterations
transf <- "no" ## Indicators of transformation
c(mcmc.samples, nchain, npar, param.name) %<-% mcmc_data(niter=niter, param.name="N") ## Read Data

## moment estimate
c(mean_mm, sd_mm, low_mm, upp_mm, time_mm) %<-% mm_analysis(data=mcmc.samples, niter=niter, npar=npar, nchain=nchain, transf=transf)

## Empirical Quantile
c(low_eq, upp_eq, time_eq) %<-% eq_analysis(data=mcmc.samples, niter=niter, npar=npar, nchain=nchain)

## MLE
set.seed(1234)
initial.p <- list(list("p"=rdirichlet(12,rep(1,4))),
                  list("p"=rdirichlet(12,rep(1,4))),
                  list("p"=rdirichlet(12,rep(1,4))))
data <- readRDS("sim/sim1.rds")
initial <- c()
for(i in 1:nchain){
  initial[i] <- data$Nw/sum(data$w*initial.p[[i]]$p[,2],data$w*initial.p[[i]]$p[,3])
}

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
full.samples.H <- full.samples[,,"N"]
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
  ylim(55000, 70000) + 
  scale_linetype_manual(values=c(1,1,1,2,3), labels = c("MLE", "MM", "Empirical Quantile", "True Interval", "NA Interval")) + 
  scale_color_manual(values=c("#F8766D","#00BA38","#619CFF","black","purple"), labels = c("MLE", "MM", "Empirical Quantile", "True Interval", "NA Interval"))
