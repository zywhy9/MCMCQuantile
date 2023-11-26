source("../function/function.R")
library(ggplot2)
library(scales)
library(Matrix)
library(rstan)

## AR(1) Simulation

niter <- 30000 ## Number of iterations
nchain <- 3
mu <- 0 
sd <- 0.5 
initial <- c(-10,0,10)

phi <- 0.99
nburnin <- 40

phi <- 0.995
nburnin <- 40

phi <- 0.9995
nburnin <- 700

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

df <- data.frame(id = rep(1:niter, nchain), x = as.vector(data), chain = as.factor(rep(1:nchain, each = niter)))
ggplot(data = df, aes(x = id, y = x, group = chain)) +
  geom_line(aes(color=chain)) +
  xlab("Iterations") +
  ylab("Sample") +
  scale_color_discrete(name = "Chain")

## Analysis
transf <- "no" ## Indicators of transformation

## moment estimate
c(mean_mm, sd_mm, low_mm, upp_mm, time_mm) %<-% mm_analysis(data=data, niter=niter, npar=1, nchain=nchain, transf=transf)

## Empirical Quantile
c(low_eq, upp_eq, time_eq) %<-% eq_analysis(data=data, niter=niter, npar=1, nchain=nchain)

## Highest posterior density interval
c(low_hpd, upp_hpd, time_hpd) %<-% hpd_analysis(data=data, niter=niter, npar=1, nchain=nchain, transf=transf)

## MLE
c(mean_mle, sd_mle, low_mle, upp_mle, time_mle, mu_mle, psi_mle, beta, gamma, delta) %<-% mle_analysis(data=data, 
                                                                                                       niter=100, 
                                                                                                       npar=1, 
                                                                                                       nchain=nchain, 
                                                                                                       transf=transf, 
                                                                                                       initial = initial,
                                                                                                       extra = T,
                                                                                                       nburnin = nburnin)

## Plot
true.sd <- sqrt(sd^2 / (1 - phi^2))
true.quantile <-  qnorm(c(0.025, 0.975), mean = mu, sd = true.sd)## True empirical quantile

# true.mean <- mean(data)
# true.sd <- sd(data)
# na.quantile <- c(true.mean - 1.96 * true.sd, true.mean + 1.96 * true.sd) ## True normal-approximated quantile

n.burnin <- round(1.3 * nburnin)
# end.iter <- 2500
end.iter <- niter
par <- 1
n.est <- 4
hex <- hue_pal()(n.est)

df <- data.frame(id = rep((n.burnin + 1):end.iter, (n.est + 1)), 
                 low = c(low_mle[(n.burnin + 1):end.iter, par], 
                         low_mm[(n.burnin + 1):end.iter, par],
                         low_eq[(n.burnin + 1):end.iter, par],
                         low_hpd[(n.burnin + 1):end.iter, par],
                         rep(true.quantile[1], (end.iter - n.burnin))),
                 upp = c(upp_mle[(n.burnin + 1):end.iter, par], 
                         upp_mm[(n.burnin + 1):end.iter, par], 
                         upp_eq[(n.burnin + 1):end.iter, par],
                         upp_hpd[(n.burnin + 1):end.iter, par],
                         rep(true.quantile[2], (end.iter - n.burnin))),
                 Type = as.factor(rep(0:n.est, each = (end.iter - n.burnin))))

ggplot(df, aes(x = id, group = Type)) + # Full interval
  geom_line(aes(y = low, col = Type, linetype = Type), lwd=0.75) +
  geom_line(aes(y = upp, col = Type, linetype = Type), lwd=0.75) +
  labs(x = "Iterations", y = "X") +
  ylim(-10, 25) + 
  scale_linetype_manual(values = c(1, 1, 1, 1, 2), labels = c("MLE", "MM", "Empirical Quantile", "HPD", "True Interval")) + 
  scale_color_manual(values = c(hex, "black"), labels = c("MLE", "MM", "Empirical Quantile", "HPD", "True Interval"))

ggplot(df, aes(x = id, group = Type)) +  # Only upper bound
  geom_line(aes(y = upp, col = Type, linetype = Type), lwd=0.75) +
  labs(x = "Iterations", y = "Upper Bound") +
  ylim(floor(min(df$upp)), ceiling(max(df$upp))) + 
  scale_linetype_manual(values = c(1, 1, 1, 1, 2), labels = c("MLE", "MM", "Empirical Quantile", "HPD", "True Interval")) + 
  scale_color_manual(values = c(hex, "black"), labels = c("MLE", "MM", "Empirical Quantile", "HPD", "True Interval"))

ggplot(df, aes(x = id, group = Type)) +  # Only lower bound
  geom_line(aes(y = low, col = Type, linetype = Type), lwd=0.75) +
  labs(x = "Iterations", y = "Lower Bound") +
  ylim(floor(min(df$low)), ceiling(max(df$low))) + 
  scale_linetype_manual(values = c(1, 1, 1, 1, 2), labels = c("MLE", "MM", "Empirical Quantile", "HPD", "True Interval")) + 
  scale_color_manual(values = c(hex, "black"), labels = c("MLE", "MM", "Empirical Quantile", "HPD", "True Interval"))

## Check absolute difference with true CI

n.burnin <- round(1.3 * nburnin)
end.iter <- niter
n.est <- 4
hex <- hue_pal()(n.est)

df <- data.frame(id = rep((n.burnin + 1):end.iter, n.est), 
                 low = c(low_mle[(n.burnin + 1):end.iter, par]-true.quantile[1], 
                         low_mm[(n.burnin + 1):end.iter, par]-true.quantile[1],
                         low_eq[(n.burnin + 1):end.iter, par]-true.quantile[1],
                         low_hpd[(n.burnin + 1):end.iter, par]-true.quantile[1]),
                 upp = c(upp_mle[(n.burnin + 1):end.iter, par]-true.quantile[2], 
                         upp_mm[(n.burnin + 1):end.iter, par]-true.quantile[2],
                         upp_eq[(n.burnin + 1):end.iter, par]-true.quantile[2],
                         upp_hpd[(n.burnin + 1):end.iter, par]-true.quantile[2]),
                 Type = as.factor(rep(1:n.est, each = (end.iter - n.burnin))))

ggplot(df, aes(x = id, group = Type)) +  # Plot for lower bound
  geom_line(aes(y = abs(low), col = Type, linetype = Type), lwd=0.75) +
  labs(x = "Iterations", y = "Absolute difference from true lower bound") +
  ylim(floor(min(abs(df$low))), ceiling(max(abs(df$low)))) + 
  scale_linetype_manual(values = c(1, 1, 1, 1), labels = c("MLE", "MM", "Empirical Quantile", "HPD")) + 
  scale_color_manual(values = hex, labels = c("MLE", "MM", "Empirical Quantile", "HPD"))

ggplot(df, aes(x = id, group = Type)) +  # Plot for upper bound
  geom_line(aes(y = abs(upp), col = Type, linetype = Type), lwd=0.75) +
  labs(x = "Iterations", y = "Absolute difference from true upper bound") +
  ylim(floor(min(abs(df$upp))), ceiling(max(abs(df$upp)))) + 
  scale_linetype_manual(values = c(1, 1, 1, 1), labels = c("MLE", "MM", "Empirical Quantile", "HPD")) + 
  scale_color_manual(values = hex, labels = c("MLE", "MM", "Empirical Quantile", "HPD"))

## ESS

ess <- rep(0, niter)

for(i in 2:niter){
  ess[i] <- ess_bulk(data[1:i,])
}










