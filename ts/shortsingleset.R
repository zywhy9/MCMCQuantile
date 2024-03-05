source("../function/function.R")
library(ggplot2)
library(scales)
library(Matrix)
library(rstan)
library(ggpubr)
library(reshape2)

niter <- 100000 ## Number of iterations
nset <- 100 ## Number of simulated sets
nchain <- 1
npar <- 1
mu <- 0 
initial <- 0

phi <- 0.9995

true_sd <- sqrt(0.5^2 / (1 - 0.9^2))
true_quantile <-  qnorm(c(0.025, 0.975), mean = mu, sd = true_sd) ## True empirical quantile


n_burnin <- 70
par <- 1
n_est <- 3
hex <- hue_pal()(n_est)

## Absolute Difference

plots <- list()

abslow <- low <- absupp <- upp <- matrix(0, ncol = n_est, nrow = niter)
abslowmle <- abslowmm <- absloweq <- absuppmle <- absuppmm <- absuppeq <- matrix(0, ncol = nset, nrow = niter)

## Calculate average values
temp <- readRDS(paste0("res9995/res1.rds")) ## Read Data
actual_iter_mmeq <- temp$actual_iter_mm  ## Actual iterations used for MM and EQ
actual_niter_mmeq <- length(actual_iter_mmeq) ## Actual number of iterations used for MM and EQ
actual_iter_mle <- temp$actual_iter_mle ## Actual iterations used for MLE
actual_niter_mle <- length(actual_iter_mle) ## Actual number of iterations used for MLE

final_id_mmeq <- which(actual_iter_mmeq>n_burnin) ## Index of the iterations larger than burn-in for MM and EQ
final_id_mle <- which(actual_iter_mle>n_burnin) ## Index of the iterations larger than burn-in for MLE
final_iter_mmeq <- actual_iter_mmeq[final_id_mmeq]  ## Values of iterations larger than burn-in for MM and EQ
final_iter_mle <- actual_iter_mle[final_id_mle] ## Values of iterations larger than burn-in for MLE
final_niter_mmeq <- length(final_iter_mmeq) ## Number of iterations larger than burn-in for MM and EQ
final_niter_mle <- length(final_iter_mle) ## Number of iterations larger than burn-in for MLE

for(set in 1:nset){
  temp <- readRDS(paste0("res9995/res",set,".rds")) ## Read Data
  
  low_mle <- temp$low_mle[final_id_mle, par]
  low_mm <- temp$low_mm[final_id_mmeq, par]
  low_eq <- temp$low_eq[final_id_mmeq, par]
  upp_mle <- temp$upp_mle[final_id_mle, par]
  upp_mm <- temp$upp_mm[final_id_mmeq, par]
  upp_eq <- temp$upp_eq[final_id_mmeq, par]
  
  abslowmle[final_iter_mle, set] <- abs(low_mle - true_quantile[1])
  abslowmm[final_iter_mmeq, set] <- abs(low_mm - true_quantile[1])
  absloweq[final_iter_mmeq, set] <- abs(low_eq - true_quantile[1])
  absuppmle[final_iter_mle, set] <- abs(upp_mle - true_quantile[2])
  absuppmm[final_iter_mmeq, set] <- abs(upp_mm - true_quantile[2])
  absuppeq[final_iter_mmeq, set] <- abs(upp_eq - true_quantile[2])
  
  abslow[final_iter_mle, 1] <- abslow[final_iter_mle, 1] + 
    ifelse(is.na(abslowmle[final_iter_mle, set]), 0, abslowmle[final_iter_mle, set]) ## Absolute difference of lower bound for MLE
  abslow[final_iter_mmeq, 2] <- abslow[final_iter_mmeq, 2] + 
    ifelse(is.na(abslowmm[final_iter_mmeq, set]), 0, abslowmm[final_iter_mmeq, set]) ## Absolute difference of lower bound for MM
  abslow[final_iter_mmeq, 3] <- abslow[final_iter_mmeq, 3] + 
    ifelse(is.na(absloweq[final_iter_mmeq, set]), 0, absloweq[final_iter_mmeq, set]) ## Absolute difference of lower bound for EQ
  absupp[final_iter_mle, 1] <- absupp[final_iter_mle, 1] + 
    ifelse(is.na(absuppmle[final_iter_mle, set]), 0, absuppmle[final_iter_mle, set]) ## Absolute difference of upper bound for MLE
  absupp[final_iter_mmeq, 2] <- absupp[final_iter_mmeq, 2] + 
    ifelse(is.na(absuppmm[final_iter_mmeq, set]), 0, absuppmm[final_iter_mmeq, set]) ## Absolute difference of upper bound for MM
  absupp[final_iter_mmeq, 3] <- absupp[final_iter_mmeq, 3] + 
    ifelse(is.na(absuppeq[final_iter_mmeq, set]), 0, absuppeq[final_iter_mmeq, set]) ## Absolute difference of upper bound for EQ
  low[final_iter_mle, 1] <- low[final_iter_mle, 1] + ifelse(is.na(low_mle), 0, low_mle) ## Lower bound for MLE
  low[final_iter_mmeq, 2] <- low[final_iter_mmeq, 2] + ifelse(is.na(low_mm), 0, low_mm) ## Lower bound for MM
  low[final_iter_mmeq, 3] <- low[final_iter_mmeq, 3] + ifelse(is.na(low_eq), 0, low_eq) ## Lower bound for EQ
  upp[final_iter_mle, 1] <- upp[final_iter_mle, 1] + ifelse(is.na(upp_mle), 0, upp_mle) ## Upper bound for MLE
  upp[final_iter_mmeq, 2] <- upp[final_iter_mmeq, 2] + ifelse(is.na(upp_mm), 0, upp_mm) ## Upper bound for MM
  upp[final_iter_mmeq, 3] <- upp[final_iter_mmeq, 3] + ifelse(is.na(upp_eq), 0, upp_eq) ## Upper bound for EQ
}

abslow_avail_nset_mle <- apply(abslowmle,1, function(x){sum(!is.na(x))}) ## Remove the set with NA for each iteration
absupp_avail_nset_mle <- apply(absuppmle,1, function(x){sum(!is.na(x))}) ## Remove the set with NA for each iteration
abslow <- abslow / cbind(abslow_avail_nset_mle, rep(nset, niter), rep(nset, niter)) ## Average absolute difference of lower bound
absupp <- absupp / cbind(absupp_avail_nset_mle, rep(nset, niter), rep(nset, niter)) ## Average absolute difference of upper bound
low <- low / cbind(abslow_avail_nset_mle, rep(nset, niter), rep(nset, niter)) ## Average lower bound
upp <- upp / cbind(absupp_avail_nset_mle, rep(nset, niter), rep(nset, niter)) ## Average upper bound

df_abs <- data.frame(id = c(final_iter_mle, final_iter_mmeq, final_iter_mmeq), 
                     low = c(abslow[which(abslow[,1]!=0),1], abslow[which(abslow[,2]!=0),2], abslow[which(abslow[,3]!=0),3]),
                     upp = c(absupp[which(absupp[,1]!=0),1], absupp[which(absupp[,2]!=0),2], absupp[which(absupp[,3]!=0),3]),
                     Type = as.factor(c(rep(1, final_niter_mle), rep(2:3, each = final_niter_mmeq))))
df_abs$avg <- (df_abs$low+df_abs$upp)/2

plot2 <- ggplot(df_abs, aes(x = id, group = Type)) +  # Plot for lower bound
  geom_line(aes(y = low, col = Type, linetype = Type), lwd=0.75) +
  geom_point(aes(y = low, col = Type)) +
  scale_x_continuous(trans="log10") +
  scale_y_continuous(trans="log10", limits = c(0.23, 2.5)) +
  labs(x = "Iterations", y = "Absolute difference from true lower bound") +
  scale_linetype_manual(values = rep(1, n_est), labels = c("MLE", "MM", "Empirical Quantile")) + 
  scale_color_manual(values = hex, labels = c("MLE", "MM", "Empirical Quantile")) +
  ggtitle(paste0("Absolute difference from true lower bound")) +
  theme(legend.position = "none", axis.title.y=element_blank())  

plot3 <- ggplot(df_abs, aes(x = id, group = Type)) +  # Plot for upper bound
  geom_line(aes(y = upp, col = Type, linetype = Type), lwd=0.75) +
  geom_point(aes(y = upp, col = Type)) +
  scale_x_continuous(trans="log10") +
  scale_y_continuous(trans="log10", limits = c(0.23, 2.5)) +
  labs(x = "Iterations", y = "Absolute difference from true upper bound") +
  scale_linetype_manual(values = rep(1, n_est), labels = c("MLE", "MM", "Empirical Quantile")) + 
  scale_color_manual(values = hex, labels = c("MLE", "MM", "Empirical Quantile")) +
  ggtitle(paste0("Absolute difference from true upper bound")) +
  theme(legend.position = "none", 
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank())  

plot4 <- ggplot(df_abs, aes(x = id, group = Type)) +  # Plot for average difference
  geom_line(aes(y = avg, col = Type, linetype = Type), lwd=0.75) +
  geom_point(aes(y = avg, col = Type)) +
  scale_x_continuous(trans="log10") +
  scale_y_continuous(trans="log10", limits = c(0.23, 2.5)) +
  labs(x = "Iterations", y = "Average absolute difference") +
  scale_linetype_manual(values = rep(1, n_est), labels = c("MLE", "MM", "Empirical Quantile")) + 
  scale_color_manual(values = hex, labels = c("MLE", "MM", "Empirical Quantile")) +
  ggtitle(paste0("Average absolute difference")) +
  theme(axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank())

df <- data.frame(id = c(final_iter_mle, final_iter_mmeq, final_iter_mmeq, (n_burnin + 1):niter), 
                 low = c(low[which(low[,1]!=0),1], low[which(low[,2]!=0),2], low[which(low[,3]!=0),3],
                         rep(true_quantile[1], niter - n_burnin)),
                 upp = c(upp[which(upp[,1]!=0),1], upp[which(upp[,2]!=0),2], upp[which(upp[,3]!=0),3],
                         rep(true_quantile[2], niter - n_burnin)),
                 Type = as.factor(c(rep(0, final_niter_mle), rep(1:2, each = final_niter_mmeq), rep(3, niter - n_burnin))))

plot1 <- ggplot(df, aes(x = id, group = Type)) + # Full interval
  geom_line(aes(y = low, col = Type, linetype = Type), lwd=0.75) +
  geom_line(aes(y = upp, col = Type, linetype = Type), lwd=0.75) +
  geom_point(aes(y = low, col = Type)) +
  geom_point(aes(y = upp, col = Type)) +
  labs(x = "Iterations", y = "Intervals") +
  # ylim(-3, 2.5) + 
  scale_linetype_manual(values = c(1, 1, 1, 2), labels = c("MLE", "MM", "Empirical Quantile", "True Interval")) + 
  scale_color_manual(values = c(hex, "black"), labels = c("MLE", "MM", "Empirical Quantile", "True Interval")) +
  ggtitle(paste0("95% intervals when phi=", phi)) +
  theme(legend.position = "none")  

plotnames <- paste0("plot",1:4)
for(j in 1:4){
  assign(plotnames[j], get(paste0("plot",j)))
}
newlist <- mget(plotnames)
plots <- c(plots, newlist)

plot.list <- ggarrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]],
                       nrow = 1, ncol = 4, common.legend = TRUE, legend="bottom")
ggexport(plot.list, filename = "plot9995.pdf", width = 30, height = 5)
