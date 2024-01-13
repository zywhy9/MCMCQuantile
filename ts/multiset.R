source("../function/function.R")
library(ggplot2)
library(scales)
library(Matrix)
library(rstan)
library(ggpubr)
library(reshape2)

niter <- 5000 ## Number of iterations
nset <- 100 ## Number of simulated sets
nchain <- 1
npar <- 1
mu <- 0 
initial <- 0

phi <- c(0.4, 0.6, 0.9)

true.sd <- sqrt(0.5^2 / (1 - 0.9^2))
true.quantile <-  qnorm(c(0.025, 0.975), mean = mu, sd = true.sd) ## True empirical quantile


n.burnin <- 70
end.iter <- niter
final.iter <- end.iter - n.burnin
par <- 1
n.est <- 3
hex <- hue_pal()(n.est)

## Absolute Difference

plots <- list()
for(i in phi){
  
  abslow <- low <- absupp <- upp <- maxlow <- maxupp <- minlow <- minupp <- matrix(0, ncol = n.est, nrow = final.iter)
  abslowmle <- abslowmm <- absloweq <- absuppmle <- absuppmm <- absuppeq <- matrix(0, ncol = nset, nrow = final.iter)
  ## Calculate average values
  for(set in 1:nset){
    temp <- readRDS(paste0("res",(i*10),"/res",set,".rds")) ## Read Data
    low_mle <- temp$low_mle[(n.burnin + 1):end.iter, par]
    low_mm <- temp$low_mm[(n.burnin + 1):end.iter, par]
    low_eq <- temp$low_eq[(n.burnin + 1):end.iter, par]
    upp_mle <- temp$upp_mle[(n.burnin + 1):end.iter, par]
    upp_mm <- temp$upp_mm[(n.burnin + 1):end.iter, par]
    upp_eq <- temp$upp_eq[(n.burnin + 1):end.iter, par]
    abslowmle[,set] <- abs(low_mle - true.quantile[1])
    abslowmm[,set] <- abs(low_mm - true.quantile[1])
    absloweq[,set] <- abs(low_eq - true.quantile[1])
    absuppmle[,set] <- abs(upp_mle - true.quantile[2])
    absuppmm[,set] <- abs(upp_mm - true.quantile[2])
    absuppeq[,set] <- abs(upp_eq - true.quantile[2])
    
    abslow[,1] <- abslow[,1] + abslowmle[,set] ## Absolute difference of lower bound for MLE
    abslow[,2] <- abslow[,2] + abslowmm[,set] ## Absolute difference of lower bound for MM
    abslow[,3] <- abslow[,3] + absloweq[,set] ## Absolute difference of lower bound for EQ
    absupp[,1] <- absupp[,1] + absuppmle[,set]## Absolute difference of upper bound for MLE
    absupp[,2] <- absupp[,2] + absuppmm[,set] ## Absolute difference of upper bound for MM
    absupp[,3] <- absupp[,3] + absuppeq[,set] ## Absolute difference of upper bound for EQ
    low[,1] <- low[,1] + low_mle ## Lower bound for MLE
    low[,2] <- low[,2] + low_mm ## Lower bound for MM
    low[,3] <- low[,3] + low_eq ## Lower bound for EQ
    upp[,1] <- upp[,1] + upp_mle ## Upper bound for MLE
    upp[,2] <- upp[,2] + upp_mm ## Upper bound for MM
    upp[,3] <- upp[,3] + upp_eq ## Upper bound for EQ
  }
  
  abslow <- abslow / nset ## Average absolute difference of lower bound
  absupp <- absupp / nset ## Average absolute difference of upper bound
  low <- low / nset ## Average lower bound
  upp <- upp / nset ## Average upper bound
  
  df.abs <- data.frame(id = rep((n.burnin + 1):end.iter, n.est), 
                       low = as.vector(abslow),
                       upp = as.vector(absupp),
                       Type = as.factor(rep(1:n.est, each = final.iter)))
  df.abs$avg <- (df.abs$low+df.abs$upp)/2
  
  plot2 <- ggplot(df.abs, aes(x = id, group = Type)) +  # Plot for lower bound
    geom_line(aes(y = low, col = Type, linetype = Type), lwd=0.75) +
    scale_x_continuous(trans="log10") +
    scale_y_continuous(trans="log10", limits = c(0.02,0.83)) +
    labs(x = "Iterations", y = "Absolute difference from true lower bound") +
    scale_linetype_manual(values = rep(1, n.est), labels = c("MLE", "MM", "Empirical Quantile")) + 
    scale_color_manual(values = hex, labels = c("MLE", "MM", "Empirical Quantile")) +
    ggtitle(paste0("Absolute difference from true lower bound")) +
    theme(legend.position = "none", axis.title.y=element_blank())  
  
  plot3 <- ggplot(df.abs, aes(x = id, group = Type)) +  # Plot for upper bound
    geom_line(aes(y = upp, col = Type, linetype = Type), lwd=0.75) +
    scale_x_continuous(trans="log10") +
    scale_y_continuous(trans="log10", limits = c(0.02,0.83)) +
    labs(x = "Iterations", y = "Absolute difference from true upper bound") +
    scale_linetype_manual(values = rep(1, n.est), labels = c("MLE", "MM", "Empirical Quantile")) + 
    scale_color_manual(values = hex, labels = c("MLE", "MM", "Empirical Quantile")) +
    ggtitle(paste0("Absolute difference from true upper bound")) +
    theme(legend.position = "none", 
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title.y=element_blank())  
  
  plot4 <- ggplot(df.abs, aes(x = id, group = Type)) +  # Plot for average difference
    geom_line(aes(y = avg, col = Type, linetype = Type), lwd=0.75) +
    scale_x_continuous(trans="log10") +
    scale_y_continuous(trans="log10", limits = c(0.02,0.83)) +
    labs(x = "Iterations", y = "Average absolute difference") +
    scale_linetype_manual(values = rep(1, n.est), labels = c("MLE", "MM", "Empirical Quantile")) + 
    scale_color_manual(values = hex, labels = c("MLE", "MM", "Empirical Quantile")) +
    ggtitle(paste0("Average absolute difference")) +
    theme(axis.text.y=element_blank(), 
          axis.ticks.y=element_blank(),
          axis.title.y=element_blank())
  
  df <- data.frame(id = rep((n.burnin + 1):end.iter, (n.est + 1)), 
                   low = c(as.vector(low),
                           rep(true.quantile[1], final.iter)),
                   upp = c(as.vector(upp),
                           rep(true.quantile[2], final.iter)),
                   Type = as.factor(rep(0:n.est, each = final.iter)))
  
  plot1 <- ggplot(df, aes(x = id, group = Type)) + # Full interval
    geom_line(aes(y = low, col = Type, linetype = Type), lwd=0.75) +
    geom_line(aes(y = upp, col = Type, linetype = Type), lwd=0.75) +
    labs(x = "Iterations", y = "Intervals") +
    ylim(-3, 2.5) + 
    scale_linetype_manual(values = c(1, 1, 1, 2), labels = c("MLE", "MM", "Empirical Quantile", "True Interval")) + 
    scale_color_manual(values = c(hex, "black"), labels = c("MLE", "MM", "Empirical Quantile", "True Interval")) +
    ggtitle(paste0("95% intervals when phi=", i)) +
    theme(legend.position = "none")  
  
  plotnames <- paste0("plot",1:4,i*10)
  for(j in 1:4){
    assign(plotnames[j], get(paste0("plot",j)))
  }
  newlist <- mget(plotnames)
  plots <- c(plots, newlist)
}


plot.list <- ggarrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]],
                       plots[[5]], plots[[6]], plots[[7]], plots[[8]],
                       plots[[9]], plots[[10]], plots[[11]], plots[[12]],
                       nrow = 3, ncol = 4, common.legend = TRUE, legend="bottom")
ggexport(plot.list, filename = "plot.pdf", width = 30, height = 14)


## Variance

low_mle <- upp_mle <- low_mm <- upp_mm <- low_eq <- upp_eq <- matrix(0, ncol = nset * length(phi), nrow = final.iter)

for(i in 1:length(phi)){ ## Read Data

  for(set in 1:nset){
    temp <- readRDS(paste0("res",(phi[i]*10),"/res",set,".rds"))
    setno <- set + (i - 1) * nset
    low_mle[,setno] <- temp$low_mle[(n.burnin + 1):end.iter, par]
    low_mm[,setno] <- temp$low_mm[(n.burnin + 1):end.iter, par]
    low_eq[,setno] <- temp$low_eq[(n.burnin + 1):end.iter, par]
    upp_mle[,setno] <- temp$upp_mle[(n.burnin + 1):end.iter, par]
    upp_mm[,setno] <- temp$upp_mm[(n.burnin + 1):end.iter, par]
    upp_eq[,setno] <- temp$upp_eq[(n.burnin + 1):end.iter, par]
  }

}

low_mle_var <- upp_mle_var <- 
  low_mm_var <- upp_mm_var <- 
  low_eq_var <- upp_eq_var <- matrix(0, ncol = length(phi), nrow = end.iter - n.burnin)
for(i in 1:length(phi)){ ## Calculate variance
  setno <- 1:nset + (i - 1) * nset
  low_mle_var[,i] <- apply(low_mle[,setno], 1, var)
  low_mm_var[,i] <- apply(low_mm[,setno], 1, var)
  low_eq_var[,i] <- apply(low_eq[,setno], 1, var)
  upp_mle_var[,i] <- apply(upp_mle[,setno], 1, var)
  upp_mm_var[,i] <- apply(upp_mm[,setno], 1, var)
  upp_eq_var[,i] <- apply(upp_eq[,setno], 1, var)
}

# plots <- list()
# for(i in 1:length(phi)){
#   df.var <- data.frame(id = rep((n.burnin + 1):end.iter, n.est),
#                        low = c(low_mle_var[,i],
#                                low_mm_var[,i],
#                                low_eq_var[,i]),
#                        upp = c(upp_mle_var[,i],
#                                upp_mm_var[,i],
#                                upp_eq_var[,i]),
#                       Type = as.factor(rep(1:n.est, each = final.iter)))
# 
#   plot1 <- ggplot(df.var, aes(x = id, group = Type)) +
#     geom_line(aes(y = low, col = Type, linetype = Type), lwd=0.75) +
#     scale_x_continuous(trans = "log10") +
#     scale_y_continuous(trans = "log10", limits = c(0.001, 0.7))+
#     labs(x = "Iteration", y = "Variance of lower bound") +
#     scale_linetype_manual(values = c(1, 1, 1), labels = c("MLE", "MM", "Empirical Quantile")) +
#     scale_color_manual(values = hex, labels = c("MLE", "MM", "Empirical Quantile")) +
#     ggtitle(paste0("Variance of lower bound when phi=", phi[i]))
# 
# 
#   plot2 <- ggplot(df.var, aes(x = id, group = Type)) +
#     geom_line(aes(y = upp, col = Type, linetype = Type), lwd=0.75) +
#     scale_x_continuous(trans = "log10") +
#     scale_y_continuous(trans = "log10", limits = c(0.001, 2.5))+
#     labs(x = "Iteration", y = "Variance of upper bound") +
#     scale_linetype_manual(values = c(1, 1, 1), labels = c("MLE", "MM", "Empirical Quantile")) +
#     scale_color_manual(values = hex, labels = c("MLE", "MM", "Empirical Quantile")) +
#     ggtitle(paste0("Variance of upper bound when phi=", phi[i]))
# 
#   if(i==1){
#     plot1 <- plot1 +
#       theme(legend.position = "none")
#     plot2 <- plot2 +
#       theme(legend.position = "none")
#   }else{
#     plot1 <- plot1 +
#       theme(legend.position = "none",
#             axis.text.y=element_blank(),
#             axis.ticks.y=element_blank(),
#             axis.title.y=element_blank())
#     plot2 <- plot2 +
#       theme(legend.position = "none",
#             axis.text.y=element_blank(),
#             axis.ticks.y=element_blank(),
#             axis.title.y=element_blank())
#   }
# 
#   plotnames <- paste0("plot",1:2,phi[i]*10)
#   for(j in 1:2){
#     assign(plotnames[j], get(paste0("plot",j)))
#   }
#   newlist <- mget(plotnames)
#   plots <- c(plots, newlist)
# }
# plot.list <- ggarrange(plots[[1]], plots[[3]], plots[[5]],
#                        plots[[2]], plots[[4]], plots[[6]],
#                        nrow = 2, ncol = 3, common.legend = TRUE, legend="bottom")
# 
# ggexport(plot.list, filename = "plotvar.pdf", width = 20, height = 10)

low_var <- rbind(low_mle_var, low_mm_var, low_eq_var)
upp_var <- rbind(upp_mle_var, upp_mm_var, upp_eq_var)
df.var <- data.frame(id = rep((n.burnin + 1):end.iter, n.est * length(phi)), 
                     low = as.vector(low_var),
                     upp = as.vector(upp_var),
                     phi = rep(phi, each = n.est * final.iter),
                     Type = as.factor(rep(rep(1:n.est, each = final.iter), length(phi))))

df.var <- melt(df.var, id.vars = c("id", "phi", "Type"), variable.name = "VarType")

VarType.labs <- c("Lower Bound", "Upper Bound")
names(VarType.labs) <- c("low", "upp")

plotvar <- ggplot(df.var, aes(x = id, y = value)) +
  geom_line(aes(col = Type, linetype = Type), lwd=0.75) +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  labs(x = "Iteration", y = "Variance") +
  scale_linetype_manual(values = c(1, 1, 1), labels = c("MLE", "MM", "Empirical Quantile")) + 
  scale_color_manual(values = hex, labels = c("MLE", "MM", "Empirical Quantile")) +
  facet_grid(vars(VarType), vars(phi), labeller = labeller(VarType = VarType.labs))

ggexport(plotvar, filename = "plotvar.pdf", width = 20, height = 10)



## Ribbon Plot

plots <- list()
for(i in phi){
  
  abslow <- low <- absupp <- upp <- maxlow <- maxupp <- minlow <- minupp <- matrix(0, ncol = n.est, nrow = final.iter)
  abslowmle <- abslowmm <- absloweq <- absuppmle <- absuppmm <- absuppeq <- matrix(0, ncol = nset, nrow = final.iter)
  ## Calculate average values
  for(set in 1:nset){
    temp <- readRDS(paste0("res",(i*10),"/res",set,".rds")) ## Read Data
    low_mle <- temp$low_mle[(n.burnin + 1):end.iter, par]
    low_mm <- temp$low_mm[(n.burnin + 1):end.iter, par]
    low_eq <- temp$low_eq[(n.burnin + 1):end.iter, par]
    upp_mle <- temp$upp_mle[(n.burnin + 1):end.iter, par]
    upp_mm <- temp$upp_mm[(n.burnin + 1):end.iter, par]
    upp_eq <- temp$upp_eq[(n.burnin + 1):end.iter, par]
    abslowmle[,set] <- abs(low_mle - true.quantile[1])
    abslowmm[,set] <- abs(low_mm - true.quantile[1])
    absloweq[,set] <- abs(low_eq - true.quantile[1])
    absuppmle[,set] <- abs(upp_mle - true.quantile[2])
    absuppmm[,set] <- abs(upp_mm - true.quantile[2])
    absuppeq[,set] <- abs(upp_eq - true.quantile[2])
    
    abslow[,1] <- abslow[,1] + abslowmle[,set] ## Absolute difference of lower bound for MLE
    abslow[,2] <- abslow[,2] + abslowmm[,set] ## Absolute difference of lower bound for MM
    abslow[,3] <- abslow[,3] + absloweq[,set] ## Absolute difference of lower bound for EQ
    absupp[,1] <- absupp[,1] + absuppmle[,set]## Absolute difference of upper bound for MLE
    absupp[,2] <- absupp[,2] + absuppmm[,set] ## Absolute difference of upper bound for MM
    absupp[,3] <- absupp[,3] + absuppeq[,set] ## Absolute difference of upper bound for EQ
    low[,1] <- low[,1] + low_mle ## Lower bound for MLE
    low[,2] <- low[,2] + low_mm ## Lower bound for MM
    low[,3] <- low[,3] + low_eq ## Lower bound for EQ
    upp[,1] <- upp[,1] + upp_mle ## Upper bound for MLE
    upp[,2] <- upp[,2] + upp_mm ## Upper bound for MM
    upp[,3] <- upp[,3] + upp_eq ## Upper bound for EQ
  }
  
  abslow <- abslow / nset ## Average absolute difference of lower bound
  absupp <- absupp / nset ## Average absolute difference of upper bound
  low <- low / nset ## Average lower bound
  upp <- upp / nset ## Average upper bound
  
  ## Find bounds of the ribbons
  for(i in 1:final.iter){
    maxlow[i,1] <- quantile(abslowmle[i,], 0.975)
    maxlow[i,2] <- quantile(abslowmm[i,], 0.975)
    maxlow[i,3] <- quantile(absloweq[i,], 0.975)
    maxupp[i,1] <- quantile(absuppmle[i,], 0.975)
    maxupp[i,2] <- quantile(absuppmm[i,], 0.975)
    maxupp[i,3] <- quantile(absuppeq[i,], 0.975)
    minlow[i,1] <- quantile(abslowmle[i,], 0.025)
    minlow[i,2] <- quantile(abslowmm[i,], 0.025)
    minlow[i,3] <- quantile(absloweq[i,], 0.025)
    minupp[i,1] <- quantile(absuppmle[i,], 0.025)
    minupp[i,2] <- quantile(absuppmm[i,], 0.025)
    minupp[i,3] <- quantile(absuppeq[i,], 0.025)
  }
  
  df.abs <- data.frame(id = rep((n.burnin + 1):end.iter, n.est), 
                       low = as.vector(abslow),
                       upp = as.vector(absupp),
                       maxlow = as.vector(maxlow),
                       maxupp = as.vector(maxupp),
                       minlow = as.vector(minlow),
                       minupp = as.vector(minupp),
                       Type = as.factor(rep(1:n.est, each = final.iter)))
  df.abs$avg <- (df.abs$low+df.abs$upp)/2
  
  plot2 <- ggplot(df.abs, aes(x = id, group = Type)) +  # Plot for lower bound
    geom_line(aes(y = low, col = Type, linetype = Type), lwd=0.75) +
    geom_ribbon(aes(ymin = minlow, ymax = maxlow, fill = Type, alpha=Type)) +
    scale_x_continuous(trans="log10") +
    scale_y_continuous(trans="log10", limits = c(0.0001,1.95)) +
    labs(x = "Iterations", y = "Absolute difference from true lower bound") +
    scale_linetype_manual(values = rep(1, n.est), labels = c("MLE", "MM", "Empirical Quantile")) + 
    scale_color_manual(values = hex, labels = c("MLE", "MM", "Empirical Quantile")) +
    scale_alpha_manual(values = c(0.3,0.3,0.3), labels = c("MLE", "MM", "Empirical Quantile")) +
    ggtitle(paste0("Absolute difference from true lower bound")) +
    theme(legend.position = "none", axis.title.y=element_blank())  
  
  plot3 <- ggplot(df.abs, aes(x = id, group = Type)) +  # Plot for upper bound
    geom_line(aes(y = upp, col = Type, linetype = Type), lwd=0.75) +
    geom_ribbon(aes(ymin = minupp, ymax = maxupp, fill = Type, alpha=Type)) +
    scale_x_continuous(trans="log10") +
    scale_y_continuous(trans="log10", limits = c(0.0001,1.95)) +
    labs(x = "Iterations", y = "Absolute difference from true upper bound") +
    # ylim(0.03, 0.35) +
    scale_linetype_manual(values = rep(1, n.est), labels = c("MLE", "MM", "Empirical Quantile")) + 
    scale_color_manual(values = hex, labels = c("MLE", "MM", "Empirical Quantile")) +
    scale_alpha_manual(values = c(0.3,0.3,0.3), labels = c("MLE", "MM", "Empirical Quantile")) +
    ggtitle(paste0("Absolute difference from true upper bound")) +
    theme(legend.position = "none", 
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title.y=element_blank())  
  
  plot4 <- ggplot(df.abs, aes(x = id, group = Type)) +  # Plot for average difference
    geom_line(aes(y = avg, col = Type, linetype = Type), lwd=0.75) +
    scale_x_continuous(trans="log10") +
    scale_y_continuous(trans="log10", limits = c(0.0001,1.95)) +
    labs(x = "Iterations", y = "Average absolute difference") +
    scale_linetype_manual(values = rep(1, n.est), labels = c("MLE", "MM", "Empirical Quantile")) + 
    scale_color_manual(values = hex, labels = c("MLE", "MM", "Empirical Quantile")) +
    ggtitle(paste0("Average absolute difference")) +
    theme(axis.text.y=element_blank(), 
          axis.ticks.y=element_blank(),
          axis.title.y=element_blank())
  
  df <- data.frame(id = rep((n.burnin + 1):end.iter, (n.est + 1)), 
                   low = c(as.vector(low),
                           rep(true.quantile[1], final.iter)),
                   upp = c(as.vector(upp),
                           rep(true.quantile[2], final.iter)),
                   Type = as.factor(rep(0:n.est, each = final.iter)))
  
  plot1 <- ggplot(df, aes(x = id, group = Type)) + # Full interval
    geom_line(aes(y = low, col = Type, linetype = Type), lwd=0.75) +
    geom_line(aes(y = upp, col = Type, linetype = Type), lwd=0.75) +
    labs(x = "Iterations", y = "Intervals") +
    ylim(-3, 2.5) + 
    scale_linetype_manual(values = c(1, 1, 1, 2), labels = c("MLE", "MM", "Empirical Quantile", "True Interval")) + 
    scale_color_manual(values = c(hex, "black"), labels = c("MLE", "MM", "Empirical Quantile", "True Interval")) +
    ggtitle(paste0("95% intervals when phi=", i)) +
    theme(legend.position = "none")  
  
  plotnames <- paste0("plot",1:4,i*10)
  for(j in 1:4){
    assign(plotnames[j], get(paste0("plot",j)))
  }
  newlist <- mget(plotnames)
  plots <- c(plots, newlist)
}


plot.list <- ggarrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]],
                       plots[[5]], plots[[6]], plots[[7]], plots[[8]],
                       plots[[9]], plots[[10]], plots[[11]], plots[[12]],
                       nrow = 3, ncol = 4, common.legend = TRUE, legend="bottom")
ggexport(plot.list, filename = "ribbonplot.pdf", width = 30, height = 14)
