source("../function/function.R")
library(ggplot2)
library(scales)
library(Matrix)
library(rstan)
library(ggpubr)

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
par <- 1
n.est <- 3
hex <- hue_pal()(n.est)

## Absolute Difference

plots <- list()
for(i in phi){
  
  abslow <- low <- absupp <- upp <- matrix(0, ncol = n.est, nrow = end.iter - n.burnin)
  ## Calculate average values
  for(set in 1:nset){
    temp <- readRDS(paste0("res",(i*10),"/res",set,".rds")) ## Read Data
    abslow[,1] <- abslow[,1] + abs(temp$low_mle[(n.burnin + 1):end.iter, par]-true.quantile[1])
    abslow[,2] <- abslow[,2] + abs(temp$low_mm[(n.burnin + 1):end.iter, par]-true.quantile[1])
    abslow[,3] <- abslow[,3] + abs(temp$low_eq[(n.burnin + 1):end.iter, par]-true.quantile[1])
    absupp[,1] <- absupp[,1] + abs(temp$upp_mle[(n.burnin + 1):end.iter, par]-true.quantile[2])
    absupp[,2] <- absupp[,2] + abs(temp$upp_mm[(n.burnin + 1):end.iter, par]-true.quantile[2])
    absupp[,3] <- absupp[,3] + abs(temp$upp_eq[(n.burnin + 1):end.iter, par]-true.quantile[2])
    low[,1] <- low[,1] + temp$low_mle[(n.burnin + 1):end.iter, par]
    low[,2] <- low[,2] + temp$low_mm[(n.burnin + 1):end.iter, par]
    low[,3] <- low[,3] + temp$low_eq[(n.burnin + 1):end.iter, par]
    upp[,1] <- upp[,1] + temp$upp_mle[(n.burnin + 1):end.iter, par]
    upp[,2] <- upp[,2] + temp$upp_mm[(n.burnin + 1):end.iter, par]
    upp[,3] <- upp[,3] + temp$upp_eq[(n.burnin + 1):end.iter, par]
  }
  
  abslow <- abslow / nset
  absupp <- absupp / nset
  low <- low / nset
  upp <- upp / nset
  
  df.abs <- data.frame(id = rep((n.burnin + 1):end.iter, n.est), 
                       low = as.vector(abslow),
                       upp = as.vector(absupp),
                       Type = as.factor(rep(1:n.est, each = (end.iter - n.burnin))))
  df.abs$logid <- log(df.abs$id)
  df.abs$loglow <- log(df.abs$low)
  df.abs$logupp <- log(df.abs$upp)
  df.abs$logavg <- log((df.abs$low+df.abs$upp)/2)
  
  plot2 <- ggplot(df.abs, aes(x = logid, group = Type)) +  # Plot for lower bound
    geom_line(aes(y = loglow, col = Type, linetype = Type), lwd=0.75) +
    labs(x = "log(Iterations)", y = "Log absolute difference from true lower bound") +
    ylim(-3.75, 0.2) + 
    scale_linetype_manual(values = rep(1, n.est), labels = c("MLE", "MM", "Empirical Quantile")) + 
    scale_color_manual(values = hex, labels = c("MLE", "MM", "Empirical Quantile")) +
    ggtitle(paste0("Log absolute difference from true lower bound")) +
    theme(legend.position = "none", axis.title.y=element_blank())  
  
  plot3 <- ggplot(df.abs, aes(x = logid, group = Type)) +  # Plot for upper bound
    geom_line(aes(y = logupp, col = Type, linetype = Type), lwd=0.75) +
    labs(x = "log(Iterations)", y = "Log absolute difference from true upper bound") +
    ylim(-3.75, 0.2) + 
    scale_linetype_manual(values = rep(1, n.est), labels = c("MLE", "MM", "Empirical Quantile")) + 
    scale_color_manual(values = hex, labels = c("MLE", "MM", "Empirical Quantile")) +
    ggtitle(paste0("Log absolute difference from true upper bound")) +
    theme(legend.position = "none", 
          axis.text.y=element_blank(), 
          axis.ticks.y=element_blank(),
          axis.title.y=element_blank())  
  
  plot4 <- ggplot(df.abs, aes(x = logid, group = Type)) +  # Plot for average difference
    geom_line(aes(y = logavg, col = Type, linetype = Type), lwd=0.75) +
    labs(x = "log(Iterations)", y = "Log average absolute difference") +
    ylim(-3.75, 0.2) + 
    scale_linetype_manual(values = rep(1, n.est), labels = c("MLE", "MM", "Empirical Quantile")) + 
    scale_color_manual(values = hex, labels = c("MLE", "MM", "Empirical Quantile")) +
    ggtitle(paste0("Log average absolute difference")) +
    theme(axis.text.y=element_blank(), 
          axis.ticks.y=element_blank(),
          axis.title.y=element_blank())
  
  df <- data.frame(id = rep((n.burnin + 1):end.iter, (n.est + 1)), 
                   low = c(as.vector(low),
                           rep(true.quantile[1], (end.iter - n.burnin))),
                   upp = c(as.vector(upp),
                           rep(true.quantile[2], (end.iter - n.burnin))),
                   Type = as.factor(rep(0:n.est, each = (end.iter - n.burnin))))
  
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

low_mle <- upp_mle <- low_mm <- upp_mm <- low_eq <- upp_eq <- matrix(0, ncol = nset * length(phi), nrow = end.iter - n.burnin)

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

plots <- list()
for(i in 1:length(phi)){
  df.var <- data.frame(id = rep((n.burnin + 1):end.iter, n.est), 
                       low = c(low_mle_var[,i], 
                               low_mm_var[,i],
                               low_eq_var[,i]),
                       upp = c(upp_mle_var[,i], 
                               upp_mm_var[,i], 
                               upp_eq_var[,i]),
                      Type = as.factor(rep(1:n.est, each = (end.iter - n.burnin))))
  
  plot1 <- ggplot(df.var, aes(x = id, group = Type)) +
    geom_line(aes(y = low, col = Type, linetype = Type), lwd=0.75) +
    labs(x = "Variance of lower bound", y = "Intervals") +
    ylim(0, 0.7) +
    scale_linetype_manual(values = c(1, 1, 1), labels = c("MLE", "MM", "Empirical Quantile")) + 
    scale_color_manual(values = hex, labels = c("MLE", "MM", "Empirical Quantile")) +
    ggtitle(paste0("Variance of lower bound when phi=", phi[i])) +
    theme(legend.position = "none")
  
  plot2 <- ggplot(df.var, aes(x = id, group = Type)) +
    geom_line(aes(y = upp, col = Type, linetype = Type), lwd=0.75) +
    labs(x = "Variance of upper bound", y = "Intervals") +
    ylim(0, 2.5) +
    scale_linetype_manual(values = c(1, 1, 1), labels = c("MLE", "MM", "Empirical Quantile")) + 
    scale_color_manual(values = hex, labels = c("MLE", "MM", "Empirical Quantile")) +
    ggtitle(paste0("Variance of upper bound when phi=", phi[i])) +
    theme(legend.position = "none")
  
  plotnames <- paste0("plot",1:2,phi[i]*10)
  for(j in 1:2){
    assign(plotnames[j], get(paste0("plot",j)))
  }
  newlist <- mget(plotnames)
  plots <- c(plots, newlist)
}
plot.list <- ggarrange(plots[[1]], plots[[2]], 
                       plots[[3]], plots[[4]],
                       plots[[5]], plots[[6]],
                       nrow = 3, ncol = 2, common.legend = TRUE, legend="bottom")
ggexport(plot.list, filename = "plotvar.pdf", width = 10, height = 10)
