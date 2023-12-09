## Plot
true.sd <- sqrt(sd^2 / (1 - phi^2))
true.quantile <-  qnorm(c(0.025, 0.975), mean = mu, sd = true.sd)## True empirical quantile

# n.burnin <- round(1.3 * nburnin)
n.burnin <- 10
# end.iter <- 2500
end.iter <- 10000
par <- 1
n.est <- 3
hex <- hue_pal()(n.est)

df <- data.frame(id = rep((n.burnin + 1):end.iter, (n.est + 1)), 
                 low = c(low_mle[(n.burnin + 1):end.iter, par], 
                         low_mm[(n.burnin + 1):end.iter, par],
                         low_eq[(n.burnin + 1):end.iter, par],
                         rep(true.quantile[1], (end.iter - n.burnin))),
                 upp = c(upp_mle[(n.burnin + 1):end.iter, par], 
                         upp_mm[(n.burnin + 1):end.iter, par], 
                         upp_eq[(n.burnin + 1):end.iter, par],
                         rep(true.quantile[2], (end.iter - n.burnin))),
                 Type = as.factor(rep(0:n.est, each = (end.iter - n.burnin))))


plot19 <- ggplot(df, aes(x = id, group = Type)) + # Full interval
  geom_line(aes(y = low, col = Type, linetype = Type), lwd=0.75) +
  geom_line(aes(y = upp, col = Type, linetype = Type), lwd=0.75) +
  labs(x = "Iterations", y = "Intervals") +
  ylim(-6, 3.5) + 
  scale_linetype_manual(values = c(1, 1, 1, 2), labels = c("MLE", "MM", "Empirical Quantile", "True Interval")) + 
  scale_color_manual(values = c(hex, "black"), labels = c("MLE", "MM", "Empirical Quantile", "True Interval")) +
  ggtitle(paste0("95% intervals when phi=", phi)) +
  theme(legend.position = "none")  


## Check absolute difference with true CI

# n.burnin <- 10
# end.iter <- niter
# n.est <- 3
# hex <- hue_pal()(n.est)

df.abs <- data.frame(id = rep((n.burnin + 1):end.iter, n.est), 
                     low = abs(c(low_mle[(n.burnin + 1):end.iter, par]-true.quantile[1], 
                             low_mm[(n.burnin + 1):end.iter, par]-true.quantile[1],
                             low_eq[(n.burnin + 1):end.iter, par]-true.quantile[1])),
                     upp = abs(c(upp_mle[(n.burnin + 1):end.iter, par]-true.quantile[2], 
                             upp_mm[(n.burnin + 1):end.iter, par]-true.quantile[2],
                             upp_eq[(n.burnin + 1):end.iter, par]-true.quantile[2])),
                     Type = as.factor(rep(1:n.est, each = (end.iter - n.burnin))))
df.abs$logid <- log(df.abs$id)
df.abs$loglow <- log(df.abs$low)
df.abs$logupp <- log(df.abs$upp)
df.abs$logavg <- log((df.abs$low+df.abs$upp)/2)

plot29 <- ggplot(df.abs, aes(x = logid, group = Type)) +  # Plot for lower bound
  geom_line(aes(y = loglow, col = Type, linetype = Type), lwd=0.75) +
  labs(x = "log(Iterations)", y = "Log absolute difference from true lower bound") +
  ylim(-15, 2.5) + 
  scale_linetype_manual(values = c(1, 1, 1), labels = c("MLE", "MM", "Empirical Quantile")) + 
  scale_color_manual(values = hex, labels = c("MLE", "MM", "Empirical Quantile")) +
  ggtitle(paste0("Log absolute difference from true lower bound")) +
  theme(legend.position = "none", axis.title.y=element_blank())  

plot39 <- ggplot(df.abs, aes(x = logid, group = Type)) +  # Plot for upper bound
  geom_line(aes(y = logupp, col = Type, linetype = Type), lwd=0.75) +
  labs(x = "log(Iterations)", y = "Log absolute difference from true upper bound") +
  ylim(-15, 2.5) + 
  scale_linetype_manual(values = c(1, 1, 1), labels = c("MLE", "MM", "Empirical Quantile")) + 
  scale_color_manual(values = hex, labels = c("MLE", "MM", "Empirical Quantile")) +
  ggtitle(paste0("Log absolute difference from true upper bound")) +
  theme(legend.position = "none", 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank())  

plot49 <- ggplot(df.abs, aes(x = logid, group = Type)) +  # Plot for upper bound
  geom_line(aes(y = logavg, col = Type, linetype = Type), lwd=0.75) +
  labs(x = "log(Iterations)", y = "Log average absolute difference") +
  ylim(-15, 2.5) + 
  scale_linetype_manual(values = c(1, 1, 1), labels = c("MLE", "MM", "Empirical Quantile")) + 
  scale_color_manual(values = hex, labels = c("MLE", "MM", "Empirical Quantile")) +
  ggtitle(paste0("Log average absolute difference")) +
  theme(axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank())


plot.list <- ggarrange(plot14, plot24, plot34, plot44,
                       plot16, plot26, plot36, plot46,
                       plot19, plot29, plot39, plot49, nrow = 3, ncol = 4)
ggexport(plot.list, filename = "plot.pdf", width = 30, height = 14)

