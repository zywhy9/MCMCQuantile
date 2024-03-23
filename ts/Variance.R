source("../function/function.R")
library(ggplot2)
library(scales)
library(Matrix)
library(rstan)
library(ggpubr)

nsim <- 5000 ## Number of iterations
nset <- 100 ## Number of simulated sets
nburnin <- 80

mmest <- mleest <- eqest <- matrix(0, ncol = nset, nrow = nsim - nburnin)
for(set in 1:nset){
  temp <- readRDS(paste0("res9/res",set,".rds")) ## Read Data
  mmest[,set] <- as.vector(temp$upp_mm[(nburnin + 1):5000,1])
  mleest[,set] <- as.vector(temp$upp_mle[(nburnin + 1):5000,1])
  eqest[,set] <- as.vector(temp$upp_eq[(nburnin + 1):5000,1])
}
mmvar <- apply(mmest, 1, var)
mlevar <- apply(mleest, 1, var)
eqvar <- apply(eqest, 1, var)

hex <- hue_pal()(3)
df <- data.frame(id = rep((nburnin + 1):nsim, 3), 
                 var = c(mmvar, mlevar, eqvar),
                 Type = as.factor(rep(1:3, each = (nsim - nburnin))))

plot <- ggplot(df, aes(x = id, group = Type)) +  # Plot for upper bound
  geom_line(aes(y = log(var), col = Type, linetype = Type), lwd=0.75) +
  labs(x = "Iterations", y = "log(Variance)") +
  scale_linetype_manual(values = rep(1, 3), labels = c("MM", "MLE", "EQ")) + 
  scale_color_manual(values = hex, labels = c("MM", "MLE", "EQ")) +
  ggtitle(paste0("Iterations vs log(Variance) of estimates across 100 datasets")) 
plot
