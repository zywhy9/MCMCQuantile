library(ggplot2)
library(boot)
library(scales)
library(gridExtra)
library(posterior)
library(MASS)
library(zeallot)
library(grid)
library(ggpubr)
library(latex2exp)

options(scipen = 6)

fulldata <- readRDS("data.rds")
# fulldata1 <- fulldata[, 1:6]
# fulldata2 <- fulldata[, 7:12]
shortdata1 <- fulldata[, 1:3]
shortdata2 <- fulldata[, 7:9]
rm(fulldata)
nburnin <- 0

## True values
true_eq <- t(matrix(c(quantile(shortdata1, c(0.025, 0.975)), quantile(shortdata2, c(0.025, 0.975))), 2, 2))
eq_diff <- apply(true_eq, 1, function(x){x[2]-x[1]})

## Function

boxcox_inv <- function(data, lambda){
  data <- (data * lambda + 1) ^ (1 / lambda)
  return(data)
}

run <- function(data, l, n) {
  len <- min(l, n)
  short <- data[(nburnin + 1):l, ]
  combined <- as.vector(short)
  thinned <- combined[seq(1, length(combined), length.out = len)]
  means <- mean(thinned)
  sds <- sd(thinned)
  lq <- quantile(thinned, 0.025)
  uq <- quantile(thinned, 0.975)
  lq_MM <- means - qnorm(0.975) * sds
  uq_MM <- means + qnorm(0.975) * sds
  qq <- qqnorm(thinned, plot.it = FALSE)
  correlation <- cor(qq[[1]], qq[[2]])
  out <- c(lq, lq_MM, uq, uq_MM, correlation)
  return(out)
}

run_boxcox <- function(data, l, n, lambda) {
  len <- min(l, n)
  short <- data[(nburnin + 1):l, ]
  combined <- as.vector(short)
  thinned <- combined[seq(1, length(combined), length.out = len)]
  means <- mean(thinned)
  sds <- sd(thinned)
  lq_boxcox <- boxcox_inv(means - qnorm(0.975) * sds, lambda)
  uq_boxcox <- boxcox_inv(means + qnorm(0.975) * sds, lambda)
  out <- c(lq_boxcox, uq_boxcox)
  return(out)
}

app_plots <- function(data, n, nburnin, par, lambda) {
  niter <- nrow(data)
  nchain <- ncol(data)
  
  thin_iter <- unique(round(seq(100, niter, length.out = 1000)))
  if(n != niter){
    thin_iter <- unique(pmax(thin_iter, n))
  }
  result <- t(sapply(thin_iter, run, n = n, data = data))
  colnames(result) <- c("lq", "lq_MM", "uq", "uq_MM", "corr")
  
  data_boxcox <- (data ^ lambda - 1) / lambda
  res_boxcox <- t(sapply(thin_iter, run_boxcox, n = n, data = data_boxcox, lambda = lambda))
  colnames(res_boxcox) <- c("lq_boxcox", "uq_boxcox")
  result <- as.data.frame(cbind(result, res_boxcox))
  
  plotlist <- list()
  thinning <- switch(as.character(format(n, scientific=F)), "2500"={"High"}, "10000"={"Medium"}, "100000"={"No"})
  
  ## CI Estimates
  df <- data.frame(
    id = rep(thin_iter, 3),
    low = c(result$lq_MM, result$lq, result$lq_boxcox),
    upp = c(result$uq_MM, result$uq, result$uq_boxcox),
    Type = as.factor(rep(c("MM", "EQ", "MM Box-Cox"), each = length(thin_iter)))
  )
  gg_cil <- ggplot(df, aes(x = id)) +
    geom_line(aes(y = low, col = Type), lwd = 0.5) +
    labs(x = "Iterations", y = "Estimates of the 0.025 Quantile") +
    scale_x_continuous(trans = "log10") +
    scale_color_manual(values = alpha(c("red", "blue", "limegreen"), 0.4), labels = c("EQ", "MM", "MM Box-Cox")) +
    geom_vline(xintercept = 5000) +
    ggtitle(paste0(thinning, " Thinning")) +
    theme_minimal() +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 14),
          plot.title = element_text(hjust = 0.5, size = 14))
  
  gg_ciu <- ggplot(df, aes(x = id)) +
    geom_line(aes(y = upp, col = Type), lwd = 0.5) +
    labs(x = "Iterations", y = "Estimates of the 0.975 Quantile") +
    scale_x_continuous(trans = "log10") +
    scale_color_manual(values = alpha(c("red", "blue", "limegreen"), 0.4), labels = c("EQ", "MM", "MM Box-Cox")) +
    geom_vline(xintercept = 5000) +
    ggtitle(paste0(thinning, " Thinning")) +
    theme_minimal() +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 14),
          plot.title = element_text(hjust = 0.5, size = 14))
  
  ## Squared Error
  ymin = min(c(
    result$lq - true_eq[par, 1],
    result$uq - true_eq[par, 2],
    result$lq_MM - true_eq[par, 1],
    result$uq_MM - true_eq[par, 2],
    result$lq_boxcox - true_eq[par, 1],
    result$uq_boxcox - true_eq[par, 2]
  ))
  ymax = max(c(
    result$lq - true_eq[par, 1],
    result$uq - true_eq[par, 2],
    result$lq_MM - true_eq[par, 1],
    result$uq_MM - true_eq[par, 2],
    result$lq_boxcox - true_eq[par, 1],
    result$uq_boxcox - true_eq[par, 2]
  ))
  ymax <- max(ymin ^ 2, ymax ^ 2)
  df <- data.frame(
    id = thin_iter,
    lq_error = (result$lq - true_eq[par, 1]) ^ 2,
    uq_error = (result$uq - true_eq[par, 2]) ^ 2,
    lq_MM_error = (result$lq_MM - true_eq[par, 1]) ^ 2,
    uq_MM_error = (result$uq_MM - true_eq[par, 2]) ^ 2,
    lq_boxcox_error = (result$lq_boxcox - true_eq[par, 1]) ^ 2,
    uq_boxcox_error = (result$uq_boxcox - true_eq[par, 2]) ^ 2
  )
  gg_error1 <- ggplot(df, aes(x = id)) +
    geom_line(aes(y = lq_error, color = "EQ Lower")) +
    geom_line(aes(y = uq_error, color = "EQ Upper")) +
    geom_line(aes(y = lq_MM_error, color = "MM Lower")) +
    geom_line(aes(y = uq_MM_error, color = "MM Upper")) +
    scale_y_continuous(limits = c(0, ymax)) +
    scale_x_log10(limits = c(100, niter)) +
    scale_color_manual(values = c("MM Lower" = alpha("blue", 0.4), "EQ Lower" = alpha("red", 0.4),
                                  "MM Upper" = alpha("darkblue", 0.4), "EQ Upper" = alpha("darkred", 0.4))) +
    labs(x = "Iterations", y = "Squared Error in Quantile Estimates", color = "Type") +
    geom_vline(xintercept = 5000) +
    ggtitle(paste0(thinning, " Thinning")) +
    theme_minimal() +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 12),
          plot.title = element_text(hjust = 0.5, size = 12))
  
  gg_error2 <- ggplot(df, aes(x = id)) +
    geom_line(aes(y = lq_error, color = "EQ Lower")) +
    geom_line(aes(y = uq_error, color = "EQ Upper")) +
    geom_line(aes(y = lq_boxcox_error, color = "MM Box-Cox Lower")) +
    geom_line(aes(y = uq_boxcox_error, color = "MM Box-Cox Upper")) +
    scale_y_continuous(limits = c(0, ymax)) +
    scale_x_log10(limits = c(100, niter)) +
    scale_color_manual(values = c("MM Box-Cox Lower" = alpha("limegreen", 0.4), "EQ Lower" = alpha("red", 0.4),
                                  "MM Box-Cox Upper" = alpha("darkgreen", 0.4), "EQ Upper" = alpha("darkred", 0.4))) +
    labs(x = "Iterations", y = "Squared Error in Quantile Estimates", color = "Type") +
    geom_vline(xintercept = 5000) +
    ggtitle(paste0(thinning, " Thinning")) +
    theme_minimal() +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 14),
          plot.title = element_text(hjust = 0.5, size = 12))
  
  plotlist <- list(gg_cil, gg_ciu, gg_error1, gg_error2)
  
  return(plotlist)
}

#### Nsuper
b <- boxcox(lm(shortdata1 ~ 1), plotit = F)
lambda <- b$x[which.max(b$y)]
## High thinning (2500), 3 chains
n <- 2500
plot_h <- app_plots(shortdata1, n, nburnin, par = 1, lambda)

## Medium thinning (10000), 3 chains
n <- 10000
plot_m <- app_plots(shortdata1, n, nburnin, par = 1, lambda)

## No thinning (100000), 3 chains
n <- 100000
plot_n <- app_plots(shortdata1, n, nburnin, par = 1, lambda)

save(plot_n, plot_m, plot_h, file = "par1.RData")

#### mean.phi
b <- boxcox(lm(shortdata2 ~ 1), plotit = F)
lambda <- b$x[which.max(b$y)]
## High thinning (2500), 3 chains
n <- 2500
plot_h <- app_plots(shortdata2, n, nburnin, par = 2, lambda)

## Medium thinning (10000), 3 chains
n <- 10000
plot_m <- app_plots(shortdata2, n, nburnin, par = 2, lambda)

## No thinning (100000), 3 chains
n <- 100000
plot_n <- app_plots(shortdata2, n, nburnin, par = 2, lambda)

save(plot_n, plot_m, plot_h, file = "par2.RData")
rm(list = ls(all.names = TRUE))


#### Six types of plots under all the scenarios
plotnames <- paste0(c("cil", "ciu", "error", "errorbc"), rep(1:2, each = 4))
load("par1.RData")
for(j in 1:4){
  temp <- list(plot_n[[j]], plot_m[[j]], plot_h[[j]])
  temp <- lapply(temp, ggplot_add, object = xlab(""))
  temp <- lapply(temp, ggplot_add, object = ylab(""))
  temp <- lapply(temp, ggplot_add, object = theme(legend.position = "none"))
  assign(plotnames[j], temp)
}

legend_grob <- get_legend(plot_n[[1]])
legend_rmse1 <- get_legend(plot_n[[3]])
legend_rmse2 <- get_legend(plot_n[[4]])


load("par2.RData")
for(j in 1:4){
  temp <- list(plot_n[[j]], plot_m[[j]], plot_h[[j]])
  temp <- lapply(temp, ggplot_add, object = xlab(""))
  temp <- lapply(temp, ggplot_add, object = ylab(""))
  temp <- lapply(temp, ggplot_add, object = theme(legend.position = "none"))
  assign(plotnames[4+j], temp)
}
rm(temp, plot_n, plot_m, plot_h)

## Tune xlimit ylimit
cil1 <- lapply(cil1, ggplot_add, object = scale_x_log10(limits = c(2500, 100000)))
cil1 <- lapply(cil1, ggplot_add, object = scale_y_continuous(limits = c(192, 195.5)))
ciu1 <- lapply(ciu1, ggplot_add, object = scale_x_log10(limits = c(2500, 100000)))
ciu1 <- lapply(ciu1, ggplot_add, object = scale_y_continuous(limits = c(217.5, 221.5)))
error1 <- lapply(error1, ggplot_add, object = scale_x_log10(limits = c(2500, 100000)))
error1 <- lapply(error1, ggplot_add, object = scale_y_continuous(limits = c(0, 4.5)))
errorbc1 <- lapply(errorbc1, ggplot_add, object = scale_x_log10(limits = c(2500, 100000)))
errorbc1 <- lapply(errorbc1, ggplot_add, object = scale_y_continuous(limits = c(0, 4.5)))

cil2 <- lapply(cil2, ggplot_add, object = scale_x_log10(limits = c(2500, 100000)))
cil2 <- lapply(cil2, ggplot_add, object = scale_y_continuous(limits = c(0.655, 0.675)))
ciu2 <- lapply(ciu2, ggplot_add, object = scale_x_log10(limits = c(2500, 100000)))
ciu2 <- lapply(ciu2, ggplot_add, object = scale_y_continuous(limits = c(0.807, 0.825)))
error2 <- lapply(error2, ggplot_add, object = scale_x_log10(limits = c(2500, 100000)))
error2 <- lapply(error2, ggplot_add, object = scale_y_continuous(limits = c(0, 0.00007)))
errorbc2 <- lapply(errorbc2, ggplot_add, object = scale_x_log10(limits = c(2500, 100000)))
errorbc2 <- lapply(errorbc2, ggplot_add, object = scale_y_continuous(limits = c(0, 0.00007)))

est1 <- c(cil1, ciu1)[c(1,4,2,5,3,6)]
est2 <- c(cil2, ciu2)[c(1,4,2,5,3,6)]
err1 <- c(error1, errorbc1)[c(1,4,2,5,3,6)]
err2 <- c(error2, errorbc2)[c(1,4,2,5,3,6)]


estplot1 <- grid.arrange(
  arrangeGrob(grobs = est1, ncol = 2, nrow = 3),
  legend_grob,
  ncol = 1,
  heights = unit(c(30, 1), c("null", "null"))
)

ggexport(
  as_ggplot(estplot1),
  filename = "pub/estplot1.pdf",
  width = 15,
  height = 14
)

estplot2 <- grid.arrange(
  arrangeGrob(grobs = est2, ncol = 2, nrow = 3),
  legend_grob,
  ncol = 1,
  heights = unit(c(30, 1), c("null", "null"))
)

ggexport(
  as_ggplot(estplot2),
  filename = "pub/estplot2.pdf",
  width = 15,
  height = 14
)


rmseplot1 <- grid.arrange(
  arrangeGrob(grobs = list(err1[[1]], err1[[3]], err1[[5]]), ncol = 1, nrow = 3),
  arrangeGrob(grobs = list(err1[[2]], err1[[4]], err1[[6]]), ncol = 1, nrow = 3),
  legend_rmse1, legend_rmse2,
  ncol = 2,
  heights = unit(c(30, 1), c("null", "null"))
)

ggexport(
  as_ggplot(rmseplot1),
  filename = "pub/rmseplot1.pdf",
  width = 15,
  height = 14
)


rmseplot2 <- grid.arrange(
  arrangeGrob(grobs = list(err2[[1]], err2[[3]], err2[[5]]), ncol = 1, nrow = 3),
  arrangeGrob(grobs = list(err2[[2]], err2[[4]], err2[[6]]), ncol = 1, nrow = 3),
  legend_rmse1, legend_rmse2,
  ncol = 2,
  heights = unit(c(30, 1), c("null", "null"))
)

ggexport(
  as_ggplot(rmseplot2),
  filename = "pub/rmseplot2.pdf",
  width = 15,
  height = 14
)

#### Posterior Density Plot

fulldata <- readRDS("data.rds")
shortdata1 <- fulldata[, 1:3]
shortdata2 <- fulldata[, 7:9]
rm(fulldata)
b1 <- boxcox(lm(shortdata1 ~ 1), plotit = F)
lambda1 <- b1$x[which.max(b1$y)]
data_boxcox1 <- (shortdata1 ^ lambda1 - 1) / lambda1
b2 <- boxcox(lm(shortdata2 ~ 1), plotit = F)
lambda2 <- b2$x[which.max(b2$y)]
data_boxcox2 <- (shortdata2 ^ lambda2 - 1) / lambda2


final_samples <- function(data, thin=F, l=5000){
  short <- data[1:l, ]
  combined <- as.vector(short)
  if(thin){
    final_data <- combined[seq(1, length(combined), length.out = 2500)]
  }else{
    final_data <- combined
  }
  return(final_data)
}

## Ns
postN1 <- final_samples(shortdata1, thin = T, l = 5000) # high thin 5000
postN2 <- final_samples(shortdata1, thin = F, l = 100000) # no thin 100000
postN3 <- final_samples(data_boxcox1, thin = T, l = 5000) # high thin Box-Cox 5000
postN4 <- final_samples(data_boxcox1, thin = F, l = 100000) # no thin Box-Cox 100000

combined_data <- data.frame(
  value = c(postN1, postN2),
  dataset = factor(rep(c("High Thinning", "No Thinning"),
                       times = c(length(postN1), length(postN2))))
)

combined_data_boxcox <- data.frame(
  value = c(postN3, postN4),
  dataset = factor(rep(c("High Thinning", "No Thinning"),
                       times = c(length(postN3), length(postN4))))
)

# Plot the density lines using ggplot2
postplot11 <- ggplot(combined_data, aes(x = value, color = dataset)) +
  geom_density() +
  labs(title = TeX(r'(Density Plot for Untransformed Samples of $N_{s}$)'),
       x = "Value",
       y = "Density",
       color = "Type") +
  scale_color_manual(values = alpha(c("blue", "red"), 0.4), labels = c("No Thinning", "High Thinning")) +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 14),
        plot.title = element_text(hjust = 0.5, size = 14))

legend_grob <- get_legend(postplot11)

postplot11 <- ggplot(combined_data, aes(x = value, color = dataset)) +
  geom_density() +
  labs(title = TeX(r'(Density Plot for Untransformed Samples of $N_{s}$)'),
       x = "",
       y = "",
       color = "Type") +
  scale_color_manual(values = alpha(c("blue", "red"), 0.4), labels = c("No Thinning", "High Thinning")) +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 14))

postplot12 <- ggplot(combined_data_boxcox, aes(x = value, color = dataset)) +
  geom_density() +
  labs(title = TeX(r'(Density Plot for Box-Cox Transformed Samples of $N_{s}$)'),
       x = "",
       y = "",
       color = "Type") +
  scale_color_manual(values = alpha(c("blue", "red"), 0.4), labels = c("No Thinning", "High Thinning")) +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 14))

## mean.phi
postphi1 <- final_samples(shortdata2, thin = T, l = 5000) # high thin
postphi2 <- final_samples(shortdata2, thin = F, l = 100000) # no thin 
postphi3 <- final_samples(data_boxcox2, thin = T, l = 5000) # high thin Box-Cox
postphi4 <- final_samples(data_boxcox2, thin = F, l = 100000) # no thin Box-Cox

combined_data <- data.frame(
  value = c(postphi1, postphi2),
  dataset = factor(rep(c("High Thinning", "No Thinning"),
                       times = c(length(postN1), length(postN2))))
)

combined_data_boxcox <- data.frame(
  value = c(postphi3, postphi4),
  dataset = factor(rep(c("High Thinning", "No Thinning"),
                       times = c(length(postN3), length(postN4))))
)

# Plot the density lines using ggplot2
postplot21 <- ggplot(combined_data, aes(x = value, color = dataset)) +
  geom_density() +
  labs(title = TeX(r'(Density Plot for Untransformed Samples of $\bar{\phi}$)'),
       x = "",
       y = "",
       color = "Type") +
  scale_color_manual(values = alpha(c("blue", "red"), 0.4), labels = c("No Thinning", "High Thinning")) +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 14))

postplot22 <- ggplot(combined_data_boxcox, aes(x = value, color = dataset)) +
  geom_density() +
  labs(title = TeX(r'(Density Plot for Box-Cox Transformed Samples of $\bar{\phi}$)'),
       x = "",
       y = "",
       color = "Type") +
  scale_color_manual(values = alpha(c("blue", "red"), 0.4), labels = c("No Thinning", "High Thinning")) +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 14))


postplot <- grid.arrange(
  arrangeGrob(
    postplot11,
    postplot21,
    postplot12,
    postplot22,
    ncol = 2,
    nrow = 2
  ),
  legend_grob,
  ncol = 1,
  heights = unit(c(30, 1), c("null", "null"))
)

ggexport(
  as_ggplot(postplot),
  filename = "pub/postplot.pdf",
  width = 15,
  height = 14
)
