# source("../../function/function.R")

full <- readRDS("full.rds")
full_data <- full$BUGSoutput$sims.array

pdf("denplots.pdf")

for(i in c(1:58,60:63)){
  name <- dimnames(full_data)[[3]][i]
  plot(density(full_data[, , i]), main=paste0("Density plot for ", name))
}

graphics.off()

data <- full_data[, , c("Nsuper", "mean.phi")]
data <- cbind(data[, , 1], data[, , 2])
saveRDS(data, "data.rds")
# niter <- 833
# nchain <- 3
# npar <- 1
# 
# mle_res <- mle_analysis(full, niter=niter, nchain=nchain, npar=npar, transf = "logit")
# eq_res <- eq_analysis(full, niter=niter, nchain=nchain, npar=npar)
# mm_res <- mm_analysis(full, niter=niter, nchain=nchain, npar=npar, transf = "logit")
# 
# res <- list(mm_res = mm_res, mle_res = mle_res, eq_res = eq_res)
# saveRDS(res, "res.rds")