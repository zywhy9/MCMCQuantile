library(rjags)
library(R2jags)
library(MCMCpack)

## Inference
ModelCode <- "
  model{
  for(week in 1:nweek){
    X[week,1] ~ dbin(p[week,1],n[week])
    for(s in 2:(nstock-1)){
      X[week,s] ~ dbin(p[week,s]/(1-sum(p[week,1:(s-1)])),n[week]-sum(X[week,1:(s-1)]))
    }
    X[week,nstock] <- n[week] - sum(X[week,1:(nstock-1)])
    for(stock in 1:nstock){
      mn[week,stock] <- max(1e-10,X[week,stock]/n[week])
    }
    mu[week,1:nstock] ~ ddirch(t[week]*mn[week,1:nstock])
  }
  # Priors
  for(week in 1:nweek){
   p[week,1:nstock] ~ ddirch(rep(1,nstock))
  }
  N <- Nw/sum(w*p[,2],w*p[,3])
  }
"
set.seed(1234)
initial.values <- list(list("p"=rdirichlet(12,rep(1,4))),
                       list("p"=rdirichlet(12,rep(1,4))),
                       list("p"=rdirichlet(12,rep(1,4))))

vars.monitor <- c("N", "p", "X")


# prior <- "
# for(stock in 1:nstock){
#   Z[1,stock] ~ dnorm(0,1/2^2)
#   expZ[1,stock] <- exp(Z[1,stock])
#   for(week in 2:nweek){
#     Z[week,stock] <- phi*Z[(week-1),stock]+epsilon[week,stock]
#     expZ[week,stock] <- exp(Z[week,stock])
#     epsilon[week,stock] ~ dnorm(0,1/(2^2*(1-phi^2)))
#   }
#   for(week in 1:nweek){
#     p[week,stock] <- expZ[week,stock]/sum(expZ[week,])
#   }
# }
# phi ~ dunif(-1,1) 
# N <- Nw/sum(w*p[,2],w*p[,3])"
# 
# initial.values <- list(list("phi"=0.5),
#                        list("phi"=-0.5),
#                        list("phi"=0))
#
# vars.monitor <- c("N", "phi")

for(set in 1:500){
  
  set.seed(set)
  
  simdata <- readRDS(paste0("sim/sim",set,".rds"))
  
  
  data <- list("Nw" = simdata$Nw,
               "mu" = simdata$mu,
               "n" = simdata$n,
               "t" = simdata$t,
               "w" = simdata$w,
               "nweek" = simdata$nweek,
               "nstock" = simdata$nstock)
  
  jagsfit <- jags(data=data, n.chains=3, inits=initial.values,
                  parameters.to.save=vars.monitor, n.iter=10000, n.burnin=0,n.thin=1,
                  DIC=F, model.file=textConnection(ModelCode))
  
  saveRDS(jagsfit$BUGSoutput$sims.array,paste0("res/res",set,".rds"))
}

