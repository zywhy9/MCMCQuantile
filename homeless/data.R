library(rjags)
library(R2jags)

ModelCode <- "model{
  #Priors
  p.s ~ dunif(0,1) #probability to be seen
  p.maybe ~ dunif(0,1) #probability a plant that says maybe they were seen
  H ~ dnorm(1000, 1/10000^2) T(0,)
  #Model
  M.maybe ~ dbin(p.maybe, M)
  M.yes ~ dbin(p.s, M-M.maybe)
  M.no <- M - M.maybe
  census.unkn ~ dbin(p.s,total)
  total <- round(H) + M.maybe #number of real homeless
}
"

for(set in 1:100){
  
  set.seed(set)
  
  simdata <- readRDS(paste0("sim/sim",set,".rds"))
  
  M.s.maybe.inits <- round(simdata$M.maybe*0.7)
  H.s.inits <- simdata$census.unkn - M.s.maybe.inits
  H.inits <- H.s.inits/0.7
  
  initial.values <- list(list("p.s" = 0.5,
                              "p.maybe" = 0.5,
                              "H" = H.inits),
                         list("p.s" = 0.6,
                              "p.maybe" = 0.3,
                              "H" = H.inits),
                         list("p.s" = 0.4,
                              "p.maybe" = 0.1,
                              "H" = H.inits))
  vars.monitor <- c("p.s", "p.maybe", "H")
  data <- list("M" = simdata$M,
               "M.yes" = simdata$M.yes,
               "M.maybe" = simdata$M.maybe,
               "census.unkn" = simdata$census.unkn)
  
  jagsfit <- jags(data=data, n.chains=3, inits=initial.values,
                  parameters.to.save=vars.monitor, n.iter=30000, n.burnin=0,n.thin=1,
                  DIC=F, model.file=textConnection(ModelCode))
  
  saveRDS(jagsfit$BUGSoutput$sims.array,paste0("res/res",set,".rds"))
}
