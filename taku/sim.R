library(simanalyse)
library(coda)
library(MCMCpack)
library(mcmcr)
library(nlist)

set.seed(1234)
load("region1.RData")

## Process data
nweek <- 12
nstock <- 4
N <- 60000
w <- w/sum(w)
regionmean <- region1mean/100
regionsd <- region1sd/100

## Find the constant
t <- rep(NA,12)

f <- function(c){
  var <- mu*(1-mu)/(c+1)
  mean(sum((var-sigma)^2))
}

for(i in 1:12){
  mu <- regionmean[i,]
  sigma <- regionsd[i,]^2
  res <- optim(10,f,method = "Brent",lower=0,upper=1000000)
  t[i] <- res$par
}

## Simulation
simcode <- "
  for(week in 1:nweek){
    X[week,1:nstock] ~ dmulti(p[week,1:nstock],n[week])
    for(stock in 1:nstock){
      mn[week,stock] <- max(1e-10,X[week,stock]/n[week])
      vr[week,stock] <- mn[week,stock]*(1-mn[week,stock])/(t[week]+1)
      sd[week,stock] <- sqrt(vr[week,stock])
      mu[week,stock] <- max(1e-10,mu_v[week,stock])
    }
    mu_v[week,1:nstock] ~ ddirch(t[week]*mn[week,1:nstock])
    mu_l[week] <- mu[week,2] + mu[week,3]
    sd_l[week] <- sqrt(vr[week,2]+vr[week,3])
  }
  Nw <- N * sum(w*p[,2],w*p[,3])
"
working_directory <- "sim/"

sims::sims_simulate(simcode, latent=NA, parameters=nlist(N=N,p=regionmean),
                    save=T,exists=NA,
                    constants = nlist(n=n,t=t,w=w,nweek=nweek,nstock=nstock),
                    nsims=500,monitor=c("mu","Nw","sd","mu_l","sd_l"),
                    path=working_directory)

file.remove("sim/.sims.rds")
current.names <- list.files(working_directory)
new.names <- paste0("sim",1:length(current.names),".rds")
file.rename(paste0(working_directory, current.names), 
            paste0(working_directory, new.names))


