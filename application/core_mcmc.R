set.seed(1234)

library(R2jags)

# Specify model in BUGS language
sink("js-tempran.jags")
cat("
model {
# Priors and constraints
for (i in 1:M){
   for (t in 1:(n.occasions-1)){
      logit(phi[i,t]) <- mean.lphi + epsilon[t]
      } #t
   for (t in 1:n.occasions){
      p[i,t] <- mean.p
      } #t
   } #i
mean.p ~ dunif(0, 1)                # Prior for mean capture
mean.phi ~ dunif(0, 1)              # Prior for mean survival
mean.lphi <- log(mean.phi / (1-mean.phi))
for (t in 1:(n.occasions-1)){
   epsilon[t] ~ dnorm(0, tau)
   }
tau <- pow(sigma, -2)
sigma ~ dunif(0, 5)              # Prior for sd of indv. variation of phi
sigma2 <- pow(sigma, 2)

for (t in 1:n.occasions){
   gamma[t] ~ dunif(0, 1)
   } #t

# Likelihood
for (i in 1:M){
   # First occasion
   # State process
   z[i,1] ~ dbern(gamma[1])
   mu1[i] <- z[i,1] * p[i,1]
   # Observation process
   y[i,1] ~ dbern(mu1[i])
   
   # Subsequent occasions
   for (t in 2:n.occasions){
      # State process
      q[i,t-1] <- 1-z[i,t-1]
      mu2[i,t] <- phi[i,t-1] * z[i,t-1] + gamma[t] * prod(q[i,1:(t-1)]) 
      z[i,t] ~ dbern(mu2[i,t])
      # Observation process
      mu3[i,t] <- z[i,t] * p[i,t]
      y[i,t] ~ dbern(mu3[i,t])
      } #t
   } #i

# Calculate derived population parameters
for (t in 1:n.occasions){
   qgamma[t] <- 1-gamma[t]
   }
cprob[1] <- gamma[1]
for (t in 2:n.occasions){
   cprob[t] <- gamma[t] * prod(qgamma[1:(t-1)])
   } #t
psi <- sum(cprob[])            # Inclusion probability
for (t in 1:n.occasions){
   b[t] <- cprob[t] / psi      # Entry probability
   } #t
for (i in 1:M){
   recruit[i,1] <- z[i,1]
   for (t in 2:n.occasions){
      recruit[i,t] <- (1-z[i,t-1]) * z[i,t]
      } #t
   } #i
for (t in 1:n.occasions){
   N[t] <- sum(z[1:M,t])        # Actual population size
   B[t] <- sum(recruit[1:M,t])  # Number of entries
   } #t
for (i in 1:M){
   Nind[i] <- sum(z[i,1:n.occasions])
   Nalive[i] <- 1-equals(Nind[i], 0)
   } #i
Nsuper <- sum(Nalive[])         # Size of superpopulation
}
",fill=TRUE)
sink()


leis <- as.matrix(read.table("leisleri.txt", sep = " ", header = FALSE))
nz <- 300
CH.aug <- rbind(leis, matrix(0, ncol = dim(leis)[2], nrow = nz))

# Bundle data
jags.data <- list(y = CH.aug,
                  n.occasions = dim(CH.aug)[2],
                  M = dim(CH.aug)[1])

# Initial values
# Good initial values for the latent state z are needed to run the model in JAGS. The simplest option that works is to give just a matrix with a 1 at all places.
z.init <- CH.aug
z.init[z.init == 0] <- 1
inits <- function() {
  list(
    mean.phi = runif(1, 0, 1),
    mean.p = runif(1, 0, 1),
    sigma = runif(1, 0, 1),
    z = z.init
  )
}

# Parameters monitored
parameters <- c("psi", "mean.p", "sigma2", "mean.phi", "N", "Nsuper", "b", "B")

# MCMC settings
ni <- 10000
nt <- 6
nb <- 5000
nc <- 3

# Call JAGS from R (BRT 127 min)
t1 <- Sys.time()
nl <- jags(
  jags.data,
  inits,
  parameters,
  "js-tempran.jags",
  n.chains = nc,
  n.thin = nt,
  n.iter = ni,
  n.burnin = nb,
  working.directory = getwd()
)
t2 <- Sys.time()

nl$tbar = difftime(t2, t1, units = "mins")

saveRDS(nl, file = "full.rds")
  
