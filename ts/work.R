set.seed(1234)
library(maotai)
library(ggplot2)
library(matlib)
library(Matrix)

## Functions

range_mean <- function(x){
  if(x%%2){
    a <- x / 2 + 1
    b <- x
  }else{
    a <- x / 2
    b <- x
  }
  return(c(a,b))
}

## Simulation

n <- 1000
d <- 10

mu <- rep(100, d)
psi <- matrix(0.9, d, d)+diag(0.1, d, d)

mu1 <- 100
mu2 <- rep(100, (d-1))

sigma11 <- psi[1,1]
sigma12 <- psi[1, 2:d]
sigma21 <- psi[2:d, 1]
sigma22 <- psi[2:d, 2:d]

cond <- function(x2){
  mubar <- mu1 + sigma12 %*% solve(sigma22) %*% (x2 - mu2)
  sigmabar <- sigma11 - sigma12 %*% solve(sigma22) %*% sigma21
  return(rnorm(1, mean = mubar, sd = sqrt(sigmabar)))
}

gibbs_sampler <- function(n_samples, d) {
  samples <- matrix(0, nrow = (n_samples+1), ncol = d)
  
  for (i in 2:(n_samples+1)) {
    for(j in 1:d){
      if(j==1){
        a <- samples[(i-1), 2:d]
      }else if(j==d){
        a <- samples[i, 1:(d-1)]
      }else{
        a <- c(samples[(i-1),(j+1):d], samples[i,1:(j-1)])
      }
      samples[i,j] <- cond(a)
    }
  }
  samples <- samples[-1,]
  return(samples)
}

res <- gibbs_sampler(n, d)

## moment estimate

mean_MM <- sd_MM <- low_MM <- upp_MM <- rep(0, n)
for(i in 1:n){
  temp <- range_mean(i)
  tempres <- res[temp[1]:temp[2],1]
  mean_MM[i] <- mean(tempres)
  sd_MM[i] <- sd(tempres)
  low_MM[i] <- quantile(tempres, 0.025)
  upp_MM[i] <- quantile(tempres, 0.975)
}

## MLE

mle_estimate <- function(samples, d, n0){
  X0 <- apply(samples[-n0,], 2, function(x){sum(x)/n0})
  X1 <- apply(samples, 2, function(x){sum(x)/n0})


  temp <- rbind(rep(0,d), samples)
  temp1 <- temp2 <- matrix(0, d, d)
  for(t in 1:n0){
    temp1 <- temp1 + (temp[(t+1),] - X1) %*% t(temp[t,] - X0)
    temp2 <- temp2 + (temp[t,] - X0) %*% t(temp[t,] - X0)
  }

  beta <- temp1 %*% solve(temp2)
  gamma <- X1 - beta %*% X0

  temp3 <- matrix(0, d, d)
  for(t in 1:n0){
    temp3 <- temp3 + (temp[(t+1),] - beta %*% temp[t,] - gamma) %*% t(temp[(t+1),] - beta %*% temp[t,] - gamma)
  }
  delta <- temp3/n0

  mu_mle <- solve(a = (diag(d)-beta), b = gamma)

  A <- solve(beta)
  B <- -t(beta)
  C <- solve(beta) %*% delta
  psi_mle <- sylvester(A, B, C)

  return(list(mu_mle,psi_mle))
  # return(mu_mle)
}

mean_mle <- sd_mle <- low_mle <- upp_mle <- rep(0,n)
for(i in 11:n){
  temp <- mle_estimate(res[1:i,], d, i)
  mean_mle[i] <- temp[[1]][1]
  tempvar <- temp[[2]][1,1]
  sd_mle[i] <- ifelse(tempvar>=0,sqrt(tempvar),0)
  low_mle[i] <- mean_mle[i] - 1.96 * sd_mle[i]
  upp_mle[i] <- mean_mle[i] + 1.96 * sd_mle[i]
}

# mle_estimate <- function(samples){
#   dm <- dim(samples)
#   if(dm[1]==2){
#     Y <- t(as.matrix(samples[2,]))
#     X <- t(as.matrix(samples[1,]))
#   }else{
#     Y <- samples[2:dm[1],]
#     X <- samples[1:(dm[1]-1),]
#   }
#   beta <- solve(t(X) %*% X) %*% t(X) %*% Y
# }


## Plot

#### Full data
df <- data.frame(id = rep(1:n,2), 
                 mean = c(mean_mle, mean_MM), 
                 low = c(low_mle, low_MM), 
                 upp = c(upp_mle, upp_MM), 
                 Type = as.factor(rep(0:1, each = n)))

#### Data with burn-in
n.burnin <- 15
df <- data.frame(id = rep((n.burnin+1):n,2), 
                 mean = c(mean_mle[(n.burnin+1):n], mean_MM[(n.burnin+1):n]), 
                 low = c(low_mle[(n.burnin+1):n], low_MM[(n.burnin+1):n]), 
                 upp = c(upp_mle[(n.burnin+1):n], upp_MM[(n.burnin+1):n]), 
                 Type = as.factor(rep(0:1, each = (n - n.burnin))))


ggplot(df, aes(x = id, group = Type)) +
  geom_line(aes(y = mean, col = Type)) + 
  geom_ribbon(aes(ymin = low, ymax = upp, fill = Type), alpha = 0.3) +
  scale_fill_discrete(name = "Type", labels = c("MLE", "MM")) + 
  scale_color_discrete(name = "Type", labels = c("MLE", "MM")) +
  labs(x = "Iterations", y = "mu")
