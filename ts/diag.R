## Check simulated data
sd <- sqrt(0.5^2 / (1 - 0.9^2) * (1 - phi^2))
data <- ts_sampler(niter, nchain, mu, sd, phi, initial, seed=47)

plot(data[1:10000])
plot(temp$mean_mle[1:10000,])
plot(temp$sd_mle[1:10000,])

## Check specific interval of bounds
lowfoc <- uppfoc <- matrix(0, nrow=200, ncol=nset)
for(set in 1:nset){
  temp <- readRDS(paste0("res995/res",set,".rds")) ## Read Data
  for(j in which(is.na(temp$mean_mle))){
    temp$low_mle[j,1] <- temp$low_mle[(j-1),1]
    temp$upp_mle[j,1] <- temp$upp_mle[(j-1),1]
  }
  lowfoc[,set] <- temp$low_mle[1301:1500, par]
  uppfoc[,set] <- temp$upp_mle[1301:1500, par]
}

plot(apply(lowfoc, 1, mean), type = "l", ylim = c(-6, 13))
lines(apply(lowfoc, 1, max), type = "l", col = "red")
lines(apply(lowfoc, 1, min), type = "l", col = "red")

for(i in 135:196){
  print(which(lowfoc[i,] == max(lowfoc[i,])))
}

## Find the dataset includes maximum value
max.id <- 1
count <- rep(0,2)
for(set in 1:nset){
  # if(set == 85 || set == 74){
  #   next
  # }
  temp <- readRDS(paste0("res995/res",set,".rds")) ## Read Data
  if(sum(is.na(temp$mean_mle))==0){
    next
  }
  temp.na <- which(is.na(temp$mean_mle))
  temp.count <- cbind(rep(set,length(temp.na)), temp.na)
  count <- rbind(count, temp.count)
  max.temp.id <- max(which(is.na(temp$mean_mle)))
  if(max.id<=max.temp.id){
    max.id <- max.temp.id
  }
}
count <- count[-1,]