

##############################################

#This script analyzes the real data of PM 2.5 
#monitoring stations. 

###############################################



library(shape)
library(MBA)
library(RColorBrewer)
library(fields)
library(MASS)
library(geoR)
library(rjags)
#library(spBayes)
library(dplyr)
library(ggplot2)
library(magrittr)
library(MASS)
library(mvtnorm)
library(invgamma)
library(truncnorm)



library(rgdal) 
pm <- read.csv("PM2.5_EasternUS.csv", header=T)
long2UTM <- function(long) {
  (floor((long + 180)/6) %% 60) + 1
}
pm$zone <- sapply(pm$long,long2UTM)

LongLatToUTM <-function(long,lat,zone){
  xy <- data.frame(Long = long, Lat = lat)
  coordinates(xy) <- c("Long", "Lat")
  proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")  ## for example
  res <- t(as.data.frame(spTransform(xy, CRS(paste("+proj=utm +zone=",zone," ellps=WGS84",sep='')))))
  return(res)
}

utm <- as.data.frame(t(mapply(LongLatToUTM, pm$longitude, pm$latitude, pm$zone)))
colnames(utm) <- c("Easting","Northing")
utm <- cbind(utm, Long=pm$longitude, Lat=pm$latitude)


# total of 546 sites, 20% test, 109 test site
numTotal <- nrow(pm)
set.seed(1127)
ind_test <- sample(1:numTotal, size=round(numTotal*0.2,0), replace=FALSE)
ind_train <- setdiff(1:numTotal, ind_test)
numSites <- length(ind_train)
sites <- as.matrix(cbind(utm$Easting/1000, utm$Northing/1000))
sites_test <- sites[ind_test,]
sites_train <- sites[ind_train,]
total <- as.data.frame(rbind(sites_train, sites_test))
total$category <- c(rep("Fitted Sites",nrow(sites_train)), rep("Hold Out Sites",nrow(sites_test)))
colnames(total)[1:2] <- c("Easting", "Northing")
y_train <- pm$log_PM_25[ind_train]
y_test <- pm$log_PM_25[ind_test]

max(rdist(sites))



coords <- as.matrix(cbind(utm$Long, utm$Lat))
coords_test <- coords[ind_test,]
coords_train <- coords[ind_train,]
df <- as.data.frame(rbind(coords_train, coords_test))
df$category <- c(rep("Fitted Sites",nrow(coords_train)), rep("Hold Out Sites",nrow(coords_test)))
colnames(df)[1:2] <- c("Longitude", "Latitude")

pdf("pm_lat.pdf", height = 5, width=5)
ggplot(df, aes(Longitude, Latitude, shape = category, color = category)) +
  scale_color_manual(values=c("black", "red")) +
  theme(text = element_text(size = 15)) +
  geom_point(size=2.5) + 
  xlab("Longitude") + 
  ylab("Latitude") +
  theme(legend.position="none")
dev.off()


pdf("pm_utm.pdf", height = 5)
ggplot(total, aes(Easting, Northing, shape = category, color = category)) +
  scale_color_manual(values=c("black", "red")) +
  theme(text = element_text(size = 15)) +
  geom_point(size=2.5) + 
  xlab("Easting (km)") + 
  ylab("Northing (km)")
dev.off()



par(mfrow=c(1,1))
pdf("surface.pdf")
surface <- mba.surf(cbind(sites_train, y_train),no.X = 100, no.Y = 100, h=5, extend = FALSE)$xyz.est
image.plot(surface,xaxs="r",yaxs="r",col = gray((0:32)/32))
contour(surface,add = T)
dev.off()




## distance matrix for x and y coordinate
Dx <- matrix(0, nrow=numSites, ncol=numSites)
Dy <- matrix(0, nrow=numSites, ncol=numSites)
for (m in 1:(numSites-1)) {
  for (n in (m+1):numSites) {
    h <- sites_train[m,] - sites_train[n,]
    Dx[m,n] <- h[1]
    Dx[n,m] <- h[1]
    Dy[m,n] <- h[2]
    Dy[n,m] <- h[2]
  }
}

makeB <- function(B){  
  b1 <- B[1,1]
  b2 <- B[2,1]
  b3 <- B[1,2]
  b4 <- B[2,2]
  pwr <- Dx^2 * b1 + Dx * Dy * (b2 + b3) + Dy^2 * b4
  return(pwr)
}






# Metropolis

# initialize
numSim <- 30000
mcmc.sigma <- rep(0.2, numSim)
mcmc.a <- rep(1.5, numSim)
mcmc.r <- rep(5, numSim)
mcmc.tau <- rep(0.18, numSim)
mcmc.mu <- rep(1.8, numSim)
mcmc.phi <- rep(0.003, numSim)
mcmc.kappa <- rep(1.5, numSim)


phi_prior0 <- 3/(0.5*max(rdist(sites)))
phi_prior1 <- 3/(0.2*max(rdist(sites)))

# Metropolis

for(i in 2:numSim){
  
  # update kappa
  var <- 2
  kappa_star <- rlnorm(1, log(mcmc.kappa[i-1]), sqrt(var))
  amin <- 1
  amax <- amin*mcmc.r[i-1]
  rotationMat <- matrix(c(cos(mcmc.a[i-1]),-sin(mcmc.a[i-1]),sin(mcmc.a[i-1]),cos(mcmc.a[i-1])),nrow=2,ncol=2)
  aMat <- matrix(c(1/amax,0,0,1/amin),nrow=2,ncol=2)
  A <- rotationMat %*% aMat
  B <- A %*% t(A)
  # covariances
  Sigma_star <- cov.spatial(sqrt(makeB(B)), cov.model= "matern", cov.pars=c(mcmc.sigma[i-1], 1/mcmc.phi[i-1]), kappa = kappa_star)
  Sigma_prev <- cov.spatial(sqrt(makeB(B)), cov.model= "matern", cov.pars=c(mcmc.sigma[i-1], 1/mcmc.phi[i-1]), kappa = mcmc.kappa[i-1])
  diag(Sigma_star) <- mcmc.sigma[i-1] + mcmc.tau[i-1]
  diag(Sigma_prev) <- mcmc.sigma[i-1] + mcmc.tau[i-1]
  # calulate ratio
  mvn_star <- dmvnorm(y_train, rep(mcmc.mu[i-1], numSites), Sigma_star, log=TRUE)
  mvn_prev <- dmvnorm(y_train, rep(mcmc.mu[i-1], numSites), Sigma_prev, log=TRUE)
  prior_star <- dunif(kappa_star, 0.1, 2)
  prior_prev <- dunif(mcmc.kappa[i-1], 0.1, 2)
  ratio <- exp(mvn_star-mvn_prev) * prior_star / prior_prev
  u <- runif(1)
  if(log(u) < log(ratio)){
    mcmc.kappa[i] <- kappa_star
  } else {
    mcmc.kappa[i] <- mcmc.kappa[i-1]
  }
  
  
  
  
  # update phi 
  var <- 2
  phi_star <- rlnorm(1, log(mcmc.phi[i-1]), sqrt(var))
  amin <- 1
  amax <- amin*mcmc.r[i-1]
  rotationMat <- matrix(c(cos(mcmc.a[i-1]),-sin(mcmc.a[i-1]),sin(mcmc.a[i-1]),cos(mcmc.a[i-1])),nrow=2,ncol=2)
  aMat <- matrix(c(1/amax,0,0,1/amin),nrow=2,ncol=2)
  A <- rotationMat %*% aMat
  B <- A %*% t(A)
  # covariances=
  Sigma_star <- cov.spatial(sqrt(makeB(B)), cov.model= "matern", cov.pars=c(mcmc.sigma[i-1], 1/phi_star), kappa = mcmc.kappa[i])
  Sigma_prev <- cov.spatial(sqrt(makeB(B)), cov.model= "matern", cov.pars=c(mcmc.sigma[i-1], 1/mcmc.phi[i-1]), kappa = mcmc.kappa[i])
  diag(Sigma_star) <- mcmc.sigma[i-1] + mcmc.tau[i-1]
  diag(Sigma_prev) <- mcmc.sigma[i-1] + mcmc.tau[i-1]
  # calulate ratio
  mvn_star <- dmvnorm(y_train, rep(mcmc.mu[i-1], numSites), Sigma_star, log=TRUE)
  mvn_prev <- dmvnorm(y_train, rep(mcmc.mu[i-1], numSites), Sigma_prev, log=TRUE)
  prior_star <- dunif(phi_star, phi_prior0, phi_prior1)
  prior_prev <- dunif(mcmc.phi[i-1], phi_prior0, phi_prior1)
  ratio <- exp(mvn_star-mvn_prev) * prior_star / prior_prev
  u <- runif(1)
  if(log(u) < log(ratio)){
    mcmc.phi[i] <- phi_star
  } else {
    mcmc.phi[i] <- mcmc.phi[i-1]
  }
  
  
  
  # update sigma squared
  
  var <- 3
  sigma_star <- rlnorm(1, log(mcmc.sigma[i-1]), sqrt(var))
  amin <- 1
  amax <- amin*mcmc.r[i-1]
  rotationMat <- matrix(c(cos(mcmc.a[i-1]),-sin(mcmc.a[i-1]),sin(mcmc.a[i-1]),cos(mcmc.a[i-1])),nrow=2,ncol=2)
  aMat <- matrix(c(1/amax,0,0,1/amin),nrow=2,ncol=2)
  A <- rotationMat %*% aMat
  B <- A %*% t(A)
  # covariances
  Sigma_star <- cov.spatial(sqrt(makeB(B)), cov.model= "matern", cov.pars=c(sigma_star, 1/mcmc.phi[i]), kappa = mcmc.kappa[i])
  Sigma_prev <- cov.spatial(sqrt(makeB(B)), cov.model= "matern", cov.pars=c(mcmc.sigma[i-1], 1/mcmc.phi[i]), kappa = mcmc.kappa[i])
  diag(Sigma_star) <- sigma_star + mcmc.tau[i-1]
  diag(Sigma_prev) <- mcmc.sigma[i-1] + mcmc.tau[i-1]
  mvn_star <- dmvnorm(y_train, rep(mcmc.mu[i-1], numSites), Sigma_star, log=TRUE)
  mvn_prev <- dmvnorm(y_train, rep(mcmc.mu[i-1], numSites), Sigma_prev, log=TRUE)
  prior_star <- dinvgamma(sigma_star, shape=1, scale=1)
  prior_prev <- dinvgamma(mcmc.sigma[i-1], shape=1, scale=1)
  ratio <- exp(mvn_star-mvn_prev) * prior_star / prior_prev
  u <- runif(1)
  if(log(u) < log(ratio)){
    mcmc.sigma[i] <- sigma_star
  } else {
    mcmc.sigma[i] <- mcmc.sigma[i-1]
  }
  
  
  
  # update tau sqaured 
  var <- 2
  tau_star <- rlnorm(1, log(mcmc.tau[i-1]), sqrt(var))
  amin <- 1
  amax <- amin*mcmc.r[i-1]
  rotationMat <- matrix(c(cos(mcmc.a[i-1]),-sin(mcmc.a[i-1]),sin(mcmc.a[i-1]),cos(mcmc.a[i-1])),nrow=2,ncol=2)
  aMat <- matrix(c(1/amax,0,0,1/amin),nrow=2,ncol=2)
  A <- rotationMat %*% aMat
  B <- A %*% t(A)
  # covariance sampled
  Sigma <- cov.spatial(sqrt(makeB(B)), cov.model= "matern", cov.pars=c(mcmc.sigma[i], 1/mcmc.phi[i]), kappa = mcmc.kappa[i])
  Sigma_prev <- Sigma
  Sigma_star <- Sigma
  diag(Sigma_star) <- mcmc.sigma[i] + tau_star
  diag(Sigma_prev) <- mcmc.sigma[i] + mcmc.tau[i-1]
  # calulate ratio
  mvn_star <- dmvnorm(y_train, rep(mcmc.mu[i-1], numSites), Sigma_star, log=TRUE)
  mvn_prev <- dmvnorm(y_train, rep(mcmc.mu[i-1], numSites), Sigma_prev, log=TRUE)
  prior_star <- dinvgamma(tau_star, shape=1, scale=1)
  prior_prev <- dinvgamma(mcmc.tau[i-1], shape=1, scale=1)
  ratio <- exp(mvn_star-mvn_prev) * prior_star / prior_prev
  u <- runif(1)
  if(log(u) < log(ratio)){
    mcmc.tau[i] <- tau_star
  } else {
    mcmc.tau[i] <- mcmc.tau[i-1]
  }
  
  
  
  
  # update angle and ratio
  var1 <- 10
  tt <- rnorm(1, mcmc.a[i-1], sqrt(var1))
  a_star <- tt%%(pi)
  var2 <- 8
  r_star <- rtruncnorm(1, a=0, b=Inf, mean = mcmc.r[i-1], sd = sqrt(var2))
  amin <- 1
  # angle ratio sampled
  amax <- amin*r_star
  rotationMat <- matrix(c(cos(a_star),-sin(a_star),sin(a_star),cos(a_star)),nrow=2,ncol=2)
  aMat <- matrix(c(1/amax,0,0,1/amin),nrow=2,ncol=2)
  A <- rotationMat %*% aMat
  B_star <- A %*% t(A)
  # angle previous
  amax <- amin*mcmc.r[i-1]
  rotationMat <- matrix(c(cos(mcmc.a[i-1]),-sin(mcmc.a[i-1]),sin(mcmc.a[i-1]),cos(mcmc.a[i-1])),nrow=2,ncol=2)
  A <- rotationMat %*% aMat
  B_prev <- A %*% t(A)
  #covariances
  Sigma_star <- cov.spatial(sqrt(makeB(B_star)), cov.model= "matern", cov.pars=c(mcmc.sigma[i], 1/mcmc.phi[i]), kappa = mcmc.kappa[i])
  Sigma_prev <- cov.spatial(sqrt(makeB(B_prev)), cov.model= "matern", cov.pars=c(mcmc.sigma[i], 1/mcmc.phi[i]), kappa = mcmc.kappa[i])
  diag(Sigma_star) <- mcmc.sigma[i] + mcmc.tau[i]
  diag(Sigma_prev) <- mcmc.sigma[i] + mcmc.tau[i]
  # calulate ratio
  mvn_star <- dmvnorm(y_train, rep(mcmc.mu[i-1], numSites), Sigma_star, log=TRUE)
  mvn_prev <- dmvnorm(y_train, rep(mcmc.mu[i-1], numSites), Sigma_prev, log=TRUE)
  prior_star <- dunif(a_star, 0, pi) * dinvgamma(r_star, shape=1, scale=1)
  prior_prev <- dunif(mcmc.a[i-1], 0, pi) * dinvgamma(mcmc.r[i-1], shape=1, scale=1)
  ratio <- exp(mvn_star-mvn_prev) * prior_star / prior_prev
  u <- runif(1)
  if(log(u) < log(ratio)){
    mcmc.a[i] <- a_star
    mcmc.r[i] <- r_star
  } else {
    mcmc.a[i] <- mcmc.a[i-1]
    mcmc.r[i] <- mcmc.r[i-1]
  }
  
  
  
  
  # update mu
  var <- 3
  mu_star <- rnorm(1, mcmc.mu[i-1], sqrt(var))
  amin <- 1
  amax <- amin*mcmc.r[i]
  rotationMat <- matrix(c(cos(mcmc.a[i]),-sin(mcmc.a[i]),sin(mcmc.a[i]),cos(mcmc.a[i])),nrow=2,ncol=2)
  aMat <- matrix(c(1/amax,0,0,1/amin),nrow=2,ncol=2)
  A <- rotationMat %*% aMat
  B <- A %*% t(A)
  # covariance sampled
  Sigma <- cov.spatial(sqrt(makeB(B)), cov.model= "matern", cov.pars=c(mcmc.sigma[i], 1/mcmc.phi[i]), kappa = mcmc.kappa[i])
  diag(Sigma) <- mcmc.sigma[i] + mcmc.tau[i]
  # calulate ratio
  mvn_star <- dmvnorm(y_train, rep(mu_star, numSites), Sigma, log=TRUE)
  mvn_prev <- dmvnorm(y_train, rep(mcmc.mu[i-1], numSites), Sigma, log=TRUE)
  prior_star <- dnorm(mu_star, 0, 100)
  prior_prev <- dnorm(mcmc.mu[i-1], 0, 100)
  ratio <- exp(mvn_star-mvn_prev) * prior_star / prior_prev
  u <- runif(1)
  if(log(u) < log(ratio)){
    mcmc.mu[i] <- mu_star
  } else {
    mcmc.mu[i] <- mcmc.mu[i-1]
  }
  
  
}

save(mcmc.sigma, mcmc.tau, mcmc.a, mcmc.r, mcmc.phi, mcmc.kappa, mcmc.mu, file="mcmc_pm_matern.Rdata")




# trace plot
pdf("trace.pdf")
par(mfrow=c(3,3))
plot(1:numSim, mcmc.a,pch = 20, type="l", xlab="iteration",ylab=expression(alpha))
plot(1:numSim, mcmc.r,pch = 20, type="l", xlab="iteration",ylab="r")
plot(1:numSim, mcmc.sigma,pch = 20, type="l", xlab="iteration",ylab=expression(sigma^2))
plot(1:numSim, mcmc.tau,pch = 20, type="l", xlab="iteration",ylab=expression(tau^2))
plot(1:numSim, mcmc.phi,pch = 20, type="l", xlab="iteration",ylab=expression(phi))
plot(1:numSim, mcmc.mu,pch = 20, type="l", xlab="iteration",ylab=expression(mu))
plot(1:numSim, mcmc.kappa,pch = 20, type="l", xlab="iteration",ylab=expression(nu))
dev.off()

pdf("runave.pdf")
par(mfrow=c(3,3))
runAve <- NULL
runAve[1] <- mcmc.a[1]
for(i in 2:length(mcmc.a)){
 runAve[i] <- (runAve[i-1] * (i-1) + mcmc.a[i]) / i
}
plot(runAve, type="l", xlab = "iteration", ylab = expression(alpha))
runAve <- NULL
runAve[1] <- mcmc.r[1]
for(i in 2:length(mcmc.r)){
 runAve[i] <- (runAve[i-1] * (i-1) + mcmc.r[i]) / i
}
plot(runAve, type="l", xlab = "iteration", ylab="r")
runAve <- NULL
runAve[1] <- mcmc.sigma[1]
for(i in 2:length(mcmc.sigma)){
 runAve[i] <- (runAve[i-1] * (i-1) + mcmc.sigma[i]) / i
}
plot(runAve, type="l", xlab = "iteration", ylab=expression(sigma^2))
runAve <- NULL
runAve[1] <- mcmc.tau[1]
for(i in 2:length(mcmc.tau)){
 runAve[i] <- (runAve[i-1] * (i-1) + mcmc.tau[i]) / i
}
plot(runAve, type="l", xlab = "iteration", ylab=expression(tau^2))
runAve <- NULL
runAve[1] <- mcmc.phi[1]
for(i in 2:length(mcmc.phi)){
  runAve[i] <- (runAve[i-1] * (i-1) + mcmc.phi[i]) / i
}
plot(runAve, type="l", xlab = "iteration", ylab=expression(phi))
runAve <- NULL
runAve[1] <- mcmc.mu[1]
for(i in 2:length(mcmc.mu)){
  runAve[i] <- (runAve[i-1] * (i-1) + mcmc.mu[i]) / i
}
plot(runAve, type="l", xlab = "iteration",ylab=expression(mu))
runAve <- NULL
runAve[1] <- mcmc.kappa[1]
for(i in 2:length(mcmc.kappa)){
 runAve[i] <- (runAve[i-1] * (i-1) + mcmc.kappa[i]) / i
}
plot(runAve, type="l", xlab = "iteration", ylab=expression(nu))

dev.off()





numBurn <- 20000
numKeep <- seq(numBurn+20, numSim, by=20)
a_samp <- mcmc.a[numKeep]
r_samp <- mcmc.r[numKeep]
sig_samp <- mcmc.sigma[numKeep]
tau_samp <- mcmc.tau[numKeep]
phi_samp <- mcmc.phi[numKeep]
mu_samp <- mcmc.mu[numKeep]
kappa_samp <- mcmc.kappa[numKeep]


#density
pdf("density.pdf")
par(mfrow=c(3,3))
plot(density(a_samp), xlab = expression(alpha),main="",y="",cex.lab=1.5)
plot(density(r_samp,from=0), xlab = "r",main="",y="",cex.lab=1.5)
plot(density(phi_samp), xlab = expression(phi),main="",y="",cex.lab=1.5)
plot(density(sig_samp,from=0), xlab = expression(sigma^2),main="",y="",cex.lab=1.5)
plot(density(tau_samp,from=0), xlab = expression(tau^2),main="",y="",cex.lab=1.5)
plot(density(mu_samp), xlab = expression(mu),main="",y="",cex.lab=1.5)
plot(density(kappa_samp), xlab = expression(nu),main="",y="",cex.lab=1.5)
dev.off()






# Krigging: marginalize out w
makePredictionAniso <- function(numSample, newSites){
  
  sigmasq <- sig_samp[numSample]
  tausq <- tau_samp[numSample]
  mu <- mu_samp[numSample]
  a <- a_samp[numSample]
  r <- r_samp[numSample]
  phi <- phi_samp[numSample]
  kappa <- kappa_samp[numSample]
  amin <- 1
  amax <- amin*r
  rotationMat <- matrix(c(cos(a),-sin(a),sin(a),cos(a)),nrow=2,ncol=2)
  aMat <- matrix(c(1/amax,0,0,1/amin),nrow=2,ncol=2)
  A <- rotationMat %*% aMat
  B <- A %*% t(A)
  SigmaObs <- cov.spatial(sqrt(makeB(B)), cov.model= "matern", cov.pars=c(sigmasq, 1/phi), kappa = kappa)
  diag(SigmaObs) <- sigmasq + tausq
  
  S <- solve(SigmaObs)
  
  predictions <- NULL
  for(j in 1:nrow(newSites)){
    newSite <- newSites[j,]
    distPO <- matrix(0, nrow=numSites, ncol = 2)
    for(i in 1:numSites){
      distPO[i,] <- newSite - sites_train[i,]
    }
    SigmaPO <- matrix(0, nrow=1, ncol=numSites)
    for(i in 1:numSites){
      SigmaPO[,i] <- cov.spatial(sqrt(distPO[i,] %*% B %*% distPO[i,]), cov.model= "matern", cov.pars=c(sigmasq, 1/phi), kappa = kappa)
    }
    mean_yPred <- mu + SigmaPO %*% S %*% (y_train - rep(mu, numSites))
    var_yPred <- sigmasq + tausq - SigmaPO %*% S %*% t(SigmaPO)
    yPred <- rnorm(1, mean_yPred, sqrt(var_yPred))
    predictions <- c(predictions, yPred)
  }
  
  return(predictions)
}
nMCMC <- length(a_samp)
predAniso <- t(sapply(1:nMCMC, makePredictionAniso, newSites=sites_test))



# Isotropy
library(spBayes)
fitIso <- spLM(y_train ~ 1, coords = sites_train,
               cov.model = "matern",
               n.samples = numSim,
               starting = list("sigma.sq" = 5,
                               "tau.sq" = 0.5,
                               "phi" = 0.03,
                               "nu" = 0.35), 
               tuning = list("phi" = 0.1,
                             "sigma.sq" = 0.1,
                             "tau.sq" = 0.1,
                             "nu" = 0.1),
               priors = list("sigma.sq.IG" = c(1, 1),
                             "tau.sq.IG" = c(1, 1),
                             "phi.Unif" = c(phi_prior0, phi_prior1),
                             "nu.Unif"=c(0.1,2)
               ), verbose = FALSE)

preds <- spPredict(fitIso, pred.covars=as.matrix(rep(1,nrow(sites_test))), pred.coords=sites_test,
                   start=numBurn+20, thin=20)
predIso <- t(preds$p.y.predictive.samples)






# 1. empirical coverage
EC <- function(mat, yObs){
  qt <- apply(mat, 2, function(x) quantile(x, probs = c(0.05, 0.95)))
  ave <- apply(mat, 2, mean)
  empCov <- data.frame(cbind(ave, t(qt), yObs))
  colnames(empCov) <- c("Mean", "Lower", "Higher", "True")
  empCov$capture <- empCov$Lower <= empCov$True & empCov$Higher >= empCov$True
  return(empCov)
}

empCovIso <- EC(predIso, y_test)
empCovAniso <- EC(predAniso, y_test)
numPred <- length(y_test)
ecIso <- sum(empCovIso$capture)/numPred
ecAniso <- sum(empCovAniso$capture)/numPred
# pdf("empCovIso.pdf")
# numPred <- length(y_test)
# ggplot(empCovAniso, aes(y=True, x=1:numPred, color=capture)) +
#   xlab("index") +
#   geom_errorbar(aes(ymax=Higher, ymin=Lower), width=0, color='black', alpha=0.3, size=2) +
#   geom_point(size=3) +
#   labs(x = "Index of New Sites", y="Predcited/Observed Value") +
#   ggtitle(paste("Empirical Coverage of Isotropic Model =",sum(empCovAniso$capture)/numPred))
# dev.off()



# 2. PMSE
mseIso <- mean((empCovIso$Mean-empCovIso$True)^2)
mseAniso <- mean((empCovAniso$Mean-empCovAniso$True)^2)


# 3. CRPS
crps_test = function(post,obs,n_pts = 1e6)
{
  F_post = ecdf(post)
  F_obs = ecdf(obs)
  
  d = c(obs,post)
  s = seq(min(d),max(d),len=n_pts)
  
  sum( (F_post(s) - F_obs(s))^2 ) * (max(d)-min(d)) / n_pts
}

crpsIso <- mean(sapply(1:ncol(predIso), function(x) crps_test(predIso[,x],y_test[x])))
crpsAniso <- mean(sapply(1:ncol(predAniso), function(x) crps_test(predAniso[,x],y_test[x])))
save()

save(mseIso, crpsIso, mseAniso, crpsAniso, predIso, predAniso, fitIso, preds, file="PM2.5.Rdata")


#load("PM2.5.Rdata")

#load("mcmc_pm_matern.Rdata")

## CRPS plot for isotropy
library(tidyr)


d_crps_iso = data_frame(
  prediction1 = predIso[,s[1]],
  prediction2 = predIso[,s[2]],
  prediction3 = predIso[,s[3]],
  prediction4 = predIso[,s[4]]
) %>% gather(dist)

yObs <- c(rep(y_test[s[1]],nrow(predIso)),rep(y_test[s[2]],nrow(predIso)),rep(y_test[s[3]],nrow(predIso)),rep(y_test[s[4]],nrow(predIso)))
name <- c("prediction1","prediction2","prediction3","prediction4")
df_yObs <- data.frame(dist=name, value=c(y_test[s[1]],y_test[s[2]],y_test[s[3]],y_test[s[4]]))
d_crps_iso <- cbind(d_crps_iso, yObs)
pmse <- d_crps_iso %>% group_by(dist) %>% summarize(value=mean(value), yObs=mean(yObs))
pmse$pmse <- round((pmse$value-pmse$yObs)^2,4)
pmse_lookup <- pmse$pmse %>% setNames(pmse$dist)
pmse_labeler = function(variable, value)
{
  paste0("pmse = ", pmse_lookup[value])
}


pdf("crps_Iso.pdf")
ggplot(d_crps_iso, aes(value, color=dist, fill=dist)) +
  geom_density(alpha=0.1) +
  facet_wrap(~dist, nrow=2, scales="free",labeller = pmse_labeler) +
  geom_vline(data=df_yObs, aes(xintercept=value), size=1) +
  theme(legend.title=element_blank()) +
  theme(text = element_text(size = 20)) + 
  theme(legend.position="none")
dev.off()

#### CRPS plot for anisotropy
d_crps_ani = data_frame(
  prediction1 = predAniso[,s[1]],
  prediction2 = predAniso[,s[2]],
  prediction3 = predAniso[,s[3]],
  prediction4 = predAniso[,s[4]]
) %>% gather(dist)

d_crps_ani <- cbind(d_crps_ani, yObs)
pmse <- d_crps_ani %>% group_by(dist) %>% summarize(value=mean(value), yObs=mean(yObs))
pmse$pmse <- round((pmse$value-pmse$yObs)^2,4)
pmse_lookup <- pmse$pmse %>% setNames(pmse$dist)
pmse_labeler = function(variable, value)
{
  paste0("pmse = ", pmse_lookup[value])
}

pdf("crps_Aniso.pdf", width=9)
ggplot(d_crps_ani, aes(value, color=dist, fill=dist)) +
  geom_density(alpha=0.1) +
  facet_wrap(~dist, nrow=2, scales="free",labeller = pmse_labeler) +
  geom_vline(data=df_yObs, aes(xintercept=value), size=1) +
  theme(legend.title=element_blank()) + 
  theme(text = element_text(size = 20)) + 
  xlim(0, 4)
dev.off()






