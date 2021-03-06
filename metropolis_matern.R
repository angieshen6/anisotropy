##### metropolis Matern



##############################################

#This script compares isotropy with geometric 
# anisotropy with matern covariance


###############################################



library(MASS)
library(mvtnorm)
library(invgamma)
library(truncnorm)
library(spBayes)
library(geoR)


# ARGS <- commandArgs(trailingOnly = TRUE)
# numSites <- as.numeric(ARGS[1])
# numTotal <- as.numeric(ARGS[2])
# r <- as.numeric(ARGS[3])
# sigmasq <- as.numeric(ARGS[4])
# tausq <- as.numeric(ARGS[5])
omega <- 60
a <- omega/180*pi
mu <- 0
phi <- 3/0.5
kappa <- 3/2

mu <- 0
phi <- 3/0.5
r <- 16
omega <- 60
a <- omega/180*pi
sigmasq <- 1
tausq <- 0.2
numSites <- 500
numTotal <- 600


# simulate data

set.seed(1127)
sites <- cbind(runif(numTotal, 0, 1), runif(numTotal, 0, 1))
amin <- 1
amax <- amin*r
rotationMat <- matrix(c(cos(a),-sin(a),sin(a),cos(a)),nrow=2,ncol=2)
aMat <- matrix(c(1/amax,0,0,1/amin),nrow=2,ncol=2)
A <- rotationMat %*% aMat
B <- A %*% t(A)
Sigma <- matrix(0,nrow=numTotal, ncol=numTotal)
for(m in 1:numTotal){
  for(n in 1:numTotal){
    h <- sites[m,] - sites[n,]
    Sigma[m,n] <- cov.spatial(sqrt(t(h) %*% B %*% h), cov.model= "matern", cov.pars=c(sigmasq, 1/phi), kappa = kappa)
  }
}
Sigma <- Sigma + diag(tausq, nrow=numTotal, ncol=numTotal)
y <- mvrnorm(1, rep(0, numTotal), Sigma)
ind_test <- sample(1:numTotal, size=numTotal-numSites)
ind_train <- setdiff(1:numTotal, ind_test)
sites_test <- sites[ind_test,]
sites_train <- sites[ind_train,]
y_test <- y[ind_test]
y_train <- y[ind_train]

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




# initialize
numSim <- 30000
mcmc.sigma <- rep(sigmasq, numSim)
mcmc.a <- rep(a, numSim)
mcmc.r <- rep(r, numSim)
mcmc.tau <- rep(tausq, numSim)
mcmc.mu <- rep(mu, numSim)
mcmc.phi <- rep(phi, numSim)
mcmc.kappa <- rep(kappa, numSim)



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
  prior_star <- dunif(phi_star, 3/0.7, 3/0.2)
  prior_prev <- dunif(mcmc.phi[i-1], 3/0.7, 3/0.2)
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

save(mcmc.sigma, mcmc.tau, mcmc.a, mcmc.r, mcmc.phi, mcmc.kappa, mcmc.mu, file="mcmc_matern_16_500.Rdata")

numBurn <- 20000
numKeep <- seq(numBurn+20, numSim, by=20)
a_samp <- mcmc.a[numKeep]
r_samp <- mcmc.r[numKeep]
sig_samp <- mcmc.sigma[numKeep]
tau_samp <- mcmc.tau[numKeep]
phi_samp <- mcmc.phi[numKeep]
mu_samp <- mcmc.mu[numKeep]
kappa_samp <- mcmc.kappa[numKeep]


# 
# par(mfrow=c(2,3))
# plot(density(a_samp, from=0), xlab = "a")
# plot(density(r_samp, from=0), xlab = "r")
# plot(density(sig_samp, from=0), xlab = "phi")
# plot(density(phi_samp, from=0), xlab = "sigma")
# plot(density(tau_samp, from=0), xlab = "tau")
# plot(density(kappa_samp, from=0), xlab = "kappa")
# 


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

fitIso <- spLM(y_train ~ 1, coords = sites_train,
               cov.model = "matern",
               n.samples = numSim,
               starting = list("sigma.sq" = sigmasq,
                               "tau.sq" = tausq,
                               "phi" = phi,
                               "nu" = kappa), 
               tuning = list("phi" = 0.1,
                             "sigma.sq" = 0.1,
                             "tau.sq" = 0.1,
                             "nu" = 0.1),
               priors = list("sigma.sq.IG" = c(1, 1),
                             "tau.sq.IG" = c(1, 1),
                             "phi.Unif" = c(3/0.7, 3/0.2),
                             "nu.Unif"=c(0.1,2)
               ), verbose = FALSE)

preds <- spPredict(fitIso, pred.covars=as.matrix(rep(1,nrow(sites_test))), pred.coords=sites_test,
                   start=numBurn+20, thin=20)
predIso <- t(preds$p.y.predictive.samples)


#plot(1:40, y_test, pch=19, cex=0.5, xlab="observed y", ylab="predicted y",ylim=c(min(y.hat), max(y.hat)))
#arrows(1:40, y.hat[2,], 1:40, y.hat[1,], angle=90, length=0.05)
#arrows(1:40, y.hat[2,], 1:40, y.hat[3,], angle=90, length=0.05)





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


name <- paste(numSites, numTotal, r, sigmasq, tausq, "matern.Rdata", sep = "_")
save(ecIso, mseIso, crpsIso, ecAniso, mseAniso, crpsAniso, file=name)

# 
# load("100_140_1.5_1_0.2_matern.Rdata")
# load("100_140_4_1_0.2_matern.Rdata")
# load("100_140_8_0.2_0.04_matern.Rdata")
# load("100_140_8_0.2_0.2_matern.Rdata")
# load("100_140_8_1_0.2_matern.Rdata")
# load("100_140_8_1_1_matern.Rdata")
# 
# load("500_600_1.5_1_0.2_matern.Rdata")
# load("500_600_4_1_0.2_matern.Rdata")
# load("500_600_8_0.2_0.04_matern.Rdata")
# load("500_600_8_0.2_0.2_matern.Rdata")
# load("500_600_8_1_0.2_matern.Rdata")
# load("500_600_8_1_1_matern.Rdata")
# 
# 
# 
# (0.221-0.213)/0.221 
# (0.051-0.046)/0.051
# (1.153-1.101)/1.153
# (0.234-0.233)/0.233
# 
# (1.99-0.195)/0.195
# (0.19-0.175)/0.19
# (0.178-0.172)/0.178
# (0.048-0.041)/0.048
# (0.949-0.933)/0.949
# (0.19-0.175)/0.19
# 
# (0.266-0.262)/0.266
# 
# (0.273-0.271)/0.273
# (0.266-0.262)/0.262
# (0.13-0.123)/0.13
# (0.587-0.577)/0.587
# (0.248-0.237)/0.248
# (0.241-0.236)/0.241
# (0.124-0.114)/0.124
# (0.551-0.546)/0.551
# 
# 
# (0.278-0.261)/0.261
