
##############################################

#This script creates the plots of ranges for 
# the PM 2.5 dataset. 

###############################################


load("PM2.5.Rdata")
save(a_samp, r_samp, mu_samp, phi_samp, sig_samp, tau_samp, kappa_samp, predAniso, file="samps_range.Rdata")
load("samps_range.Rdata")

kappa <- mean(kappa_samp)
#phi <- mean(phi_samp)
a <- mean(a_samp)
r <- mean(r_samp)



# function to calculate c given direction
calcRange <- function(direction){
  v <- c(cos(direction), sin(direction))
  amin <- 1
  amax <- amin*r
  rotationMat <- matrix(c(cos(a),-sin(a),sin(a),cos(a)),nrow=2,ncol=2)
  aMat <- matrix(c(1/amax,0,0,1/amin),nrow=2,ncol=2)
  A <- rotationMat %*% aMat
  B <- A %*% t(A)
  # exp(- r * sqrt(v %*% B %*% v)) = 0.05
  f <- function(x) {
    abs(as.vector(cov.spatial(x*sqrt(v %*% B %*% v), cov.model= "matern", cov.pars=c(1, 1), kappa = kappa) - 0.05))
  }
  
  golden.ratio = 2/(sqrt(5) + 1)
  
  lower.bound = 0
  upper.bound = 50
  tolerance = 0.0001
  ### Use the golden ratio to set the initial test points
  x1 = upper.bound - golden.ratio*(upper.bound - lower.bound)
  x2 = lower.bound + golden.ratio*(upper.bound - lower.bound)
  
  ### Evaluate the function at the test points
  f1 = f(x1)
  f2 = f(x2)
  
  iteration = 0
  
  while (abs(upper.bound - lower.bound) > tolerance)
  {
    iteration = iteration + 1
    if (f2 > f1)

    {

      upper.bound = x2
      x2 = x1
      f2 = f1
      
      x1 = upper.bound - golden.ratio*(upper.bound - lower.bound)
      f1 = f(x1)
    } 
    else 
    {
      lower.bound = x1
      x1 = x2
      
      f1 = f2
      
      x2 = lower.bound + golden.ratio*(upper.bound - lower.bound)
      f2 = f(x2)
    }
  }
  
  estimated.minimizer = (lower.bound + upper.bound)/2
  cat('Estimated Minimizer =', estimated.minimizer, '\n')
  
  
  return(estimated.minimizer)
}


# cDirect <- function(d){
#   c <- mapply(calcRange, direction=d, a=a_samp, r=r_samp, phi=phi_samp, kappa=kappa_samp)
#   res <- data.frame(lq=quantile(c, 0.025), ave=mean(c), hq=quantile(c, 0.975))
#   return(res)
# }


angles <- seq(0, 360, by=1)
radians <- angles/180*pi
ranges <- NULL
for(i in radians){
  ranges <- c(ranges, calcRange(i))
}


pdf("ranges.pdf",width=10,height=6)
par(mfrow=c(1,2))
plot(angles, ranges, ylim=c(min(ranges),max(ranges)), type="l", xlab="degrees", ylab="range",cex.lab=1.5)
#lines(angles, df$lq, col="red")
#lines(angles, df$hq, col="red")
plot(ranges*cos(radians), ranges*sin(radians), type="l", xlab="x", ylab="y", cex.lab=1.5,asp=1)
dev.off()



  