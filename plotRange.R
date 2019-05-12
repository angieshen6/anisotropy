

load("samps_range.Rdata")

# function to calculate c given direction
calcRange <- function(direction, a, r, phi){
  v <- c(cos(direction), sin(direction))
  amin <- 1
  amax <- amin*r
  rotationMat <- matrix(c(cos(a),-sin(a),sin(a),cos(a)),nrow=2,ncol=2)
  aMat <- matrix(c(1/amax,0,0,1/amin),nrow=2,ncol=2)
  A <- rotationMat %*% aMat
  B <- A %*% t(A)
  # exp(- r * sqrt(v %*% B %*% v)) = 0.05
  c <- (-log(0.05)) / (phi * sqrt(t(v) %*% B %*% v))
  return(c)
}
  

cDirect <- function(d){
  c <- mapply(calcRange, direction=d, a=a_samp, r=r_samp, phi=phi_samp)
  res <- data.frame(lq=quantile(c, 0.025), ave=mean(c), hq=quantile(c, 0.975))
  return(res)
}


angles <- seq(0, 360, by=5)
radians <- angles/180*pi
df <- NULL
for(i in radians){
  df <- rbind(df, cDirect(i))
}
rownames(df) <- angles


pdf("ranges.pdf",width=10,height=6)
par(mfrow=c(1,2))
plot(angles, df$ave, ylim=c(min(df$ave),max(df$ave)), type="l", xlab="degrees", ylab="range",cex.lab=1.5)
#lines(angles, df$lq, col="red")
#lines(angles, df$hq, col="red")
plot(df$ave*cos(radians), df$ave*sin(radians), type="l", xlab="x", ylab="y", cex.lab=1.5,asp=1)
dev.off()

save(a_samp, r_samp, mu_samp, phi_samp, sig_samp, tau_samp, kappa_samp, predAniso, file="samps_range.Rdata")


  