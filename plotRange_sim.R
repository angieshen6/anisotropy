######

library(tikzDevice)


phi <- 3/0.1

calcRange <- function(direction, a, r){
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





olist <- c(60, 150, 60, 150)
alist <- olist/180*pi
rlist <- c(1.5, 1.5, 4, 4)
angles <- seq(0, 360, by=1)
radians <- angles/180*pi

pdf("ranges_sim.pdf",width=10, height=12)
#par(mfrow=c(4,2),mar=c(5,6,4,4))
par(mfrow=c(4,2))
for(i in 1:4){
  ranges <- sapply(radians, calcRange, a=alist[i], r=rlist[i])
  plot(angles, ranges, type="l", 
       xlab=paste0("angle=",olist[i],", ratio=", rlist[i]), ylim=c(0,2),ylab="range",cex.lab=2, cex.axis = 2)
  plot(ranges*cos(radians), ranges*sin(radians), type="l", 
       xlab=paste0("angle=",olist[i],", ratio=", rlist[i]), ylab="range",cex.lab=2, cex.axis = 2,asp=1)
  
}
dev.off()







#pdf("ranges_sim_prod.pdf",width=10, height=12)
tikz('ranges_sim_prob.tex', standAlone = TRUE, width=10, height=10)

calcRange <- function(direction, phiX, phiY){
  
  v <- c(cos(direction), sin(direction))
  c <- (-log(0.05)) / (phiX * abs(v[1]) + phiY * abs(v[2]))
  return(c)
}


phiX_list <- c(5, 10, 10, 5)
phiY_list <- c(10, 5, 10, 5)
d_list <- seq(0, 1, by=0.01)
angles <- seq(0, 360, by=1)
radians <- angles/180*pi
angles <- seq(0, 360, by=1)

par(mfrow=c(4,2),mar=c(5,6,4,4))
for(i in 1:4){
  ranges <- sapply(radians, calcRange, phiX=phiX_list[i], phiY=phiY_list[i])
  plot(angles, ranges, type="l", 
       xlab=paste0("$\\phi_x=$",phiX_list[i],", $\\phi_y=$", phiY_list[i]), ylim=c(0,1),ylab="range",cex.axis = 2, cex.lab=2)
  plot(ranges*cos(radians), ranges*sin(radians), type="l", 
       xlab=paste0("$\\phi_x=$",phiX_list[i],", $\\phi_y=$", phiY_list[i]), ylab="range",cex.lab=2, cex.axis = 2, asp=1)
  
}
dev.off()


# Compile the tex file
tools::texi2dvi('ranges_sim_prob.tex',pdf=T)





