setwd("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/")            ## Server## breast cancer data
library(oro.nifti)

setwd("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/")            ## Server

# simulated data
load("save_simulation_123vpComp.Rdata")
# voxelwise fit
  name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/voxelwise/sim_Theta_voxelwise_021.Rdata"
  load(name)

I <- 75 # numer of rows
J <- 75 # numer of columns
N <- I*J # number of pixels

# ------------------------------------------------------------------------------
# functions
source("onePixFunctions.r")
source("spatialFunctions.r")

# Xmatrix
Matrix <- Cp(time, model="tofts.kermode") # First curve: plasma concentration
#Matrix <- c()     # without Cp
# basis function for different kep values
for (kep in exp(seq(-3,3,0.1)))
#for (kep in exp(seq(-5,3,0.05)))
{
    ft <- model.weinmann(time,kep,model="tofts.kermode")
    Matrix <- cbind(Matrix,ft)
}

# Define dimensions
T <- length(time)
Q <- ncol(Matrix) # number of basis functions, number of parameters

# ------------------------------------------------------------------------------
# data for nnent
sdX <- apply(Matrix,2,sd)
scaledX <- scale(Matrix,center=F,scale=sdX)

Response <- array(NA, c(I*J,T))
Fit <- array(NA, c(I*J,T))
Theta <- array(NA, c(I*J,Q))   # 2D array of Parameter values per Pixel

for(i in 1:I)
{
  for(j in 1:J)
  {
    Response[(i-1)*J+j,] <- CONC_SIM_NOISE[i,j,]
  }
}


# ------------------------------------------------------------------------------
# plot  maps
  
  require("fields")
  
  z_range <- c(-0.001,0.05)
  
  name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/voxelwise/Theta_map_021.pdf"
  
  pdf(name, width = 32, height = 32)
  par(mfrow=c(8,8), mar=c(1, 1, 1, 1), cex=1)
  #for(k in c(1:3,55:150))
  for(k in 1:Q)
  {
      Thetak_matrix <- t(matrix(Theta[,k], nrow=J))
      image(Thetak_matrix, axes=FALSE,zlim=z_range,col=gray(50:0/50),main=paste("k=",k))
  }
  dev.off()
     
  sse <- array(NA, I*J)
  error <- Response - Fit
  for(n in 1:N)
  {
  	 sse[n] <- sum(error[n,]^2)
  }
  
  name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/voxelwise/sse_voxelwise_021.pdf"
  pdf(name, width = 4, height = 4)
  par(mfrow=c(1,1), mar=c(1, 1, 1, 1), cex=1)
  
      sse_matrix <- t(matrix(sse, nrow=J))
      z_range <- c(0,0.15)
      image.plot(sse_matrix, axes=FALSE,zlim=z_range,col=gray(0:50/50),main="sse")
  
  dev.off()

# reduced Fit  
  name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/voxelwise/Theta2_map_021.pdf"
  
  pdf(name, width = 32, height = 32)
  par(mfrow=c(8,8), mar=c(1, 1, 1, 1), cex=1)
  #for(k in c(1:3,55:150))
  for(k in 1:Q)
  {
      Thetak_matrix <- t(matrix(Theta2[,k], nrow=J))
      image(Thetak_matrix, axes=FALSE,zlim=z_range,col=gray(50:0/50),main=paste("k=",k))
  }
  dev.off()
  
    
  sse <- array(NA, I*J)
  error <- Response - Fit2
  for(n in 1:N)
  {
  	 sse[n] <- sum(error[n,]^2)
  }
  
  name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/voxelwise/sse2_voxelwise_021.pdf"
  pdf(name, width = 4, height = 4)
  par(mfrow=c(1,1), mar=c(1, 1, 1, 1), cex=1)
  
      sse_matrix <- t(matrix(sse, nrow=J))
      z_range <- c(0,0.25)
      image(sse_matrix, axes=FALSE,zlim=z_range,col=gray(0:50/50))
      box()
        
  dev.off()
  
  name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/voxelwise/sse_legend.pdf"
  pdf(name, width = 1, height = 4)
  par(mfrow=c(1,1), mar=c(0.1, 0.1, 0.1, 0.1),  cex=1.4)
  
      sse_matrix <- t(matrix(sse, nrow=J))
      z_range <- c(0,0.25)
      image.plot(sse_matrix, legend.only = TRUE, axes=FALSE,zlim=z_range,col=gray(0:50/50))
        
  dev.off()
  
  
  
  # number of compartments
  number <- array(NA, I*J)
  number2 <- array(NA, I*J)
  number3 <- array(NA, I*J)
  vp_array <- array(NA, I*J)
  vp_array2 <- array(NA, I*J)
  for(n in 1:N)
  {
  	 #number[n] <- length(which(spatFit$selected[[n]]!=1))
  	 number[n] <- length(which(which(Theta[n,]>0)!=1))
  	 
  	 number2[n] <- length(which(which(Theta[n,]>10^(-5))!=1))
  
  	 vp_array[n] <- length(which(which(Theta[n,]!=0)==1))
  	 
  	 
  	 number3[n] <- length(which(which(Theta2[n,]>0)!=1))  
  	 vp_array2[n] <- length(which(which(Theta2[n,]!=0)==1))
  
  }
  
  name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/voxelwise/q_estimates_voxelwise_021.pdf" 
  pdf(name, width = 4, height = 4)
  par(mfrow=c(1,1), mar=c(1, 1, 1, 1), cex=1)
  
      number_matrix <- t(matrix(number, nrow=J))
      z_range <- c(0,50)
      image.plot(number_matrix, axes=FALSE,zlim=z_range,main="q",col=heat.colors(n=50, alpha = 1))
  
  dev.off()
  
  name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/voxelwise/q_rounded_voxelwise_021.pdf"
  pdf(name, width = 4, height = 4)
  par(mfrow=c(1,1), mar=c(1, 1, 1, 1), cex=1)
  
      number2_matrix <- t(matrix(number2, nrow=J))
      z_range <- c(0,10)
      image.plot(number2_matrix, axes=FALSE,zlim=z_range,main="q",col=heat.colors(n=50, alpha = 1))
  
  dev.off()
  
  name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/voxelwise/vp_estimates_voxelwise_021.pdf"
  pdf(name, width = 4, height = 4)
  par(mfrow=c(1,1), mar=c(1, 1, 1, 1), cex=1)
  
      vp_matrix <- t(matrix(vp_array, nrow=J))
      z_range <- c(0,1)
      image.plot(vp_matrix, axes=FALSE,zlim=z_range,main="vp",col=gray(0:50/50))
      box()
  
  dev.off()
  
  name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/voxelwise/q2_estimates_voxelwise_021.pdf"  
  pdf(name, width = 4, height = 4)
  par(mfrow=c(1,1), mar=c(1, 1, 1, 1), cex=1)
        
    number_matrix <- t(matrix(number3, nrow=J))
    z_range <- c(0,4)
    image(number_matrix, axes=FALSE,zlim=z_range,col=gray.colors(4, start=0.95, end=0))
    box()   

  dev.off()
  
  
  name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/voxelwise/vp2_estimates_voxelwise_021.pdf"
  pdf(name, width = 4, height = 4)
  par(mfrow=c(1,1), mar=c(1, 1, 1, 1), cex=1)
  
      vp_matrix <- t(matrix(vp_array2, nrow=J))
      z_range <- c(0,1)
      image(vp_matrix, axes=FALSE,zlim=z_range,col=gray.colors(2, start=1, end=0.3))
      box()
  
  dev.off()
  

load("save_simulation_123vpComp.Rdata")
PAR_VP[which(is.na(PAR_VP)==TRUE)] <- 0
  
 ########################
# Fit, selection...
name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/voxelwise/Fit_Selection_voxelwise_021.pdf"
pdf(name, width = 12, height = 8)
par(mfrow=c(2,3), mar=c(4, 4, 1, 1), cex=1)
y_range <- c(-0.01,0.16)

i <- 5  # q=1
j <- 20 # without vp

plot(1:Q, Theta[indx(c(i,j),I,J),], ylab=expression(theta),xlab=expression(k), type="h", ylim=y_range, main=paste("voxel",i,j))
true_k <- c(1,c(15,38,45)+1)
points(true_k, c(PAR_VP[i,j],PAR_KTRANS1[i,j],PAR_KTRANS2[i,j],PAR_KTRANS3[i,j])*sdX[true_k], type="h", col="blue",lwd=1.5)
points(1:Q,Theta2[indx(c(i,j),I,J),],col="gray", type="h",lwd=1.8)     

legend("topright", # x-y coordinates for location of the legend  
		cex=0.8,		# character expansion factor
		legend=c("True","Final estimates","First estimation step"),     # Legend labels  
		col=c("blue","gray","black"),   # Color of points or lines  
		lty=c(1,1,1),                    # Line type  
		lwd=c(1.5,1.8,1))  

#legend("topright", # x-y coordinates for location of the legend  
#		cex=0.8,		# character expansion factor
#		legend=c("True","Estimated"),     # Legend labels  
#		col=c("blue","black"),   # Color of points or lines  
#		lty=c(1,1),                    # Line type  
#		lwd=c(1.5,1))                    # Line width 

i <- 38  # q=2
j <- 20 # without vp

plot(1:Q, Theta[indx(c(i,j),I,J),], ylab=expression(theta),xlab=expression(k), type="h", ylim=y_range, main=paste("voxel",i,j))
true_k <- c(1,c(15,38,45)+1)
points(true_k, c(PAR_VP[i,j],PAR_KTRANS1[i,j],PAR_KTRANS2[i,j],PAR_KTRANS3[i,j])*sdX[true_k], type="h", col="blue",lwd=1.5)     
points(1:Q,Theta2[indx(c(i,j),I,J),],col="gray", type="h",lwd=1.8)     

legend("topright", # x-y coordinates for location of the legend  
		cex=0.8,		# character expansion factor
		legend=c("True","Final estimates","First estimation step"),     # Legend labels  
		col=c("blue","gray","black"),   # Color of points or lines  
		lty=c(1,1,1),                    # Line type  
		lwd=c(1.5,1.8,1))  

i <- 70  # q=3
j <- 20 # without vp

plot(1:Q, Theta[indx(c(i,j),I,J),], ylab=expression(theta),xlab=expression(k), type="h", ylim=y_range, main=paste("voxel",i,j))
true_k <- c(1,c(15,38,45)+1)
points(1:Q,Theta2[indx(c(i,j),I,J),],col="gray", type="h",lwd=1.8)     

legend("topright", # x-y coordinates for location of the legend  
		cex=0.8,		# character expansion factor
		legend=c("True","Final estimates","First estimation step"),     # Legend labels  
		col=c("blue","gray","black"),   # Color of points or lines  
		lty=c(1,1,1),                    # Line type  
		lwd=c(1.5,1.8,1))  


y_range <- c(-0.01,0.75)
i <- 5  # q=1
j <- 20 # without vp

plot(time, Response[indx(c(i,j),I,J),], ylab="C(t)", xlab="t",ylim=y_range)
points(time, CONC_SIM[i,j,], col="blue", type="l", lty=1, lwd=1.8)
points(time,Fit[indx(c(i,j),I,J),],col="black", type="l", lwd=1.8)
points(time,Fit2[indx(c(i,j),I,J),],col="gray", type="l", lwd=1.8)


i <- 38  # q=2
j <- 20 # without vp

plot(time, Response[indx(c(i,j),I,J),], ylab="C(t)", xlab="t",ylim=y_range)
points(time, CONC_SIM[i,j,], col="blue", type="l", lty=1, lwd=1.8)
points(time,Fit[indx(c(i,j),I,J),],col="black", type="l", lwd=1.8)
points(time,Fit2[indx(c(i,j),I,J),],col="gray", type="l", lwd=1.8)



i <- 70  # q=3
j <- 20 # without vp

plot(time, Response[indx(c(i,j),I,J),], ylab="C(t)", xlab="t", ylim=y_range)
points(time, CONC_SIM[i,j,], col="blue", type="l", lty=1, lwd=1.8)
points(time,Fit[indx(c(i,j),I,J),],col="black", type="l", lwd=1.8)
points(time,Fit2[indx(c(i,j),I,J),],col="gray", type="l", lwd=1.8)

dev.off()


################
# AIC

# einlesen spatial Fit
name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/spatFit_365.Rdata"
load(name)  
Fit2 <- spatFit$reducedFit
Theta2 <- spatFit$reducedTheta

# einlesen voxelwise Fit
#name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/voxelwise/sim_Theta_voxelwise_021.Rdata"
#load(name)

  q_voxelwise <- c()
  sse_voxelwise <- c()
  aic_voxelwise <- c()
  bic_voxelwise <- c()

  sse <- array(NA, I*J)
  error <- Response - Fit2
  number3 <- array(NA, I*J)
  for(n in 1:N)
  {
  	 sse[n] <- sum(error[n,]^2,na.rm=TRUE)
  	 #number3[n] <- length(which(which(Theta2[n,]>0)!=1))
  	 
  	 number3[n] <- length(which(Theta2[n,]>0))
  }
  
  true_sigma <- 0.05  # for simulated data
    
  q <- sum(number3,na.rm=TRUE) 
  aic <- sum(sse,na.rm=TRUE)/(true_sigma^2) + 2*q  
  bic <- sum(sse,na.rm=TRUE)/(true_sigma^2) + log(N*T) * q
  
  q_voxelwise <- c(q_voxelwise,q)
  sse_voxelwise <- c(sse_voxelwise,sum(sse,na.rm=TRUE))
  aic_voxelwise <- c(aic_voxelwise,aic)
  bic_voxelwise <- c(bic_voxelwise,bic)
  
  print(q_voxelwise)
  print(round(sse_voxelwise))
  print(round(aic_voxelwise))
  print(round(bic_voxelwise))
