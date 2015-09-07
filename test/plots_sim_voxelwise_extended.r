setwd("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/")            ## Server## breast cancer data
library(oro.nifti)

setwd("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/")            ## Server

# simulated data
load("save_simulation_123vpComp_extended.Rdata")
# voxelwise fit
  name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/voxelwise/sim_Theta_voxelwise_extended_1.Rdata"
  load(name)

I <- 75 # numer of rows
J <- 75 # numer of columns
N <- I*J # number of pixels
K <- 30  # number of simulated images

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



# ------------------------------------------------------------------------------
# plot  maps
  
  require("fields")

  msse <- array(NA, I*J)   
  sse <- array(NA, c(I*J,K))
  
  error <- array(NA,c(I*J,T,K))
  for(k in 1:K)
  {
     error[,,k] <- Response[,,k] - Fit2[,,k]
  }
  
  #error <- Response - spatFit$SpatialFit
  for(n in 1:N)
  {
    	 msse[n] <-  sum(error[n,,]^2)/K
  }
  
  name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/voxelwise/sse_k_voxelwise_extended_1.pdf"
  pdf(name, width = 24, height = 24)  
  par(mfrow=c(6,6), mar=c(1, 1, 1, 1), cex=1)
  for(k in 1:K)
  {
    for(n in 1:N)
    {
    	 sse[n,k] <- sum(error[n,,k]^2)
    }
         
      sse_matrix <- t(matrix(sse[,k], nrow=J))
      z_range <- c(0,0.14)
      image.plot(sse_matrix, axes=FALSE,zlim=z_range,col=gray(0:50/50),main="sse")
      
  }
  dev.off()
  
  name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/voxelwise/msse_voxelwise_extended_1.pdf"
  pdf(name, width = 4, height = 4)
  par(mfrow=c(1,1), mar=c(1, 1, 1, 1), cex=1)
  
      sse_matrix <- t(matrix(msse, nrow=J))
      z_range <- c(0,0.25)
      image(sse_matrix, axes=FALSE,zlim=z_range,col=gray(0:50/50))
      box()
  
  dev.off()
  


#############  
   
  msse <- array(NA, I*J)   
  sse <- array(NA, c(I*J, K))
  error <- Response - Fit
  for(n in 1:N)
  {
    	 msse[n] <-  sum(error[n,,]^2)/K
  }
  
  name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/voxelwise/sse_k_voxelwise_extended.pdf"
  pdf(name, width = 24, height = 24)  
  par(mfrow=c(6,6), mar=c(1, 1, 1, 1), cex=1)
  for(k in 1:K)
  {
    for(n in 1:N)
    {
    	 sse[n,k] <- sum(error[n,,k]^2)
    }
         
      sse_matrix <- t(matrix(sse[,k], nrow=J))
      z_range <- c(0,0.38)
      image.plot(sse_matrix, axes=FALSE,zlim=z_range,col=gray(0:50/50),main="sse")
      
  }
  dev.off()
  
  name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/voxelwise/msse_voxelwise_extended.pdf"
  pdf(name, width = 4, height = 4)
  par(mfrow=c(1,1), mar=c(1, 1, 1, 1), cex=1)
  
      sse_matrix <- t(matrix(msse, nrow=J))
      z_range <- c(0,0.38)
      image.plot(sse_matrix, axes=FALSE,zlim=z_range,col=gray(0:50/50),main="msse")
  
  dev.off()


  
  
  # number of compartments
  number <- array(NA, c(I*J,K))
  mnumber <- array(NA, I*J)
  number2 <- array(NA, c(I*J,K))
  number3 <- array(NA, c(I*J,K))
  vp_array <- array(NA, c(I*J,K))
  mvp_array <- array(NA, I*J)
  vp_array2 <- array(NA, c(I*J,K))
  for(n in 1:N)
  {
    for(k in 1:K)
    {
    	 #number[n] <- length(which(spatFit$selected[[n]]!=1))
    	 number[n,k] <- length(which(which(Theta[n,,k]>0)!=1)) 
    	 
    	 number2[n,k] <- length(which(which(Theta[n,,k]>10^(-5))!=1))  
    
    	 vp_array[n,k] <- length(which(which(Theta[n,,k]!=0)==1))  
    	   	 
    	 number3[n,k] <- length(which(which(Theta2[n,,k]>0)!=1))   
    	 vp_array2[n,k] <- length(which(which(Theta2[n,,k]!=0)==1))  
    }
    mnumber[n] <- sum(number3[n,])/K
    mvp_array[n] <- sum(vp_array2[n,])/K
  }
  
  name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/voxelwise/mq2_estimates_voxelwise_extended_1.pdf" 
  pdf(name, width = 4, height = 4)
  par(mfrow=c(1,1), mar=c(1, 1, 1, 1), cex=1)
  
      number_matrix <- t(matrix(mnumber, nrow=J))
      z_range <- c(0,4)
      image(number_matrix, axes=FALSE,zlim=z_range,col=gray.colors(4, start=0.95, end=0))
      box()
  
  dev.off()
  
  name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/voxelwise/mvp2_estimates_voxelwise_extended_1.pdf"
  pdf(name, width = 4, height = 4)
  par(mfrow=c(1,1), mar=c(1, 1, 1, 1), cex=1)
  
      vp_matrix <- t(matrix(mvp_array, nrow=J))
      z_range <- c(0,1)
      image(vp_matrix, axes=FALSE,zlim=z_range,col=gray.colors(2, start=1, end=0.3))
      box()
  
  dev.off()
  
  
##############  
  name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/voxelwise/q_rounded_voxelwise_021.pdf"
  pdf(name, width = 4, height = 4)
  par(mfrow=c(1,1), mar=c(1, 1, 1, 1), cex=1)
  
      number2_matrix <- t(matrix(number2, nrow=J))
      z_range <- c(0,10)
      image.plot(number2_matrix, axes=FALSE,zlim=z_range,main="q",col=heat.colors(n=50, alpha = 1))
  
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
name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/voxelwise/Fit_Selection_voxelwise_extended.pdf"
pdf(name, width = 12, height = 8)
par(mfrow=c(2,3), mar=c(4, 4, 1, 1), cex=1)
y_range <- c(-0.01,0.16)

i <- 5  # q=1
j <- 20 # without vp

plot(1:Q, Theta[indx(c(i,j),I,J),,1], ylab=expression(theta),xlab=expression(k), type="h", ylim=y_range, main=paste("voxel",i,j))
true_k <- c(1,c(15,38,45)+1)
points(true_k, c(PAR_VP[i,j],PAR_KTRANS1[i,j],PAR_KTRANS2[i,j],PAR_KTRANS3[i,j])*sdX[true_k], type="h", col="blue",lwd=1.5)
points(1:Q,Theta2[indx(c(i,j),I,J),,1],col="gray", type="h",lwd=1.8)     

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

plot(1:Q, Theta[indx(c(i,j),I,J),,1], ylab=expression(theta),xlab=expression(k), type="h", ylim=y_range, main=paste("voxel",i,j))
true_k <- c(1,c(15,38,45)+1)
points(true_k, c(PAR_VP[i,j],PAR_KTRANS1[i,j],PAR_KTRANS2[i,j],PAR_KTRANS3[i,j])*sdX[true_k], type="h", col="blue",lwd=1.5)     
points(1:Q,Theta2[indx(c(i,j),I,J),,1],col="gray", type="h",lwd=1.8)     

legend("topright", # x-y coordinates for location of the legend  
		cex=0.8,		# character expansion factor
		legend=c("True","Final estimates","First estimation step"),     # Legend labels  
		col=c("blue","gray","black"),   # Color of points or lines  
		lty=c(1,1,1),                    # Line type  
		lwd=c(1.5,1.8,1))  

i <- 70  # q=3
j <- 20 # without vp

plot(1:Q, Theta[indx(c(i,j),I,J),,1], ylab=expression(theta),xlab=expression(k), type="h", ylim=y_range, main=paste("voxel",i,j))
true_k <- c(1,c(15,38,45)+1)
points(1:Q,Theta2[indx(c(i,j),I,J),,1],col="gray", type="h",lwd=1.8)     

legend("topright", # x-y coordinates for location of the legend  
		cex=0.8,		# character expansion factor
		legend=c("True","Final estimates","First estimation step"),     # Legend labels  
		col=c("blue","gray","black"),   # Color of points or lines  
		lty=c(1,1,1),                    # Line type  
		lwd=c(1.5,1.8,1))  


y_range <- c(-0.01,0.75)
i <- 5  # q=1
j <- 20 # without vp

plot(time, Response[indx(c(i,j),I,J),,1], ylab="C(t)", xlab="t",ylim=y_range)
points(time, CONC_SIM[i,j,], col="blue", type="l", lty=1, lwd=1.8)
points(time,Fit[indx(c(i,j),I,J),,1],col="black", type="l", lwd=1.8)
points(time,Fit2[indx(c(i,j),I,J),,1],col="gray", type="l", lwd=1.8)


i <- 38  # q=2
j <- 20 # without vp

plot(time, Response[indx(c(i,j),I,J),,1], ylab="C(t)", xlab="t",ylim=y_range)
points(time, CONC_SIM[i,j,], col="blue", type="l", lty=1, lwd=1.8)
points(time,Fit[indx(c(i,j),I,J),,1],col="black", type="l", lwd=1.8)
points(time,Fit2[indx(c(i,j),I,J),,1],col="gray", type="l", lwd=1.8)



i <- 70  # q=3
j <- 20 # without vp

plot(time, Response[indx(c(i,j),I,J),,1], ylab="C(t)", xlab="t", ylim=y_range)
points(time, CONC_SIM[i,j,], col="blue", type="l", lty=1, lwd=1.8)
points(time,Fit[indx(c(i,j),I,J),,1],col="black", type="l", lwd=1.8)
points(time,Fit2[indx(c(i,j),I,J),,1],col="gray", type="l", lwd=1.8)

dev.off()



###########################
### für Disputation
########################
# Fit, selection...
name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/voxelwise/Fit_Selection_365_single.pdf"
pdf(name, width = 4, height = 3)
par(mfrow=c(1,1), mar=c(4, 4, 1, 1), cex=1)
y_range <- c(-0.001,0.035)

i <- 38  # q=2
j <- 20 # without vp

plot(1:Q, Theta[indx(c(i,j),I,J),,1], ylab=expression(theta),xlab=expression(k), type="h", ylim=y_range)

sel <- rep(NA,Q)
sel[1] <- Theta[indx(c(i,j),I,J),1,1]
sel[18] <- Theta[indx(c(i,j),I,J),18,1]
sel[42] <- Theta[indx(c(i,j),I,J),42,1]
points(1:Q, sel, col="blue", type="h",lwd=1.8)
true_k <- c(1,c(15,38,45)+1)
#points(true_k, c(PAR_VP[i,j],PAR_KTRANS1[i,j],PAR_KTRANS2[i,j],PAR_KTRANS3[i,j])*sdX[true_k], type="h", col="blue",lwd=1.5)     
#points(1:Q,Theta2[indx(c(i,j),I,J),,1],col="gray", type="h",lwd=1.8)     

#legend("topright", # x-y coordinates for location of the legend  
#		cex=0.8,		# character expansion factor
#		legend=c("True","Final estimates","First estimation step"),     # Legend labels  
#		col=c("blue","gray","black"),   # Color of points or lines  
#		lty=c(1,1,1),                    # Line type  
#		lwd=c(1.5,1.8,1))  

dev.off()


################
# AIC

# einlesen spatial Fit
load("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/spatFit_extended_1.Rdata")

# einlesen voxelwise Fit
  #name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/voxelwise/sim_Theta_voxelwise_021.Rdata"
  #load(name)
name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/voxelwise/sim_Theta_voxelwise_extended_1.Rdata"
load(name)


q <- array(NA,K)
aic <- array(NA,K)
bic <- array(NA,K)

for(k in 1:K)
{
    # bei spatial Fit
    Fit2_k <- spatFit_list[[k]]$SpatialFit
    Theta2_k <- spatFit_list[[k]]$reducedTheta
  
    # bei voxelwise Fit
    #Fit2_k <- Fit2[,,k]
    #Theta2_k <- Theta2[,,k]
  
  
#    q_voxelwise <- c()
#    sse_voxelwise <- c()
#    aic_voxelwise <- c()
#    bic_voxelwise <- c()
  
    sse <- array(NA, I*J)
    error <- Response[,,k] - Fit2_k
    number3 <- array(NA, I*J)
    for(n in 1:N)
    {
    	 sse[n] <- sum(error[n,]^2,na.rm=TRUE)
    	 #number3[n] <- length(which(which(Theta2[n,]>0)!=1))
    	 
    	 number3[n] <- length(which(Theta2_k[n,]>0))
    }
    
    true_sigma <- 0.05  # for simulated data
      
    q[k] <- sum(number3,na.rm=TRUE) 
    aic[k] <- sum(sse,na.rm=TRUE)/(true_sigma^2) + 2*q[k]  
    bic[k] <- sum(sse,na.rm=TRUE)/(true_sigma^2) + log(N*T) * q[k]
}
  
#  q_voxelwise <- c(q_voxelwise,q)
#  sse_voxelwise <- c(sse_voxelwise,sum(sse,na.rm=TRUE))
#  aic_voxelwise <- c(aic_voxelwise,aic)
#  bic_voxelwise <- c(bic_voxelwise,bic)
#  
#  print(q_voxelwise)
#  print(round(sse_voxelwise))
#  print(round(aic_voxelwise))
#  print(round(bic_voxelwise))
#