setwd("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/")            ## Server## breast cancer data
library(oro.nifti)


for(s in c(1,2,3,4,5,6,7,8,9,10,11,12))
#for(s in 6)
{
    var_s <- as.name(s)


    name <- paste("in scan ", var_s, sep=":")
    cat(name, fill=TRUE)

    ## read scan data / image
    name <- paste("~/Dateien/Daten_DCEMRI/breast/scan", var_s, sep="")
    name <- paste(name,".img",sep="")

    #x <- read.img(name, verbose=FALSE, warn = -1)
    x <- readANALYZE(name, verbose=FALSE, warn = -1)

    ## read time array
    name <- paste("~/Dateien/Daten_DCEMRI/breast/time", var_s, sep="")
    name <- paste(name,".txt",sep="")

    time <- scan(name)

    ## fill mask-array for image, such that True only for pixels in the ROI
    #name <- paste("~/Dateien/Daten_DCEMRI/breast/rois_spatEnet/roi", var_s, sep="")
    name <- paste("~/Dateien/Daten_DCEMRI/breast/rois/roi", var_s, sep="")
    name <- paste(name,".txt",sep="")

    roi <- scan(name)

    dimx <- dim(x)
    mask <- array(FALSE, c(dimx[1], dimx[2], dimx[3]))

    xmin <- min(roi[c(1,3,5,7)])
    xmax <- max(roi[c(1,3,5,7)])
    ymin <- min(roi[c(2,4,6,8)])
    ymax <- max(roi[c(2,4,6,8)])

    for(i in xmin : xmax) {
    	for(j in ymin : ymax) {
    		for(k in 1:dimx[3]){
    			mask[i,j,k] <- TRUE
    		}
    	}
    }

    conc <- x[xmin : xmax, ymin : ymax,2,]
    mask <- mask[xmin : xmax, ymin : ymax,2]

    print(sum(is.na(conc)))


    I <- nrow(conc)
    J <- ncol(conc)
    N <- I*J        # number of pixels


    # ------------------------------------------------------------------------------
    # functions
    source("onePixFunctions.r")
    source("spatialFunctions.r")

    # Xmatrix
    Matrix <- Cp(time, model="tofts.kermode") # First curve: plasma concentration
    #Matrix <- c()     # without Cp
    # basis function for different kep values
    #for (kep in exp(seq(-2,2,1)))
    for (kep in exp(seq(-3,3,0.1)))
    {
        ft <- model.weinmann(time,kep,model="tofts.kermode")
        Matrix <- cbind(Matrix,ft)
    }

    # Define dimensions
    T <- length(time)
    Q <- ncol(Matrix) # number of basis functions, number of parameters

    # ------------------------------------------------------------------------------
    # data for nnent
    p <- ncol(Matrix)

    sdX <- apply(Matrix,2,sd)
    scaledX <- scale(Matrix,center=F,scale=sdX)

    xypix <- cbind(sort(rep(1:I,J)),rep(1:J,I))
    colnames(xypix) <- c("x","y")
    xysample <- xypix

    # checkerboard pattern
    xyblack <- xysample[which((xysample[,1]+xysample[,2])%%2==0, arr.ind=TRUE),1:2]
    xywhite <- xysample[which((xysample[,1]+xysample[,2])%%2==1, arr.ind=TRUE),1:2]

    xyblack_list <- c()
    xywhite_list <- c()
    for(i in 1:length(xywhite[,1]))
    {
           xyblack_list[[i]] <- xyblack[i,]
           xywhite_list[[i]] <- xywhite[i,]
    }


    # Response array (two-dimensional index)
    Response <- array(NA, c(I*J,T))
    Fit <- array(NA, c(I*J,T))
    Theta <- array(NA, c(I*J,Q))   # 2D array of Parameter values per Pixel

    for(i in 1:I)
    {
      for(j in 1:J)
      {
        Response[(i-1)*J+j,] <- conc[i,j,]
      }
    }

  # voxelwise fit
  name <- paste("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/voxelwise/voxelwiseFit_scan", var_s, sep="")
  name <- paste(name,"_10.Rdata",sep="")
  load(name)
  
  require("fields")
  
  z_range <- c(-0.001,0.1)
  
  name <- paste("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/voxelwise/Theta_map_scan", var_s, sep="")
  name <- paste(name,"_10.pdf",sep="")
  
  pdf(name, width = 32, height = 32)
  par(mfrow=c(8,8), mar=c(1, 1, 1, 1), cex=1)
  #for(k in c(1:3,55:150))
  for(k in 1:Q)
  {
      Thetak_matrix <- t(matrix(Theta[,k], nrow=J))
      image(Thetak_matrix, axes=FALSE,zlim=z_range,col=gray(50:0/50),main=paste("k=",k-1))
  }
  dev.off()
#    
#  sse <- array(NA, I*J)
#  error <- Response - Fit
#  for(n in 1:N)
#  {
#  	 sse[n] <- sum(error[n,]^2)
#  }
#  
#  name <- paste("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/voxelwise/sse_voxelwise_scan", var_s, sep="")
#  name <- paste(name,"_10.pdf",sep="")
#  pdf(name, width = 4, height = 4)
#  par(mfrow=c(1,1), mar=c(1, 1, 1, 1), cex=1)
#  
#      sse_matrix <- t(matrix(sse, nrow=J))
#      z_range <- c(0,0.15)
#      image.plot(sse_matrix, axes=FALSE,zlim=z_range,col=gray(0:50/50),main="sse")
#  
#  dev.off()
  
#  # two-step
#  z_range <- c(-0.001,0.05)
#  
#  name <- paste("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/voxelwise/Theta_map_scan", var_s, sep="")
#  name <- paste(name,"_10.pdf",sep="")
#  
#  pdf(name, width = 32, height = 32)
#  par(mfrow=c(8,8), mar=c(1, 1, 1, 1), cex=1)
#  #for(k in c(1:3,55:150))
#  for(k in 1:Q)
#  {
#      Thetak_matrix <- t(matrix(Theta2[,k], nrow=J))
#      image(Thetak_matrix, axes=FALSE,zlim=z_range,col=gray(50:0/50),main=paste("k=",k))
#  }
#  dev.off()
#    
  sse <- array(NA, I*J)
  error <- Response - Fit2
  for(n in 1:N)
  {
  	 sse[n] <- sum(error[n,]^2)
  }
  
  name <- paste("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/voxelwise/sse_voxelwise_scan", var_s, sep="")
  name <- paste(name,"_10.pdf",sep="")
  pdf(name, width = 4, height = 4)
  par(mfrow=c(1,1), mar=c(1, 1, 1, 1), cex=1)
  
      sse_matrix <- t(matrix(sse, nrow=J))
      z_range <- c(0,0.15)
      image.plot(sse_matrix, axes=FALSE,zlim=z_range,col=gray(0:50/50),main="sse")
  
  dev.off()
  
  
  # number of compartments
  number <- array(NA, I*J)
  number2 <- array(NA, I*J)
  vp_array <- array(NA, I*J)
  number3 <- array(NA, I*J)
  vp_array2 <- array(NA, I*J)
  for(n in 1:N)
  {
  	 #number[n] <- length(which(spatFit$selected[[n]]!=1))
  	 number[n] <- length(which(which(Theta[n,]>0)!=1))
  	 
  	 number2[n] <- length(which(which(Theta[n,]>10^(-5))!=1))
  
  	 vp_array[n] <- length(which(which(Theta[n,]!=0)==1))
  	 
  	 number3[n] <- length(which(which(Theta2[n,]>0)!=1))  
  	 vp_array2[n] <- length(which(which(Theta2[n,]!=0)==1))
  
  	 #vp_array[n] <-  length(which(spatFit$selected[[n]]==1))
  }
  
  
  ################
  # AIC
  
  # AIC = -2log-likeli+2q 
  sse <- array(NA, I*J)
  error <- Response - Fit2
  number3 <- array(NA, I*J)
  for(n in 1:N)
  {
  	 sse[n] <- sum(error[n,]^2)
  	 number3[n] <- length(which(which(Theta2[n,]>0)!=1))
  }
  
#  # load fullest model -> estimate for sigma^2
#  name <- paste("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/voxelwise/voxelwiseFit_scan", var_s, sep="")
#  name <- paste(name,"_10.Rdata",sep="")
#  load(name)
#  sse_fullest <- array(NA, I*J)
#  error_fullest <- Response - Fit
#  for(n in 1:N)
#  {
#  	 sse_fullest[n] <- sum(error_fullest[n,]^2)
#  }
#  
#  true_sigma <- 0.05
#  sigma2_esti <- sum(sse_fullest)/(N*T)
#  aic <- sum(sse)/(N*T*sigma2_esti) + 2*sum(number3)/(N*T)
#  
#  print("AIC")
#  print(aic)


  
#  name <- paste("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/voxelwise/q_estimates_voxelwise_scan", var_s, sep="")
#  name <- paste(name,"_10.pdf",sep="")
#  pdf(name, width = 4, height = 4)
#  par(mfrow=c(1,1), mar=c(1, 1, 1, 1), cex=1)
#  
#      number_matrix <- t(matrix(number, nrow=J))
#      z_range <- c(0,50)
#      image.plot(number_matrix, axes=FALSE,zlim=z_range,main="q",col=heat.colors(n=50, alpha = 1))
#  
#  dev.off()
#  
#  name <- paste("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/voxelwise/q_rounded_voxelwise_scan", var_s, sep="")
#  name <- paste(name,"_10.pdf",sep="")
#  pdf(name, width = 4, height = 4)
#  par(mfrow=c(1,1), mar=c(1, 1, 1, 1), cex=1)
#  
#      number2_matrix <- t(matrix(number2, nrow=J))
#      z_range <- c(0,10)
#      image.plot(number2_matrix, axes=FALSE,zlim=z_range,main="q",col=heat.colors(n=50, alpha = 1))
#  
#  dev.off()
#  
#  name <- paste("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/voxelwise/vp_estimates_voxelwise_scan", var_s, sep="")
#  name <- paste(name,"_10.pdf",sep="")
#  pdf(name, width = 4, height = 4)
#  par(mfrow=c(1,1), mar=c(1, 1, 1, 1), cex=1)
#  
#      vp_matrix <- t(matrix(vp_array, nrow=J))
#      z_range <- c(0,1)
#      image.plot(vp_matrix, axes=FALSE,zlim=z_range,main="vp",col=gray(0:50/50))
#      box()
#  
#  dev.off()
  
  name <- paste("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/voxelwise/q_estimates_voxelwise_scan", var_s, sep="")
  name <- paste(name,"_10.pdf",sep="")
  pdf(name, width = 4, height = 4)
  par(mfrow=c(1,1), mar=c(1, 1, 1, 1), cex=1)
  
      number_matrix <- t(matrix(number3, nrow=J))
      z_range <- c(0,4)
      image(number_matrix, axes=FALSE,zlim=z_range,col=gray.colors(4, start=0.95, end=0))
      box()     
      #z_range <- c(0,5)
      #image.plot(number_matrix, axes=FALSE,zlim=z_range,main="q",col=heat.colors(n=50, alpha = 1))
  
  dev.off()
  
  name <- paste("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/voxelwise/vp_estimates_voxelwise_scan", var_s, sep="")
  name <- paste(name,"_10.pdf",sep="")
  pdf(name, width = 4, height = 4)
  par(mfrow=c(1,1), mar=c(1, 1, 1, 1), cex=1)
  
      vp_matrix <- t(matrix(vp_array2, nrow=J))
      
      z_range <- c(0,1)
      image(vp_matrix, axes=FALSE,zlim=z_range,col=gray.colors(2, start=1, end=0.3))
      box()
      #z_range <- c(0,1)
      #image.plot(vp_matrix, axes=FALSE,zlim=z_range,main="vp",col=gray(0:50/50))
      #box()
  
  dev.off()
    
#  ########################
#  # Fit, selection...
#  name <- paste("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/voxelwise/Fit_voxelwise_scan", var_s, sep="")
#  name <- paste(name,"_10.pdf",sep="")
#  pdf(name, width = 12, height = 8)
#  par(mfrow=c(2,3), mar=c(4, 4, 1, 1), cex=1)
#  y_range <- c(-0.01,0.1)
#  
#  # scan01 normal tissue sourounding tumor
#  i <- 24
#  j <- 22 
#  
#  plot(1:Q, Theta[indx(c(i,j),I,J),], ylab=expression(theta),xlab=expression(k), type="h", ylim=y_range, main=paste("voxel",i,j))
#  points(1:Q,Theta2[indx(c(i,j),I,J),],col="gray", type="h",lwd=1.8)
#
#  legend("topright", # x-y coordinates for location of the legend  
#		cex=0.8,		# character expansion factor
#		legend=c("Final estimates","First estimation step"),     # Legend labels  
#		col=c("gray","black"),   # Color of points or lines  
#		lty=c(1,1),                    # Line type  
#		lwd=c(1.8,1))  
#
#  
#  # scan01 tumor edge
#  i <- 40
#  j <- 36 # oder 37
#  
#  plot(1:Q, Theta[indx(c(i,j),I,J),], ylab=expression(theta),xlab=expression(k), type="h", ylim=y_range, main=paste("voxel",i,j))
#  points(1:Q,Theta2[indx(c(i,j),I,J),],col="gray", type="h",lwd=1.8)
#  
#  legend("topright", # x-y coordinates for location of the legend  
#		cex=0.8,		# character expansion factor
#		legend=c("Final estimates","First estimation step"),     # Legend labels  
#		col=c("gray","black"),   # Color of points or lines  
#		lty=c(1,1),                    # Line type  
#		lwd=c(1.8,1))  
#  
#  # scan01 inside the tumor
#  i <- 45
#  j <- 27
#  
#  plot(1:Q, Theta[indx(c(i,j),I,J),], ylab=expression(theta),xlab=expression(k), type="h", ylim=y_range, main=paste("voxel",i,j))
#  points(1:Q,Theta2[indx(c(i,j),I,J),],col="gray", type="h",lwd=1.8)
#  legend("topright", # x-y coordinates for location of the legend  
#		cex=0.8,		# character expansion factor
#		legend=c("Final estimates","First estimation step"),     # Legend labels  
#		col=c("gray","black"),   # Color of points or lines  
#		lty=c(1,1),                    # Line type  
#		lwd=c(1.8,1))  
#  
#  
#  y_range <- c(-0.01,0.55)
#  # scan01 normal tissue sourounding tumor
#  i <- 24
#  j <- 22 
#  
#  plot(time, Response[indx(c(i,j),I,J),], ylab="C(t)", xlab="t",ylim=y_range)
#  points(time,Fit[indx(c(i,j),I,J),],col="black", type="l", lwd=1.8)
#  points(time,Fit2[indx(c(i,j),I,J),],col="gray", type="l", lwd=1.8)
#  
#  
#  # scan01 tumor edge
#  i <- 40
#  j <- 36 # oder 37
#  
#  plot(time, Response[indx(c(i,j),I,J),], ylab="C(t)", xlab="t",ylim=y_range)
#  points(time,Fit[indx(c(i,j),I,J),],col="black", type="l", lwd=1.8)
#  points(time,Fit2[indx(c(i,j),I,J),],col="gray", type="l", lwd=1.8)
#  
#  
#  # scan01 inside the tumor
#  i <- 45
#  j <- 27
#  
#  plot(time, Response[indx(c(i,j),I,J),], ylab="C(t)", xlab="t",ylim=y_range)
#  points(time,Fit[indx(c(i,j),I,J),],col="black", type="l", lwd=1.8)
#  points(time,Fit2[indx(c(i,j),I,J),],col="gray", type="l", lwd=1.8)
#  
#  dev.off()


}
