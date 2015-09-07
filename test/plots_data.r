setwd("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/")            ## Server## breast cancer data
library(oro.nifti)

# quantiles vp estimates
for(s in c(1,2,3,4,5,6,7,8,9,10,11,12))
{ 

    var_s <- as.name(s)

    name <- paste("in scan ", var_s, sep=":")
    cat(name, fill=TRUE)

    # spatial fit
    name <- paste("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/spatFit_scan", var_s, sep="")
    name <- paste(name,"_47.Rdata",sep="")    
    load(name)
  
    print( quantile(spatFit$reducedTheta[which(spatFit$reducedTheta[,1]>0),1]/sdX[1],na.rm=TRUE, probs=c(0,0.25,0.5,0.75,0.9,0.95,1)))
        
}


for(s in c(1,2,3,4,5,6,7,8,9,10,11,12))
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

    # spatial fit
    name <- paste("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/spatFit_scan", var_s, sep="")
    name <- paste(name,"_47.Rdata",sep="")    
    load(name)
    
    require("fields")
    
    #z_range <- c(-0.001,0.05)
    #
    #name <- paste("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/Theta_map_scan", var_s, sep="")
    #name <- paste(name,"_47.pdf",sep="")
    #
    #pdf(name, width = 32, height = 32)
    #par(mfrow=c(8,8), mar=c(1, 1, 1, 1), cex=1)
    ##for(k in c(1:3,55:150))
    #for(k in 1:Q)
    #{
    #    Thetak_matrix <- t(matrix(spatFit$SpatialTheta[,k], nrow=J))
    #    image(Thetak_matrix, axes=FALSE,zlim=z_range,col=gray(50:0/50),main=paste("k=",k))
    #}
    #dev.off()
    
    
    z_range <- c(-0.001,0.1)
    
    name <- paste("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/reducedTheta_map_scan", var_s, sep="")
    name <- paste(name,"_47.pdf",sep="")
    pdf(name, width = 32, height = 32)
    par(mfrow=c(8,8), mar=c(1, 1, 1, 1), cex=1.8)
    for(k in 1:Q)
    {
        Thetak_matrix <- t(matrix(spatFit$reducedTheta[,k], nrow=J))
        image(Thetak_matrix, axes=FALSE,zlim=z_range,col=gray(50:0/50),main=paste("k=",k-1))
        box()
    }
    dev.off()
    
    
    z_range <- c(-0.001,0.1)
    
    name <- paste("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/reducedTheta_map2_scan", var_s, sep="")
    name <- paste(name,"_47.pdf",sep="")
    pdf(name, width = 28, height = 36)
    par(mfrow=c(9,7), mar=c(1, 1, 1, 1), cex=1.9)
    for(k in 1:Q)
    {
        Thetak_matrix <- t(matrix(spatFit$reducedTheta[,k], nrow=J))
        image(Thetak_matrix, axes=FALSE,zlim=z_range,col=gray(50:0/50),main=paste("k=",k-1))
        box()
    }
    dev.off()
    
    
    ##load("save_simulation_123vpComp.Rdata")
    ##PAR_VP[which(is.na(PAR_VP)==TRUE)] <- 0
    #name <- paste("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/Basis_Selection_scan", var_s, sep="")
    #name <- paste(name,"_47.pdf",sep="")
    #pdf(name, width = 20, height = 12)
    #par(mfrow=c(3,5), mar=c(4, 4, 1, 1), cex=1)
    #y_range <- c(-0.01,0.2)
    #for(i in c(5,35,60))
    #{
    #  for(j in c(5,20,40,50,60))
    #  {
    #
    #      plot(1:Q, spatFit$SpatialTheta[indx(c(i,j),I,J),], ylab=expression(K^{trans}), type="h", ylim=y_range, main=paste("voxel",i,j))
    #
    #      points(1:Q,spatFit$reducedTheta[indx(c(i,j),I,J),]/sdX,col="red", type="h")
    #      points(1:Q, spatFit$SpatialTheta[indx(c(i,j),I,J),], type="h") # nochmal drüber
    #
    #      #points(Selected_array[[indx(c(i,j),I,J)]],spatFit$SpatialTheta[indx(c(i,j),I,J),Selected_array[[indx(c(i,j),I,J)]]],col="red", type="h")
    #  }
    #}
    #dev.off()
    
    #sse <- array(NA, I*J)
    #error <- Response - spatFit$SpatialFit
    #for(n in 1:N)
    #{
    #	 sse[n] <- sum(error[n,]^2)
    #}
    #
    #name <- paste("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/sse_Fit_scan", var_s, sep="")
    #name <- paste(name,"_47.pdf",sep="")
    #pdf(name, width = 4, height = 4)
    #par(mfrow=c(1,1), mar=c(1, 1, 1, 1), cex=1)
    #
    #    sse_matrix <- t(matrix(sse, nrow=J))
    #    z_range <- c(0,0.15)
    #    image.plot(sse_matrix, axes=FALSE,zlim=z_range,col=gray(0:50/50),main="sse")
    #
    #dev.off()
    
    
    sse <- array(NA, I*J)
    error <- Response - spatFit$reducedFit
    for(n in 1:N)
    {
    	 sse[n] <- sum(error[n,]^2)
    }
    
    name <- paste("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/sse_scan", var_s, sep="")
    name <- paste(name,"_47.pdf",sep="")
    pdf(name, width = 4, height = 4)
    par(mfrow=c(1,1), mar=c(1, 1, 1, 1), cex=1)
    
        sse_matrix <- t(matrix(sse, nrow=J))
        z_range <- c(0,0.15)
        image.plot(sse_matrix, axes=FALSE,zlim=z_range,col=gray(0:50/50),main="sse")
    
    dev.off()
    
    
    name <- paste("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/1Comp/lsFit_scan", var_s, sep="")
    name <- paste(name,"_3.Rdata",sep="")    
    load(name)
    
    name <- paste("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/sse_compare_scan", var_s, sep="")
    name <- paste(name,"_47.pdf",sep="")
    
    z_range <- c(0,0.15)
    
    pdf(name, width = 4, height = 4)
    par(mfrow=c(1,1), mar=c(1, 1, 1, 1), cex=1)
        #image.plot(Fit1Comp$sse[,,1], axes=FALSE,zlim=z_range,col=gray(0:50/50),main="1Comp")        
        #image.plot(sse_matrix, axes=FALSE,zlim=z_range,col=gray(0:50/50),main="spatial")
        
        #z_range <- range(Fit1Comp$sse[,,1]-sse_matrix,na.rm=TRUE)
        
        z_range <- c(-0.05,0.07)
        
        image(Fit1Comp$sse[,,1]-sse_matrix, axes=FALSE,zlim=z_range,col=gray(0:50/50))
        box()
        
    dev.off()
    
#    name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/sse_compare_legend.pdf"    
#    pdf(name, width = 1, height = 4)
#    z_range <- c(-0.05,0.07)
#    #par(mfrow=c(1,1), mar=c(-0.1, -0.1, -0.1, 2),  cex=1.4)
#    par(mfrow=c(1,1), oma=c(0.1, 0.1, 0.1, 2),  cex=1.4)
#    
#    image.plot(sse_matrix, legend.only = TRUE, axes=FALSE, zlim=z_range,col=gray(0:50/50))
#    
#    dev.off()
    
    
    
    # number of compartments
    number <- array(NA, I*J)
    vp_array <- array(NA, I*J)
    for(n in 1:N)
    {
    	 #number[n] <- length(which(spatFit$selected[[n]]!=1))
    	 number[n] <- length(which(which(spatFit$reducedTheta[n,]>0)!=1))
    
    	 vp_array[n] <- length(which(which(spatFit$reducedTheta[n,]!=0)==1))
    
    	 #vp_array[n] <-  length(which(spatFit$selected[[n]]==1))
    }
    
    name <- paste("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/q_estimates_scan", var_s, sep="")
    name <- paste(name,"_47.pdf",sep="")
    pdf(name, width = 4, height = 4)
    par(mfrow=c(1,1), mar=c(1, 1, 1, 1), cex=1)
    
        number_matrix <- t(matrix(number, nrow=J))
        
        z_range <- c(0,4)
        image(number_matrix, axes=FALSE,zlim=z_range,col=gray.colors(4, start=0.95, end=0))
        box()
    
    dev.off()
    
    name <- paste("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/vp_estimates_scan", var_s, sep="")
    name <- paste(name,"_47.pdf",sep="")
    pdf(name, width = 4, height = 4)
    par(mfrow=c(1,1), mar=c(1, 1, 1, 1), cex=1)
    
        vp_matrix <- t(matrix(vp_array, nrow=J))
        z_range <- c(0,1)
        image(vp_matrix, axes=FALSE,zlim=z_range,col=gray.colors(2, start=1, end=0.3))
        box()
    
    dev.off()
        
    
    
#    #####   scan1
#    ## conc
#    name <- paste("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/conc_scan", var_s, sep="")
#    name <- paste(name,"_47.pdf",sep="")
#    pdf(name, width = 4, height = 4)
#    par(mfrow=c(1,1), mar=c(4, 4, 4, 4), cex=1)
#            
#        image(1:I,1:J,conc[,,10], main="",xlab="i",ylab="j")
#        
#        # scan01 normal tissue sourounding tumor
#        i <- 24
#        j <- 22    
#        #abline(h=j)
#        #abline(v=i)
#        points(i,j,cex=1)  
#        
#        # scan01 tumor edge
#        i <- 40
#        j <- 36 # oder 37
#        #abline(h=j)
#        #abline(v=i)
#        points(i,j,cex=1) 
#        
#        # scan01 inside the tumor
#        i <- 45
#        j <- 27 
#        #abline(h=j)
#        #abline(v=i)
#        points(i,j,cex=1) 
#    
#    dev.off()
#        
#    
#    ########################
#    # Fit, selection...
#    name <- paste("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/Fit_Selection_scan", var_s, sep="")
#    name <- paste(name,"_47.pdf",sep="")
#    pdf(name, width = 12, height = 8)
#    par(mfrow=c(2,3), mar=c(4, 4, 1, 1), cex=1)
#    y_range <- c(-0.01,0.1)
#    
#    # scan01 normal tissue sourounding tumor
#    i <- 24
#    j <- 22 
#    
#    plot(1:Q, spatFit$SpatialTheta[indx(c(i,j),I,J),], ylab=expression(theta),xlab=expression(k), type="h", ylim=y_range, main=paste("voxel",i,j))
#    points(1:Q,spatFit$reducedTheta[indx(c(i,j),I,J),],col="gray", type="h",lwd=1.8)
#    points(1:Q, spatFit$SpatialTheta[indx(c(i,j),I,J),], type="h") # nochmal drüber
#    legend("topright", # x-y coordinates for location of the legend  
#    		cex=0.8,		# character expansion factor
#    		legend=c("Final estimates","First estimation step"),     # Legend labels  
#    		col=c("gray","black"),   # Color of points or lines  
#    		lty=c(1,1),                    # Line type  
#    		lwd=c(1.8,1))                    # Line width 
#    
#    # scan01 tumor edge
#    i <- 40
#    j <- 36 # oder 37
#    
#    plot(1:Q, spatFit$SpatialTheta[indx(c(i,j),I,J),], ylab=expression(theta),xlab=expression(k), type="h", ylim=y_range, main=paste("voxel",i,j))
#    points(1:Q,spatFit$reducedTheta[indx(c(i,j),I,J),],col="gray", type="h",lwd=1.8)
#    points(1:Q, spatFit$SpatialTheta[indx(c(i,j),I,J),], type="h") # nochmal drüber
#    legend("topright", # x-y coordinates for location of the legend  
#    		cex=0.8,		# character expansion factor
#    		legend=c("Final estimates","First estimation step"),     # Legend labels  
#    		col=c("gray","black"),   # Color of points or lines  
#    		lty=c(1,1),                    # Line type  
#    		lwd=c(1.8,1))                    # Line width 
#    
#    # scan01 inside the tumor
#    i <- 45
#    j <- 27
#    
#    plot(1:Q, spatFit$SpatialTheta[indx(c(i,j),I,J),], ylab=expression(theta),xlab=expression(k), type="h", ylim=y_range, main=paste("voxel",i,j))
#    points(1:Q,spatFit$reducedTheta[indx(c(i,j),I,J),],col="gray", type="h",lwd=1.8)
#    points(1:Q, spatFit$SpatialTheta[indx(c(i,j),I,J),], type="h") # nochmal drüber
#    legend("topright", # x-y coordinates for location of the legend  
#    		cex=0.8,		# character expansion factor
#    		legend=c("Final estimates","First estimation step"),     # Legend labels  
#    		col=c("gray","black"),   # Color of points or lines  
#    		lty=c(1,1),                    # Line type  
#    		lwd=c(1.8,1))                    # Line width 
#    
#    
#    y_range <- c(-0.01,0.55)
#    # scan01 normal tissue sourounding tumor
#    i <- 24
#    j <- 22 
#    
#    plot(time, Response[indx(c(i,j),I,J),], ylab="C(t)", xlab="t",ylim=y_range)
#    points(time,spatFit$SpatialFit[indx(c(i,j),I,J),],col="gray", type="l", lwd=1.8)
#    
#    
#    # scan01 tumor edge
#    i <- 40
#    j <- 36 # oder 37
#    
#    plot(time, Response[indx(c(i,j),I,J),], ylab="C(t)", xlab="t",ylim=y_range)
#    points(time,spatFit$SpatialFit[indx(c(i,j),I,J),],col="gray", type="l", lwd=1.8)
#    
#    
#    # scan01 inside the tumor
#    i <- 45
#    j <- 27
#    
#    plot(time, Response[indx(c(i,j),I,J),], ylab="C(t)", xlab="t", ylim=y_range)
#    points(time,spatFit$SpatialFit[indx(c(i,j),I,J),],col="gray", type="l", lwd=1.8)
#    
#    dev.off()
    
    
####   scan03
 #####
    ## conc
    name <- paste("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/conc2_scan", var_s, sep="")
    name <- paste(name,"_47.pdf",sep="")
    pdf(name, width = 4, height = 4)
    par(mfrow=c(1,1), mar=c(4, 4, 0.1, 0.1), cex=1.45)
            
        image(1:I,1:J,conc[,,10], main="",xlab="i",ylab="j", col=gray.colors(50, start=1, end=0))
        
        # scan01 normal tissue sourounding tumor
        i <- 21
        j <- 20    
        #abline(h=j)
        #abline(v=i)
        points(i,j,cex=1)  
        
        # scan01 tumor edge
        i <- 35
        j <- 33 # oder 37
        #abline(h=j)
        #abline(v=i)
        points(i,j,cex=1) 
        
        # scan01 inside the tumor
        i <- 51
        j <- 41 
        #abline(h=j)
        #abline(v=i)
        points(i,j,cex=1) 
    
    dev.off()
    
    
    name <- paste("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/conc2_scan", var_s, sep="")
    name <- paste(name,"_47_legend.pdf",sep="")
    pdf(name, width = 1.0, height = 4.1)
	  par(mar=c(0.1, 0.1, 0.1, 0.1),  cex=1.45)
	  image.plot(1:I,1:J,conc[,,10],axes=FALSE,col=gray.colors(50, start=1, end=0), legend.only=TRUE)
	  dev.off()
    
# für Disputation
   name <- paste("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/conc3_scan", var_s, sep="")
    name <- paste(name,"_47.pdf",sep="")
    pdf(name, width = 5.5, height = 4)
    par(mfrow=c(1,1), mar=c(0.1, 0.1, 0.1, 0.1), cex=1)
        
        image.plot(1:I,1:J,conc[,,10], axes=FALSE, col=gray.colors(50, start=0, end=1))
    
    dev.off()
    

   
     ########################
    # Fit, selection...
    name <- paste("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/Fit2_Selection_scan", var_s, sep="")
    name <- paste(name,"_47.pdf",sep="")
    pdf(name, width = 12, height = 8)
    par(mfrow=c(2,3), mar=c(4, 4, 1, 1), cex=1.2)
    y_range <- c(0.001,0.18)
    
    # scan01 normal tissue sourounding tumor
    i <- 21
    j <- 20 
    
    plot(1:Q, spatFit$SpatialTheta[indx(c(i,j),I,J),], ylab=expression(theta),xlab=expression(k), col="grey",type="h", ylim=y_range, main=paste("voxel",i,j) ,lwd=1.3)
    points(1:Q,spatFit$reducedTheta[indx(c(i,j),I,J),],col="black", type="h",lwd=1.8)
    #points(1:Q, spatFit$SpatialTheta[indx(c(i,j),I,J),],col="grey", type="h") # nochmal drüber
    abline(0,0)
    legend("topright", # x-y coordinates for location of the legend  
    		cex=0.8,		# character expansion factor
    		legend=c("Final estimates","First estimation step"),     # Legend labels  
    		col=c("black","gray"),   # Color of points or lines  
    		lty=c(1,1),                    # Line type  
    		lwd=c(1.8,1.3))                    # Line width 
    
    # scan01 tumor edge
    i <- 35
    j <- 33 # oder 37
    
    plot(1:Q, spatFit$SpatialTheta[indx(c(i,j),I,J),], ylab=expression(theta),xlab=expression(k), col="grey", type="h", ylim=y_range, main=paste("voxel",i,j),lwd=1.3)
    points(1:Q,spatFit$reducedTheta[indx(c(i,j),I,J),],col="black", type="h",lwd=1.8)
    #points(1:Q, spatFit$SpatialTheta[indx(c(i,j),I,J),],col="grey", type="h") # nochmal drüber
    abline(0,0)
    legend("topright", # x-y coordinates for location of the legend  
    		cex=0.8,		# character expansion factor
    		legend=c("Final estimates","First estimation step"),     # Legend labels  
    		col=c("black","gray"),   # Color of points or lines  
    		lty=c(1,1),                    # Line type  
    		lwd=c(1.8,1.3))                    # Line width 
        
    # scan01 inside the tumor
    i <- 51
    j <- 41
    
   plot(1:Q, spatFit$SpatialTheta[indx(c(i,j),I,J),], ylab=expression(theta),xlab=expression(k), col="grey", type="h", ylim=y_range, main=paste("voxel",i,j),lwd=1.3)
    points(1:Q,spatFit$reducedTheta[indx(c(i,j),I,J),],col="black", type="h",lwd=1.8)
    #points(1:Q, spatFit$SpatialTheta[indx(c(i,j),I,J),],col="grey", type="h") # nochmal drüber
    abline(0,0)
    legend("topright", # x-y coordinates for location of the legend  
    		cex=0.8,		# character expansion factor
    		legend=c("Final estimates","First estimation step"),     # Legend labels  
    		col=c("black","gray"),   # Color of points or lines  
    		lty=c(1,1),                    # Line type  
    		lwd=c(1.8,1.3))                    # Line width 
        
    
    y_range <- c(-0.01,0.85)
    # scan01 normal tissue sourounding tumor
    i <- 21
    j <- 20 
    
    plot(time, Response[indx(c(i,j),I,J),], ylab="C(t)", col="black",xlab="t",ylim=y_range)
    points(time,spatFit$SpatialFit[indx(c(i,j),I,J),],col="black", type="l", lwd=1.8)
    
    
    # scan01 tumor edge
    i <- 35
    j <- 33 # oder 37
    
    plot(time, Response[indx(c(i,j),I,J),], ylab="C(t)", col="black",xlab="t",ylim=y_range)
    points(time,spatFit$SpatialFit[indx(c(i,j),I,J),],col="black", type="l", lwd=1.8)
    
    
    # scan01 inside the tumor
    i <- 51
    j <- 41
    
    plot(time, Response[indx(c(i,j),I,J),], ylab="C(t)", col="black",xlab="t",ylim=y_range)
    points(time,spatFit$SpatialFit[indx(c(i,j),I,J),],col="black", type="l", lwd=1.8)
    
    dev.off()
    
    
    ########################
    # Fit, selection...
    name <- paste("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/Fit2_Selection_simple_scan", var_s, sep="")
    name <- paste(name,"_47.pdf",sep="")
    pdf(name, width = 12, height = 8)
    par(mfrow=c(2,3), mar=c(4, 4, 1, 1), cex=1)
    y_range <- c(0,0.18)
    
    # scan01 normal tissue sourounding tumor
    i <- 21
    j <- 20 
    
    plot(1:Q, spatFit$reducedTheta[indx(c(i,j),I,J),], ylab=expression(theta),xlab=expression(k), type="h", ylim=y_range, main=paste("voxel",i,j),,lwd=1.8)
    abline(0,0)
    legend("topright", # x-y coordinates for location of the legend  
    		cex=0.8,		# character expansion factor
    		legend=c("Estimates"),     # Legend labels  
    		col=c("black"),   # Color of points or lines  
    		lty=c(1),                    # Line type  
    		lwd=c(1.8))                    # Line width 
    
    # scan01 tumor edge
    i <- 35
    j <- 33 # oder 37
    
    plot(1:Q, spatFit$reducedTheta[indx(c(i,j),I,J),], ylab=expression(theta),xlab=expression(k), type="h", ylim=y_range, main=paste("voxel",i,j),,lwd=1.8)
    abline(0,0)
    legend("topright", # x-y coordinates for location of the legend  
    		cex=0.8,		# character expansion factor
    		legend=c("Estimates"),     # Legend labels  
    		col=c("black"),   # Color of points or lines  
    		lty=c(1),                    # Line type  
    		lwd=c(1.8))                    # Line width 
        
    # scan01 inside the tumor
    i <- 51
    j <- 41
    
    plot(1:Q, spatFit$reducedTheta[indx(c(i,j),I,J),], ylab=expression(theta),xlab=expression(k), type="h", ylim=y_range, main=paste("voxel",i,j),,lwd=1.8)
    abline(0,0)
    legend("topright", # x-y coordinates for location of the legend  
    		cex=0.8,		# character expansion factor
    		legend=c("Estimates"),     # Legend labels  
    		col=c("black"),   # Color of points or lines  
    		lty=c(1),                    # Line type  
    		lwd=c(1.8))                    # Line width 
        
    
    y_range <- c(-0.01,0.85)
    # scan01 normal tissue sourounding tumor
    i <- 21
    j <- 20 
    
    plot(time, Response[indx(c(i,j),I,J),], ylab="C(t)", xlab="t",ylim=y_range)
    points(time,spatFit$SpatialFit[indx(c(i,j),I,J),],col="black", type="l", lwd=1.8)
    
    
    # scan01 tumor edge
    i <- 35
    j <- 33 # oder 37
    
    plot(time, Response[indx(c(i,j),I,J),], ylab="C(t)", xlab="t",ylim=y_range)
    points(time,spatFit$SpatialFit[indx(c(i,j),I,J),],col="black", type="l", lwd=1.8)
    
    
    # scan01 inside the tumor
    i <- 51
    j <- 41
    
    plot(time, Response[indx(c(i,j),I,J),], ylab="C(t)", xlab="t", ylim=y_range)
    points(time,spatFit$SpatialFit[indx(c(i,j),I,J),],col="black", type="l", lwd=1.8)
    
    dev.off()
    
    
    
    
    
    
    
    ####################
    ###
    ## scan06: large q estimate
    ######
    ### conc
    #name <- paste("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/conc_scan", var_s, sep="")
    #name <- paste(name,"_47.pdf",sep="")
    #pdf(name, width = 4, height = 4)
    #par(mfrow=c(1,1), mar=c(4, 4, 4, 4), cex=1)
    #        
    #    image(1:I,1:J,conc[,,10], main="",xlab="i",ylab="j")
    #    
    #    # scan06: q=15
    #    i <- 51
    #    j <- 31    
    #    #abline(h=j)
    #    #abline(v=i)
    #    points(i,j,cex=1)  
    #    
    #    # scan06: q=10
    #    i <- 52
    #    j <- 30    
    #    #abline(h=j)
    #    #abline(v=i)
    #    points(i,j,cex=1)
    #    
    #    # scan06: q=0
    #    i <- 50
    #    j <- 32 # oder 37
    #    #abline(h=j)
    #    #abline(v=i)
    #    points(i,j,cex=1) 
    #    
    #dev.off()
    #
    #
    
    
#    ########################
#    # Fit, selection...
#    name <- paste("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/Fit_Selection_scan", var_s, sep="")
#    name <- paste(name,"_47.pdf",sep="")
#    pdf(name, width = 12, height = 8)
#    par(mfrow=c(2,3), mar=c(4, 4, 1, 1), cex=1)
#    y_range <- c(-0.05,0.05)
#    
#    # scan08: q=5
#        i <- 4
#        j <- 88   
#    
#    plot(1:Q, spatFit$SpatialTheta[indx(c(i,j),I,J),], ylab=expression(theta),xlab=expression(k), type="h", ylim=y_range, main=paste("voxel",i,j))
#    points(1:Q,spatFit$reducedTheta[indx(c(i,j),I,J),],col="gray", type="h",lwd=1.8)
#    points(1:Q, spatFit$SpatialTheta[indx(c(i,j),I,J),], type="h") # nochmal drüber
#    legend("topright", # x-y coordinates for location of the legend  
#    		cex=0.8,		# character expansion factor
#    		legend=c("Final estimates","First estimation step"),     # Legend labels  
#    		col=c("gray","black"),   # Color of points or lines  
#    		lty=c(1,1),                    # Line type  
#    		lwd=c(1.8,1))                    # Line width 
#    
#    # scan08: q=5
#        i <- 70
#        j <- 64    
#    
#    plot(1:Q, spatFit$SpatialTheta[indx(c(i,j),I,J),], ylab=expression(theta),xlab=expression(k), type="h", ylim=y_range, main=paste("voxel",i,j))
#    points(1:Q,spatFit$reducedTheta[indx(c(i,j),I,J),],col="gray", type="h",lwd=1.8)
#    points(1:Q, spatFit$SpatialTheta[indx(c(i,j),I,J),], type="h") # nochmal drüber
#    legend("topright", # x-y coordinates for location of the legend  
#    		cex=0.8,		# character expansion factor
#    		legend=c("Final estimates","First estimation step"),     # Legend labels  
#    		col=c("gray","black"),   # Color of points or lines  
#    		lty=c(1,1),                    # Line type  
#    		lwd=c(1.8,1))                    # Line width 
#    
#    # scan08: q=5
#        i <- 7
#        j <- 72 # oder 37
#    
#    plot(1:Q, spatFit$SpatialTheta[indx(c(i,j),I,J),], ylab=expression(theta),xlab=expression(k), type="h", ylim=y_range, main=paste("voxel",i,j))
#    points(1:Q,spatFit$reducedTheta[indx(c(i,j),I,J),],col="gray", type="h",lwd=1.8)
#    points(1:Q, spatFit$SpatialTheta[indx(c(i,j),I,J),], type="h") # nochmal drüber
#    legend("topright", # x-y coordinates for location of the legend  
#    		cex=0.8,		# character expansion factor
#    		legend=c("Final estimates","First estimation step"),     # Legend labels  
#    		col=c("gray","black"),   # Color of points or lines  
#    		lty=c(1,1),                    # Line type  
#    		lwd=c(1.8,1))                    # Line width 
#    
#    
#    y_range <- c(-0.2,0.25)
#    # scan06: q=15
#        i <- 51
#        j <- 31   
#    
#    plot(time, Response[indx(c(i,j),I,J),], ylab="C(t)", xlab="t",ylim=y_range)
#    points(time,spatFit$SpatialFit[indx(c(i,j),I,J),],col="gray", type="l", lwd=1.8)
#    
#    
#    # scan06: q=10
#        i <- 52
#        j <- 30    
#    
#    plot(time, Response[indx(c(i,j),I,J),], ylab="C(t)", xlab="t",ylim=y_range)
#    points(time,spatFit$SpatialFit[indx(c(i,j),I,J),],col="gray", type="l", lwd=1.8)
#    
#    
#    # scan06: q=0
#        i <- 50
#        j <- 32 # oder 37
#    
#    plot(time, Response[indx(c(i,j),I,J),], ylab="C(t)", xlab="t", ylim=y_range)
#    points(time,spatFit$SpatialFit[indx(c(i,j),I,J),],col="gray", type="l", lwd=1.8)
#    
#    dev.off()
    
    
}

#########################################################################################

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
    name <- paste("~/Dateien/Daten_DCEMRI/breast/rois/roi", var_s, sep="")
    #name <- paste("~/Dateien/Daten_DCEMRI/breast/rois/roi", var_s, sep="")
    name <- paste(name,".txt",sep="")

    roi <- scan(name)

    dimx <- dim(x)
    mask <- array(FALSE, c(dimx[1], dimx[2], dimx[3]))

#    xmin <- min(roi[c(1,3,5,7)])
#    xmax <- max(roi[c(1,3,5,7)])
#    ymin <- min(roi[c(2,4,6,8)])
#    ymax <- max(roi[c(2,4,6,8)])
#
#    for(i in xmin : xmax) {
#    	for(j in ymin : ymax) {
#    		for(k in 1:dimx[3]){
#    			mask[i,j,k] <- TRUE
#    		}
#    	}
#    }
#
#    conc <- x[xmin : xmax, ymin : ymax,2,]
#    mask <- mask[xmin : xmax, ymin : ymax,2]
    
    conc <- x[,,2,]
    mask <- mask[,,2]

    print(sum(is.na(conc)))


    I <- nrow(conc)
    J <- ncol(conc)
    N <- I*J        # number of pixels



name <- paste("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/conc_scan", var_s, sep="")
    name <- paste(name,".pdf",sep="")
    pdf(name, width = 16, height = 8)
    par(mfrow=c(2,4), mar=c(1, 1, 1, 1), cex=1.8)
                    
        z_range <- c(-0.1,0.5)
        
        image(conc[,,1], axes=FALSE,zlim=z_range,col=gray(0:50/50),main=paste("t =",round(time[1]*60),"seconds"))
        box()
        
        image(conc[,,2], axes=FALSE,zlim=z_range,col=gray(0:50/50),main=paste("t =",round(time[2]*60),"seconds"))
        box()
        
        image(conc[,,3], axes=FALSE,zlim=z_range,col=gray(0:50/50),main=paste("t =",round(time[3]*60),"seconds"))
        box()
        
        image(conc[,,4], axes=FALSE,zlim=z_range,col=gray(0:50/50),main=paste("t =",round(time[4]*60),"seconds"))
        box()
        
        
        image(conc[,,8], axes=FALSE,zlim=z_range,col=gray(0:50/50),main=paste("t =",round(time[8]*60),"seconds"))
        box()
        
        image(conc[,,16], axes=FALSE,zlim=z_range,col=gray(0:50/50),main=paste("t =",round(time[16]*60),"seconds"))
        box()
        
        image(conc[,,24], axes=FALSE,zlim=z_range,col=gray(0:50/50),main=paste("t =",round(time[24]*60),"seconds"))
        box()
        
        image(conc[,,32], axes=FALSE,zlim=z_range,col=gray(0:50/50),main=paste("t =",round(time[32]*60),"seconds"))
        box()
                
dev.off()


name <- paste("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/conc_scan", var_s, sep="")
name <- paste(name,"_legend.pdf",sep="")
  pdf(name, width = 1, height = 4)
  par(mfrow=c(1,1), mar=c(0.1, 0.1, 0.1, 0.1),  cex=1.4)
  
      z_range <- c(-0.1,0.5)
      image.plot(conc[,,1], legend.only = TRUE, axes=FALSE,zlim=z_range,col=gray(0:50/50))
        
  dev.off()