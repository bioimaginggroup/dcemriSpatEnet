setwd("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/")            ## Server## breast cancer data
library(oro.nifti)

q_spatial <- c()
sse_spatial <- c()
aic_spatial <- c()
bic_spatial <- c()

q_1Comp <- c()
sse_1Comp <- c()
aic_1Comp <- c()
bic_1Comp <- c()

q_voxelwise <- c()
sse_voxelwise <- c()
aic_voxelwise <- c()
bic_voxelwise <- c()
N_list <- c()
N2_list <- c()


#for(s in c(1,2,3))
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
    
    # number of voxels with measurement (not NaN)
    N2 <- N-sum(is.na(Response[,2]))


  # load fit to compute AIC
  # spatial
#  name <- paste("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatFit_scan", var_s, sep="")
#  name <- paste(name,"_12.Rdata",sep="")
#  load(name)
#
#  Fit2 <- spatFit$reducedFit
#  Theta2 <- spatFit$reducedTheta
  
  
  
  # load fullest model -> estimate for sigma^2
  name <- paste("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/voxelwise/voxelwiseFit_scan", var_s, sep="")
  name <- paste(name,"_8.Rdata",sep="")
  load(name)
  

  sse_fullest <- array(NA, I*J)
  error_fullest <- Response - Fit
  for(n in 1:N)
  {
  	 sse_fullest[n] <- sum(error_fullest[n,]^2,na.rm=TRUE)
  }
  sigma2_esti <- sum(sse_fullest,na.rm=TRUE)/(N2*T)
  #true_sigma <- 0.05  # for simulated data
  


  #######################################
  # voxelwise
  name <- paste("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/voxelwise/voxelwiseFit_scan", var_s, sep="")
  name <- paste(name,"_10.Rdata",sep="")
  load(name)

  sse <- array(NA, I*J)
  error <- Response - Fit2
  number3 <- array(NA, I*J)
  for(n in 1:N)
  {
  	 sse[n] <- sum(error[n,]^2)
  	 #number3[n] <- length(which(which(Theta2[n,]>0)!=1))
  	 number3[n] <- length(which(Theta2[n,]>0))
  }
 
  # AIC = -2log-likeli+2q
  #aic <- sum(sse,na.rm=TRUE)/(N*T*sigma2_esti) + 2*sum(number3,na.rm=TRUE)/(N*T)
     
  q <- sum(number3,na.rm=TRUE)
  aic <- sum(sse,na.rm=TRUE)/(sigma2_esti) + 2*q  
  bic <- sum(sse,na.rm=TRUE)/(sigma2_esti) + log(N2*T) * q
  
  q_voxelwise <- c(q_voxelwise,q)
  sse_voxelwise <- c(sse_voxelwise,sum(sse,na.rm=TRUE))
  aic_voxelwise <- c(aic_voxelwise,aic)
  bic_voxelwise <- c(bic_voxelwise,bic)
  
  N_list <- c(N_list,N)
  N2_list <- c(N2_list,N2)

  
#  #######################################
#  # extended 1Comp
#  name <- paste("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/1Comp/lsFit_scan", var_s, sep="")
#  name <- paste(name,"_3.Rdata",sep="")
#  load(name)
#
#  sse <- Fit1Comp$sse[,,1]
#  
# 
#  # AIC = -2log-likeli+2q  
#  q <- 2*N2 # fixed model, vp + ktrans1 per voxel
#  aic <- sum(sse,na.rm=TRUE)/(sigma2_esti) + 2*q  
#  bic <- sum(sse,na.rm=TRUE)/(sigma2_esti) + log(N2*T) * q
#  
#  q_1Comp <- c(q_1Comp,q)
#  sse_1Comp <- c(sse_1Comp,sum(sse,na.rm=TRUE))
#  aic_1Comp <- c(aic_1Comp,aic)
#  bic_1Comp <- c(bic_1Comp,bic)
#  N_list <- c(N_list,N)
#  N2_list <- c(N2_list,N2)
  
  
  
#  #######################################
#  # spatial
#  name <- paste("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/spatFit_scan", var_s, sep="")
#  name <- paste(name,"_47.Rdata",sep="")
#  load(name)
#
#  sse <- array(NA, I*J)
#  error <- Response - spatFit$reducedFit
#  number3 <- array(NA, I*J)
#  for(n in 1:N)
#  {
#  	 sse[n] <- sum(error[n,]^2)
#  	 #number3[n] <- length(which(which(spatFit$reducedTheta[n,]>0)!=1))
#  	 number3[n] <- length(which(spatFit$reducedTheta[n,]>0))
#  }
# 
#  # AIC = -2log-likeli+2q  
#  q <- sum(number3,na.rm=TRUE)
#  aic <- sum(sse,na.rm=TRUE)/(sigma2_esti) + 2*q  
#  bic <- sum(sse,na.rm=TRUE)/(sigma2_esti) + log(N2*T) * q
#  
#  q_spatial <- c(q_spatial,q)
#  sse_spatial <- c(sse_spatial,sum(sse,na.rm=TRUE))
#  aic_spatial <- c(aic_spatial,aic)
#  bic_spatial <- c(bic_spatial,bic)
#  N_list <- c(N_list,N)
#  N2_list <- c(N2_list,N2)
#  
  
  
    #print("AIC")
    #  print(aic)
    #  print(sigma2_esti)
    #  print(sum(number3)/(N))
    #print("BIC")
    #print(bic)
    
}

#print(q_spatial)
#print(round(sse_spatial))
#print(round(aic_spatial))
#print(round(bic_spatial))
#print(N_list)
#

print(q_voxelwise)
print(round(sse_voxelwise))
print(round(aic_voxelwise))
print(round(bic_voxelwise))
print(N_list)
print(N2_list)


#print(q_1Comp)
#print(round(sse_1Comp))
#print(round(aic_1Comp))
#print(round(bic_1Comp))


