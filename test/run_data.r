setwd("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/")            ## Server## breast cancer data
library(oro.nifti)
#for(s in c(9,10))
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
    
#    count <- 0
#    for(n in 1:N)
#    {
#    
#      if(sum(Response[n,]<0)>T/2)
#      {
#          #print("negative concentration in voxel")
#          #print(n)
#          count <- count + 1
#          
#          Response[n,] <- rep(NA,T)
#      }
#    
#    }
#
### Response
#name <- paste("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/Response_scan", var_s, sep="")
#name <- paste(name,"_na.pdf",sep="")
#pdf(name, width = 4, height = 4)
#par(mfrow=c(1,1), mar=c(4, 4, 4, 4), cex=1)
#    
#    
#    Response_matrix <- t(matrix(Response[,10], nrow=J))    
#    image(1:I,1:J,Response_matrix[,], main="",xlab="i",ylab="j")
#
#dev.off()
    
    
    
    # ------------------------------------------------------------------------------
    # Voxelwise Fit
    Theta <- array(NA, c(I*J,Q))
    Fit <- array(NA, c(I*J,T))
    Theta2 <- array(0, c(I*J,Q))
    Fit2 <- array(NA, c(I*J,T))
    
    for (pix in 1:length(xysample[,1]))
    {
        single_index <- xysample[pix,]
         
        if(!is.na(Response[indx(single_index,I,J),1]))
        {
          test <- nnenet2(scaledX, Response[indx(single_index,I,J),], lambda=10^(-7), 5, firstpen=TRUE)
          Theta[indx(single_index,I,J),] <- test$beta
          Fit[indx(single_index,I,J),] <- test$fit
          Theta2[indx(single_index,I,J),test$selected] <- test$reducedBeta
          Fit2[indx(single_index,I,J),] <- test$reducedFit
          
        }  
    }
    
    name <- paste("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/voxelwise/voxelwiseFit_scan", var_s, sep="")
    name <- paste(name,"_10.Rdata",sep="")
    
    save(Response,scaledX,Theta,Fit,Theta2,Fit2,file=name)
    
   
    
#    # ------------------------------------------------------------------------------
#    # starting value for spatial Fit
#    # Mean Response on ROI level
#    
#    conc_zeros <- conc
#    conc_zeros[which(is.na(conc))] <- 0
#    
#    summed_curve <- array(0, T)
#    for(i in 1:I)
#    {
#      for(j in 1:J)
#      {
#         summed_curve <- summed_curve +  conc_zeros[i,j,]
#      }
#    }
#    MeanResponse <- array(NA, c(I*J,T))
#    for(i in 1:I)
#    {
#      for(j in 1:J)
#      {
#          MeanResponse[(i-1)*J+j,] <- summed_curve/N
#      }
#    }
#    
#    test <- nnenet(scaledX, MeanResponse[1,], lambda=0.0001, 0.125, firstpen=FALSE)
#    
#    MeanTheta <- array(NA, c(I*J,Q))
#    MeanFit <- array(NA, c(I*J,T))
#    for (pix in 1:length(xysample[,1]))
#    {
#          MeanTheta[pix,] <- test$beta
#          MeanFit[pix,] <- test$fit
#    }
#    
#    #save(Response,scaledX,Theta,Fit,file="~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/scan03_MeanTheta_startingvalue_11.Rdata")
#        
#    # ------------------------------------------------------------------------------
#    # spatial fit for fixed lambda and sv
#    #
#    library(multicore,lib="~/R/x86_64-pc-linux-gnu-library/2.12/")
#    options(cores=10)
#    multicore=TRUE
#    
#    sv2 <- 5
#    lambda2 <- 10^(-10)
#    spatFit <- spatialFit(xyblack_list, xywhite_list, I=I, J=J, x=scaledX, Response=Response, time=time,
#                Theta=MeanTheta, lambda=lambda2, sv=sv2, firstpen=TRUE, centersOnly=TRUE)
#    
#    
#    name <- paste("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/spatFit_scan", var_s, sep="")
#    name <- paste(name,"_47.Rdata",sep="")
#                
#    save(Response,scaledX,time,spatFit,file=name)
#

}

#    ## ------------------------------------------------------------------------------
#    # cross-validation
#    
#    library(multicore,lib="~/R/x86_64-pc-linux-gnu-library/2.12/")
#    options(cores=10)
#    multicore=TRUE
#    
#    #load("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/sim_MeanTheta_startingvalue_015.Rdata")
#    #s_values <- c(0.01,0.02,0.03,0.04,0.05,seq(0.08,0.3,by=0.03))
#    #s_values <- c(0.02,0.05,0.07,seq(0.1,0.3,by=0.1))
#    #lambda_values <- c(1,10,100,1000)
#    #s_values <- c(0.05,0.1,0.15,0.2,0.25,0.3)
#    
#    lambda_values <- c(0.01,0.1,1,10)
#    s_values <- 5
#    
#    K <- 5  
#    n <- length(time)
#    nK <- floor(n/K)
#    neworder <- sample(1:n,n)
#    
#    # Matrix test data sets
#    test_data <- c()
#    for (k in 1:K)
#    {
#        if (k < K)
#        {
#            ik <- neworder[(k-1)*nK + (1:nK)]
#        }
#        else
#        {
#            ik <- neworder[((K-1)*nK+1):n]
#        }
#        test_data[[k]] <- ik
#    }
#    
#    # use the same test_data as in...
#    load("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/cv_034.Rdata")
#    
#    cv <- matrix(0,length(s_values),length(lambda_values))
#    colnames(cv) <- lambda_values
#    rownames(cv) <- s_values
#    
#    j <- 1
#        for (lam in lambda_values)
#        {
#            i <- 1
#            for (sv in s_values)
#            {
#                for(k in 1:K)
#                {
#                    
#                    ik <- test_data[[k]]
#                    
#                    # fit with training set
#                    spatFit <- spatialFit(xyblack_list, xywhite_list, I=I, J=J, x=scaledX[-ik,], 
#                                  Response=Response[,-ik], time=time[-ik], Theta=MeanTheta, 
#                                  lambda=lam, sv=sv, firstpen=TRUE, centersOnly=TRUE)
#                                                                  
#                    # prediction error on test data set     
#                    N <- I*J
#                    fit <- array(NA,c(N,length(ik)))
#                    for(n in 1:N)
#                    {
#                        fit[n,] <-  spatFit$reducedTheta[n,spatFit$selected[[n]]]%*%t(scaledX[ik,spatFit$selected[[n]]])
#                    
#                    }
#                                             
#                    prediction_error <- sum((Response[,ik] - fit[,])^2,na.rm=TRUE)
#    
#                                                              
#                    cv[i,j] <- cv[i,j] + prediction_error 
#                                   
#                }   
#                #name <- paste("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatFit_030_k",k,"_i",i,"_j",j,".Rdata", sep="")
#                #save(ik, lam, sv, spatFit,prediction_error,file=name)  
#                
#                print("cv for")
#                print(c(lam,sv))
#                print(cv[i,j])         
#                
#                i <- i+1
#            }
#            j <- j+1
#        }
#    
#    name <- paste("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/cv_scan", var_s, sep="")
#    name <- paste(name,"_11.Rdata",sep="")
#    
#    save(test_data,cv,file=name)
#
#}


