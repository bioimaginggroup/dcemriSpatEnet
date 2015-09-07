#setwd("C://Users//Julia//Desktop//Homeoffice_201203//Jan-Gertheis//Code_spatEnet//")
#setwd("C://Users//Julia//Desktop//Homeoffice_201203//Jan-Gertheis//")
setwd("Z://Dateien//Projekte//Jan-Gertheis//Code_spatEnet//")      ## Windows
setwd("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/")            ## Server

## breast cancer data
library(oro.nifti)
s <- 3
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
	
	
	I <- nrow(conc)
	J <- ncol(conc)
	N <- I*J        # number of pixels

	

#################################
# simulated data
load("~/Dateien/Software/R-Skripte/Output/SpatialDCEMRI/save_simulation_spatial_blocks_sigma2_09.Rdata")

load("save_simulation_spatial_blocks_withoutnoise.Rdata")
load("Z://Dateien//Software//R-Skripte//Output//SpatialDCEMRI//save_simulation_spatial_blocks_sigma2_09.Rdata")

load("save_simulation_spatial_blocks_withoutnoise.Rdata")
I <- 25 # numer of rows
J <- 25 # numer of columns

load("save_simulation_spatial_blocks_withoutnoise_100.Rdata")
I <- 100 #numer of rows
J <- 100 # numer of columns

#I <- 25 # numer of rows
#J <- 25 # numer of columns

#I <- 10 # number of pixels in x-direction
#J <- 8 # number of pixels in y-direction
#I <- 30 # number of pixels in x-direction
#J <- 30 # number of pixels in y-direction

N <- I*J # number of pixels


# data
#load("Z:/Research/DCE-MRI/Response.RData")
#load("Z://Dateien//Projekte//Jan-Gertheis//Response.RData")
#load("Response.RData")
#time <- seq(0,9.1,length=46)


#########################################
# functions
#source("Z:/Research/DCE-MRI/onePixFunctions.r")
source("onePixFunctions.r")
source("spatialFunctions.r")

# Xmatrix
Matrix <- Cp(time, model="tofts.kermode") # First curve: plasma concentration
#Matrix <- c()     # without Cp
# basis function for different kep values
#for (kep in exp(seq(-2,2,1)))
for (kep in exp(seq(-5,3,0.05)))
{    
    ft <- model.weinmann(time,kep,model="tofts.kermode")
    Matrix <- cbind(Matrix,ft)
}

# Define dimensions
T <- length(time)
Q <- ncol(Matrix) # number of basis functions, number of parameters



# True, simulated Theta-Matrix
Theta_true <- array(0, c(I*J,Q))   # 2D array of Parameter values per Pixel
Theta_Q <- exp(seq(-5,3,0.05))

for(i in 1:I)
{
  for(j in 1:J) 
  {
      for(k in 1:(Q-1))
      {
          #print(print(c(i,j,k)))
          
          if((Theta_Q[k] < PAR_KEP1[i,j]) && (PAR_KEP1[i,j] <= Theta_Q[k+1]))
          {
              #print("True 1")
              Theta_true[(i-1)*J+j,k] <- PAR_KTRANS1[i,j]
          }
          
          if((Theta_Q[k] < PAR_KEP2[i,j]) && (PAR_KEP2[i,j] <= Theta_Q[k+1]))
          {
              #print("True 2")
              Theta_true[(i-1)*J+j,k] <- PAR_KTRANS2[i,j]
          }   
          if((Theta_Q[k] < PAR_KEP3[i,j]) && (PAR_KEP3[i,j] <= Theta_Q[k+1]))
          {
              #print("True 3")
              Theta_true[(i-1)*J+j,k] <- PAR_KTRANS3[i,j]
          }         
      }        
  }
}

save(Theta_true,file="sim_Theta_true_123Comp.Rdata")

# true
z_range <- c(-0.0001,2.1)
name <- "sim_map_Ktrans_true_123Comp.pdf"
pdf(name, width = 40, height = 40)
par(mfrow=c(10,10), mar=c(1, 1, 1, 1), cex=1)
for(k in 51:150)
{
    Thetak_matrix <- t(matrix(Theta_true[,k], nrow=I))  
    image(Thetak_matrix, axes=FALSE,zlim=z_range,col=gray(100:0/100),main=paste("k=",k))
}
dev.off()

load(file="sim_Theta_true.Rdata")

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

# as lists
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
  
    #print((i-1)*J+j)
    #Response[(i-1)*J+j,] <-  as.vector(MRIdata[(MRIdata$x==i)&(MRIdata$y==j),1])
    #Response[(i-1)*J+j,] <- CONC_SIM_NOISE[i,j,]
    Response[(i-1)*J+j,] <- conc[i,j,]
        
  }
}
 
# Mean Response on ROI level 
summed_curve <- array(0, T)
for(i in 1:I)
{
  for(j in 1:J) 
  {
     #summed_curve <- summed_curve +  as.vector(MRIdata[(MRIdata$x==i)&(MRIdata$y==j),1])
     summed_curve <- summed_curve +  CONC_SIM_NOISE[i,j,]
  }
}
MeanResponse <- array(NA, c(I*J,T))
for(i in 1:I)
{
  for(j in 1:J) 
  {
      MeanResponse[(i-1)*J+j,] <- summed_curve/N
  }
}

#######################
## Voxelwise Fit with MeanResponse
#######################

test <- nnenet(scaledX, MeanResponse[1,], lambda=0.0001, 0.125, firstpen=FALSE)

MeanTheta <- array(NA, c(I*J,Q))
MeanFit <- array(NA, c(I*J,T))
for (pix in 1:length(xysample[,1]))
{
      MeanTheta[pix,] <- test$beta
      MeanFit[pix,] <- test$fit      
}

save(Response,scaledX,Theta,Fit,file="scan03_MeanTheta_voxelwise_lambda0k0001_TK.Rdata")

      
#######################
## Voxelwise Fit 
#######################
Theta <- array(NA, c(I*J,Q))
Fit <- array(NA, c(I*J,T))
for (pix in 1:length(xysample[,1]))
{
      single_index <- xysample[pix,]

      test <- nnenet(scaledX, Response[indx(single_index,I,J),], lambda=10, 0.125, firstpen=FALSE)
      Theta[indx(single_index,I,J),] <- test$beta
      Fit[indx(single_index,I,J),] <- test$fit      
}

save(Response,scaledX,Theta,Fit,file="sim_Theta_voxelwise_lambda10_TK_withoutnoise.Rdata")

load("sim_Theta_voxelwise_lambda0k001_TK.Rdata")

#load("save_Theta_voxelwise_lambda0k1_TK.Rdata")
#load("save_Theta_voxelwise.RData")



#######################
## Mean coefficients (Alternative 2) 
#######################
MeanTheta <- array(NA, c(I*J,Q))
for(q in 1:Q)
{
      MeanTheta[1:N,q] <- mean(Theta[,q])
}

name <- "plotCoefficients_voxelwise_lambda10_sim_true_withoutnoise.pdf"
pdf(name, width = 16, height = 16)
par(mfrow=c(8,8), mar=c(4, 1, 1, 1), cex=1)
for(i in 1:8)
{
  for(j in 1:8)
  {
      #plot(seq(-2,2,by=0.01),(Theta2[indx(c(i,j),I,J),]/sdX)[-1],type="h",xlab=expression(log(k[ep])),ylab=expression(K^{trans}), lwd=1.8)
      #points(seq(-2,2,by=0.01),(Theta[indx(c(i,j),I,J),]/sdX)[-1],type="h", col="red")
      plot(seq(-5,3,0.05),(Theta[indx(c(i,j),I,J),]/sdX)[-1],type="h",xlab=expression(log(k[ep])),ylab=expression(K^{trans}), ylim=c(0,0.1), col="red")
      abline(v=c(log(PAR_KEP1[i,j]),log(PAR_KEP2[i,j])), col="blue")
      points(seq(-5,3,0.05),(Theta_true[indx(c(i,j),I,J),]/sdX)[-1],type="h", ylim=c(0,0.1), col="grey")
            
  }
}
dev.off()

name <- "plotFits1_voxelwise_lambda0k01_sim.pdf"
pdf(name, width = 16, height = 16)
par(mfrow=c(8,8), mar=c(4, 1, 1, 1), cex=1)
for(i in 1:8)
{
  for(j in 1:8)
  {
      plot(time, Response[indx(c(i,j),I,J),])
      points(time,Fit[indx(c(i,j),I,J),],col="black", type="l", lwd=1.8)
      #points(time,MeanFit[indx(c(i,j),I,J),], col="red", type="l")
  }
}
dev.off()


##################################
## Spatial fit: parallel update
#################################
library(multicore,lib="~/R/x86_64-pc-linux-gnu-library/2.12/")
#require(multicore)
options(cores=8)
multicore=TRUE

#Theta2 <- Theta
#Fit2 <- array(NA, c(I*J,T))

#Theta2 <- MeanTheta
#load("sim_Theta2_lambda0k01_TK_Theta_True.Rdata")
load("sim_Theta_voxelwise_lambda10_TK_withoutnoise.Rdata")
#Theta2 <- Theta
#Theta2 <- Theta_true
Theta2 <- MeanTheta

Fit2 <- array(NA, c(I*J,T))

S <- 50 # number of iterations
sv2 <- 0.75
#sv2 <- 0.125
lambda2 <- 100
                         
for (s in 1:S)
{
    # alle schwarzen Pixel parallel updaten
    test3 <- mclapply(xyblack_list, FUN = callSpatEnet, I=I, J=J, x=scaledX, Response=Response, 
                    Theta=Theta2, lambda=lambda2, sv=sv2,firstpen=FALSE)
    # aktualisiere Werte Theta2    
    for(i in 1:length(xyblack_list))
    {
         Theta2[indx(xyblack_list[[i]],I,J),] <- test3[[i]]$beta
         Fit2[indx(xyblack_list[[i]],I,J),] <- test3[[i]]$fit 
    }
    
    # alle weißen Pixel parallel updaten   
    test3 <- mclapply(xywhite_list, FUN = callSpatEnet, I=I, J=J, x=scaledX, Response=Response, 
                    Theta=Theta2, lambda=lambda2, sv=sv2,firstpen=FALSE)
                    
    # aktualisiere Werte Theta2    
    for(i in 1:length(xywhite_list))
    {
         Theta2[indx(xywhite_list[[i]],I,J),] <- test3[[i]]$beta
         Fit2[indx(xywhite_list[[i]],I,J),] <- test3[[i]]$fit 
    }

}  

save(Response,scaledX,Theta, Theta2, Fit, Fit2, file="sim_Theta2_lambda100_TK_Theta_mean_withoutnoise.Rdata")

#### simulated annealing
load("sim_Theta2_lambda100_TK_Theta_mean_withoutnoise.Rdata")

Theta2 <- MeanTheta
Fit2 <- array(NA, c(I*J,T))

S <- 50 # number of iterations
sv2 <- 0.75

lambda2 <- 1000
for (s in 1:5)
{
    lambda2 <- lambda2/10 
    for (r in 1:S)
    {
    
        # alle schwarzen Pixel parallel updaten
        test3 <- mclapply(xyblack_list, FUN = callSpatEnet, I=I, J=J, x=scaledX, Response=Response, 
                        Theta=Theta2, lambda=lambda2, sv=sv2,firstpen=FALSE)
        # aktualisiere Werte Theta2    
        for(i in 1:length(xyblack_list))
        {
             Theta2[indx(xyblack_list[[i]],I,J),] <- test3[[i]]$beta
             Fit2[indx(xyblack_list[[i]],I,J),] <- test3[[i]]$fit 
        }
        
        # alle weißen Pixel parallel updaten   
        test3 <- mclapply(xywhite_list, FUN = callSpatEnet, I=I, J=J, x=scaledX, Response=Response, 
                        Theta=Theta2, lambda=lambda2, sv=sv2,firstpen=FALSE)
                        
        # aktualisiere Werte Theta2    
        for(i in 1:length(xywhite_list))
        {
             Theta2[indx(xywhite_list[[i]],I,J),] <- test3[[i]]$beta
             Fit2[indx(xywhite_list[[i]],I,J),] <- test3[[i]]$fit 
        }    
    }
        
    name <- paste("scan03_map_Ktrans_lambda", lambda2, "_Theta_mean.Rdata", sep="")
    save(Response,scaledX,Theta, Theta2, Fit, Fit2, file=name)
    
    z_range <- c(-0.001,0.05)

    name <- paste("scan03_map_Ktrans_lambda", lambda2,"_Theta_mean.pdf", sep="")
    pdf(name, width = 40, height = 40)
    par(mfrow=c(10,10), mar=c(1, 1, 1, 1), cex=1)
    for(k in 51:150)
    {
        Thetak_matrix <- t(matrix(Theta2[,k], nrow=J))  
        image(Thetak_matrix, axes=FALSE,zlim=z_range,col=gray(0:50/50),main=paste("k=",k))
    }
    dev.off()
    
}  

# only plot maps
lambda2 <- 1000
for (s in 1:5)
{
    lambda2 <- lambda2/10 
    
    name <- paste("scan03_map_Ktrans_lambda", lambda2, "_Theta_mean.Rdata", sep="")
    load(name)
    
    z_range <- c(-0.001,0.05)

    name <- paste("scan03_map_Ktrans_lambda", lambda2,"_Theta_mean.pdf", sep="")
    pdf(name, width = 40, height = 40)
    par(mfrow=c(10,10), mar=c(1, 1, 1, 1), cex=1)
    for(k in 51:150)
    {
        Thetak_matrix <- t(matrix(Theta2[,k], nrow=J))  
        image(Thetak_matrix, axes=FALSE,zlim=z_range,col=gray(0:50/50),main=paste("k=",k))
    }
    dev.off()

}


######################
##  Visualize results
#####################
load("save_simulation_123Comp.Rdata")
name <- "sim_123Comp_True.pdf"
pdf(name, width = 12, height = 4)
par(mfrow=c(1,3), mar=c(4, 1, 1, 1), cex=1)
for(i in c(1,26,51))
{
  for(j in 1:1)
  {
      #plot(time, Response[indx(c(i,j),I,J),], ylim=c(-0.05,0.3))
      plot(time, Response[indx(c(i,j),I,J),])
      #points(time,Fit2[indx(c(i,j),I,J),],col="black", type="l", lwd=1.8)
      #points(time,MeanFit[indx(c(i,j),I,J),], col="red", type="l")
  }
}
dev.off()


name <- "sim_plotCoefficients1_123Comp_True.pdf"
pdf(name, width = 48, height = 48)
par(mfrow=c(24,24), mar=c(4, 1, 1, 1), cex=1)
for(i in 1:24)
{
  for(j in 1:24)
  {
      #plot(seq(-2,2,by=0.01),(Theta2[indx(c(i,j),I,J),]/sdX)[-1],type="h",xlab=expression(log(k[ep])),ylab=expression(K^{trans}), lwd=1.8)
      #points(seq(-2,2,by=0.01),(Theta[indx(c(i,j),I,J),]/sdX)[-1],type="h", col="red")
      plot(seq(-5,3,0.05),(Theta2[indx(c(i,j),I,J),]/sdX)[-1],type="h",xlab=expression(log(k[ep])),ylab=expression(K^{trans}), ylim=c(0,0.07))
      #points(seq(-5,3,0.05),(MeanTheta[indx(c(i,j),I,J),]/sdX)[-1],type="h", col="grey")
      points(seq(-5,3,0.05),(Theta_true[indx(c(i,j),I,J),]/sdX)[-1],type="h", col="grey")
      points(seq(-5,3,0.05),(Theta[indx(c(i,j),I,J),]/sdX)[-1],type="h", col="red")
      points(seq(-5,3,0.05),(Theta2[indx(c(i,j),I,J),]/sdX)[-1],type="h", col="black")
                  
  }
}
dev.off()

# Parameter map

# Theta2 als Matrix
#k <- 85
#Thetak_matrix <- t(matrix(Theta2[,k], nrow=I))  

#z_range <- range(Theta2)

z_range <- c(-0.001,0.05)
name <- "sim_map_Ktrans_lambda100_Theta_mean_withoutnoise_100.pdf"
pdf(name, width = 40, height = 40)
par(mfrow=c(10,10), mar=c(1, 1, 1, 1), cex=1)
for(k in 51:150)
{
    Thetak_matrix <- t(matrix(Theta2[,k], nrow=I))  
    image(Thetak_matrix, axes=FALSE,zlim=z_range,col=gray(0:50/50),main=paste("k=",k))
}
dev.off()



# voxelwise
#save(Response,scaledX,Theta,Fit,file="sim_Theta_voxelwise_lambda0k0001.Rdata")
load("sim_Theta_voxelwise_lambda0k0001.Rdata")
name <- "sim_map_Ktrans_voxelwise_lambda0k0001.pdf"
pdf(name, width = 40, height = 40)
par(mfrow=c(10,10), mar=c(1, 1, 1, 1), cex=1)
for(k in 51:150)
{
    Thetak_matrix <- t(matrix(Theta[,k], nrow=I))  
    image(Thetak_matrix, axes=FALSE,zlim=z_range,col=gray(50:0/50),main=paste("k=",k))
}
dev.off()

sse <- array(NA, I*J)
error <- Response - Fit
for(n in 1:N)
{
	 sse[n] <- sum(error[n,]^2)
}

name <- paste("sim_sse_voxelwise_lambda_lambda0k0001", lambda2,"_002.pdf", sep="")
pdf(name, width = 4, height = 4)
par(mfrow=c(1,1), mar=c(1, 1, 1, 1), cex=1)

sse_matrix <- t(matrix(sse, nrow=J))
image.plot(sse_matrix, axes=FALSE,col=gray(0:50/50),main=sse)

dev.off()


# meanTheta
name <- "sim_map_lambda0k0001_MeanTheta_fromMeanConc_100.pdf"
pdf(name, width = 40, height = 40)
par(mfrow=c(10,10), mar=c(1, 1, 1, 1), cex=1)
for(k in 51:150)
{
    Thetak_matrix <- t(matrix(MeanTheta[,k], nrow=I))  
    image(Thetak_matrix, axes=FALSE,zlim=z_range,col=gray(0:50/50),main=paste("k=",k))
}
dev.off()

# true
name <- "sim_map_Ktrans_true_100.pdf"
pdf(name, width = 40, height = 40)
par(mfrow=c(10,10), mar=c(1, 1, 1, 1), cex=1)
for(k in 51:150)
{
    Thetak_matrix <- t(matrix(Theta_true[,k], nrow=I))  
    image(Thetak_matrix, axes=FALSE,zlim=z_range,col=gray(0:50/50),main=paste("k=",k))
}
dev.off()

###    true without noise
name <- "sim_map_Ktrans_true_withoutnoise.pdf"
pdf(name, width = 40, height = 40)
par(mfrow=c(10,10), mar=c(1, 1, 1, 1), cex=1)
for(k in 51:150)
{
    Thetak_matrix <- t(matrix(Theta_true[,k], nrow=I))  
    image(Thetak_matrix, axes=FALSE,zlim=z_range,col=gray(0:50/50),main=paste("k=",k))
}

dev.off()



#postscript("Z:/Research/DCE-MRI/Paper/plotVox1.ps",height=6,width=6,
#horizontal=F)
name <- "plotBasisFunctions_Toftskermode.pdf"
pdf(name, width = 8, height = 8)
plot(Matrix[,1], type="l", ylim=c(0,5.3), xlab="time point", ylab="response / predictors")
for (i in 2:ncol(Matrix))
lines(Matrix[,i])
lines(Response,col="blue",lwd=2)
#title(paste("x =", xp, ", y =", yp))
dev.off()



#Theta0 <- Theta
#Fit0 <- Fit
#
## Werte nahe Null auf Null setzen
#for(n in 1:N)
#{
#    for(q in 1:Q)
#    {
#       if(abs(Theta0[n,q])< 10^(-5))    
#       {
#          Theta0[n,q] <- 0
#          
#          
#       }   
#       
#    }
#    
#    Fit0[n,] <- scaledX %*% Theta0[n,]
#}
#

#test <- lapply(xysample, FUN= callNnenet, ...)

plot(seq(-2,2,by=0.01),Theta[indx(c(5,5),I,J),][-1],type="h",xlab=expression(log(k[ep])),
ylab=expression(K^{trans}))

plot(time, Response[indx(c(5,5),I,J),], xlab="observed concentrations",ylab="fitted concentrations")
points(time,Fit[indx(c(5,5),I,J),], col="red")
 
#plot(seq(-2,2,by=0.01),(nnenet1$beta/sdX)[-1],type="h",xlab=expression(log(k[ep])),
#ylab=expression(K^{trans}))

# räumlicher Fit nur für einen Pixel -> für lambda=1 wird die Kurve überschätzt, lambda=0.01 passt
single_index <- c(5,5)
Theta2 <- Theta
Fit2 <- array(NA, c(I*J,T))

test2 <- callSpatEnet(single_index, I, J, scaledX, Response, Theta2, lambda=1, sv=0.125, firstpen=FALSE)
      Theta2[indx(single_index,I,J),] <- test2$beta
      Fit2[indx(single_index,I,J),] <- test2$fit 

plot(time, Response[indx(c(5,5),I,J),], xlab="observed concentrations",ylab="fitted concentrations")
points(time,Fit[indx(c(5,5),I,J),], col="red")
points(time,Fit2[indx(c(5,5),I,J),], col="blue")

plot(seq(-2,2,by=1),Theta2[indx(c(5,5),I,J),][-1],type="h",xlab=expression(log(k[ep])),
ylab=expression(K^{trans}))
points(seq(-2,2,by=1),Theta[indx(c(5,5),I,J),][-1],type="h",col="red")

##  Tests für einzelnen Pixel
x <- scaledX
lambda <- 10^(-5)
sv <- 0.125

y <- Response[indx(single_index,I,J),]
neighborList <- getNeighbors(single_index,I,J)
Xi <- getXi(Theta, neighborList, I, J)
D <- getD(ncol(x),neighborList)

test4 <- spatEnet(x, y, lambda, sv, firstpen=T, Xi, D)


# Test mit lambda=1
# läuft ewig, selbst für kleinen Bildausschnitt
Theta2 <- Theta
Fit2 <- array(NA, c(I*J,T))
for (s in 1:10)
{
  for(pix in 1:length(xyblack[,1]))
  {
      single_index <- xyblack[pix,]
      
      test2 <- callSpatEnet(single_index, I, J, scaledX, Response, Theta2, lambda=10^(-5), sv=0.125, firstpen=FALSE)
      Theta2[indx(single_index,I,J),] <- test2$beta
      Fit2[indx(single_index,I,J),] <- test2$fit    
  }
  for(pix in 1:length(xywhite[,1]))
  {
      single_index <- xywhite[pix,]
      
      test2 <- callSpatEnet(single_index, I, J, scaledX, Response, 
                            Theta2, lambda=10^(-5), sv=0.125, firstpen=FALSE)
      Theta2[indx(single_index,I,J),] <- test2$beta
      Fit2[indx(single_index,I,J),] <- test2$fit    
  }
} 

save(Response,scaledX,Theta, Theta2, Fit, Fit2, file="save_Theta_mitCpFPF_lambda10hm5.Rdata")




     
                              
################################################################################
# Ergebnisse grafisch
theta_12 <- array(Theta[,2], c(I,J))
theta_22 <- array(Theta2[,2], c(I,J))
image(theta_12)
X11()
image(theta_22)
# -> Ergebnisse fast identisch mit nicht-räumlich für  lambda=10^(-5), 10^(-1)
# -> Ergebnisse glatter für lambda=100, 10,

name <- "plotFit_55.pdf"
pdf(name, width = 5, height = 5)
plot(time, Response[55,], xlab="observed concentrations",ylab="fitted concentrations")
points(time,Fit2[55,])
points(time,Fit[55,], col="red")
dev.off()



name <- "plotCoefficients_46.pdf"
index <- xysample[46,]
neighbor_list <- getNeighbors(index, I,J)

pdf(name, width = 8, height = 8)
par(mfrow=c(3,3), mar=c(4, 1, 1, 1), cex=1)

plot(1, type="n", axes=F, xlab="", ylab="")

# oben Mitte
plot(seq(-2,2,by=0.01),(Theta2[indx(neighbor_list[[3]],I,J),]/sdX)[-1],type="h",xlab=expression(log(k[ep])),ylab=expression(K^{trans}))
points(seq(-2,2,by=0.01),(Theta[indx(neighbor_list[[3]],I,J),]/sdX)[-1],type="h", col="red")
plot(1, type="n", axes=F, xlab="", ylab="")

# Mitte links
plot(seq(-2,2,by=0.01),(Theta2[indx(neighbor_list[[1]],I,J),]/sdX)[-1],type="h",xlab=expression(log(k[ep])),ylab=expression(K^{trans}))
points(seq(-2,2,by=0.01),(Theta[indx(neighbor_list[[1]],I,J),]/sdX)[-1],type="h", col="red")
# Mitte mitte
plot(seq(-2,2,by=0.01),(Theta2[indx(index,I,J),]/sdX)[-1],type="h",xlab=expression(log(k[ep])),ylab=expression(K^{trans}))
points(seq(-2,2,by=0.01),(Theta[indx(index,I,J),]/sdX)[-1],type="h", col="red")
# Mitte rechts
plot(seq(-2,2,by=0.01),(Theta2[indx(neighbor_list[[2]],I,J),]/sdX)[-1],type="h",xlab=expression(log(k[ep])),ylab=expression(K^{trans}))
points(seq(-2,2,by=0.01),(Theta[indx(neighbor_list[[2]],I,J),]/sdX)[-1],type="h", col="red")

plot(1, type="n", axes=F, xlab="", ylab="")
# unten Mitte
plot(seq(-2,2,by=0.01),(Theta2[indx(neighbor_list[[4]],I,J),]/sdX)[-1],type="h",xlab=expression(log(k[ep])),ylab=expression(K^{trans}))
points(seq(-2,2,by=0.01),(Theta[indx(neighbor_list[[4]],I,J),]/sdX)[-1],type="h", col="red")
plot(1, type="n", axes=F, xlab="", ylab="")

dev.off()




#######################
## Tests konstruiert
#######################
# Starting values for Theta
Theta <- array(runif(I*J*Q,0.01,0.2), c(I*J,Q))

single_index <- xysample[71,]

# zum Vergleich: nicht-räumlicher Fit
test1 <- nnenet(scaledX, Response[indx(single_index,I,J),], 10^(-5), 0.125, firstpen=T) 

# räumlicher Fit
test2 <- callSpatEnet(single_index, I, J, scaledX, Response, Theta, lambda=10^(-5), sv=0.125, firstpen=TRUE)

plot(time, Response[indx(single_index,I,J),])
points(time, test1$fit, col="blue")
points(time, test2$fit, col="red")


## Belegung zum Test spatEnet
x <- scaledX
lambda <- 10^(-5)
sv <- 0.125

y <- Response[indx(single_index,I,J),]
neighborList <- getNeighbors(single_index,I,J)
Xi <- getXi(Theta, neighborList, I, J)
D <- getD(ncol(x),neighborList)

test4 <- spatEnet(x, y, lambda, sv, firstpen=T, Xi, D)
# liefert genau das gleiche Ergebnis wie Aufruf von callSpatEnet



#### Beispiel konstruiert
single_index <- xysample[71,]
y <- Response[indx(single_index,I,J),]
Xi <- c(0,0.03,0.5,0,0,0, 0,0.03,0.5,0,0,0, 0,0.03,0.5,0,0,0)
E <- diag(1,nrow=6)
D <- rbind(E,E,E)

test5 <- spatEnet(x, y, lambda, sv, firstpen=T, Xi, D)

test6 <- spatEnet(x, y, lambda=0.01, sv, firstpen=T, Xi, D)



