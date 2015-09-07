setwd("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/")            ## Server

# simulated data
#load("save_simulation_spatial_blocks_withoutnoise.Rdata")
#I <- 25 # numer of rows
#J <- 25 # numer of columns
load("save_simulation_123vpComp.Rdata")
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
#for (kep in exp(seq(-2,2,1)))

#for (kep in exp(seq(-5,3,0.05)))
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
    Response[(i-1)*J+j,] <- CONC_SIM_NOISE[i,j,]
  }
}

# ------------------------------------------------------------------------------
# starting value for spatial Fit
# Mean Response on ROI level
summed_curve <- array(0, T)
for(i in 1:I)
{
  for(j in 1:J)
  {
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

test <- nnenet(scaledX, MeanResponse[1,], lambda=0.0001, 0.125, firstpen=FALSE)

MeanTheta <- array(NA, c(I*J,Q))
MeanFit <- array(NA, c(I*J,T))
for (pix in 1:length(xysample[,1]))
{
      MeanTheta[pix,] <- test$beta
      MeanFit[pix,] <- test$fit
}

#save(MeanTheta,file="~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/sim_MeanTheta_startingvalue_303.Rdata")




## ------------------------------------------------------------------------------
## Voxelwise Fit
#Theta <- array(NA, c(I*J,Q))
#Fit <- array(NA, c(I*J,T))
#Theta2 <- array(0, c(I*J,Q))
#Fit2 <- array(NA, c(I*J,T)) 
#Selected_array <- vector("list", N)
#
#for (pix in 1:length(xysample[,1]))
#{
#      single_index <- xysample[pix,]
#
#      #test <- nnenet(scaledX, Response[indx(single_index,I,J),], lambda=10^(-1), 0.75, firstpen=TRUE)
#      
#      #print(single_index)
#      
#      test <- nnenet2(scaledX, Response[indx(single_index,I,J),], lambda=10^(-12), 1, firstpen=TRUE)
#      
#      Theta[indx(single_index,I,J),] <- test$beta
#      Fit[indx(single_index,I,J),] <- test$fit
#      Theta2[indx(single_index,I,J),test$selected] <- test$reducedBeta
#      Fit2[indx(single_index,I,J),] <- test$reducedFit
#      Selected_array[[indx(single_index,I,J)]] <- test$selected
#}
#
#save(Response,scaledX,Theta,Fit,Theta2,Fit2,file="~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/voxelwise/sim_Theta_voxelwise_023.Rdata")
#
#

 ------------------------------------------------------------------------------
# spatial fit for fixed lambda and sv

#library(multicore,lib="~/R/x86_64-pc-linux-gnu-library/2.12/")
#library(multicore)
require(multicore)
options(cores=8)
multicore=TRUE

#load("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/sim_MeanTheta_startingvalue_015.Rdata")
#load("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/sim_MeanTheta_startingvalue_202.Rdata")

sv2 <- 0.5
lambda2 <- 200
spatFit <- spatialFit(xyblack_list, xywhite_list, I=I, J=J, x=scaledX, Response=Response, time=time,
            Theta=MeanTheta, lambda=lambda2, sv=sv2, firstpen=TRUE, centersOnly=TRUE)

save(Response,scaledX,time,spatFit,sv2,lambda2,file="~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/spatFit_367.Rdata")
  

### ------------------------------------------------------------------------------
## cross-validation
#
#library(multicore,lib="~/R/x86_64-pc-linux-gnu-library/2.12/")
#options(cores=6)
#multicore=TRUE
#
##load("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/sim_MeanTheta_startingvalue_015.Rdata")
##s_values <- c(0.01,0.02,0.03,0.04,0.05,seq(0.08,0.3,by=0.03))
##s_values <- c(0.02,0.05,0.07,seq(0.1,0.3,by=0.1))
##lambda_values <- c(1,10,100,1000)
##s_values <- c(0.05,0.1,0.15,0.2,0.25,0.3)
#
#lambda_values <- c(500,1000,1500)
#s_values <- 0.2
#
#K <- 5  
#n <- length(time)
#nK <- floor(n/K)
#neworder <- sample(1:n,n)
#
#cv <- matrix(0,length(s_values),length(lambda_values))
#colnames(cv) <- lambda_values
#rownames(cv) <- s_values
#
## Matrix test data sets
#test_data <- c()
#for (k in 1:K)
#{
#    if (k < K)
#    {
#        ik <- neworder[(k-1)*nK + (1:nK)]
#    }
#    else
#    {
#        ik <- neworder[((K-1)*nK+1):n]
#    }
#    test_data[[k]] <- ik
#}
#
## use the same test_data as in...
#load("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/cv_034.Rdata")
#
#j <- 1
#    for (lam in lambda_values)
#    {
#        i <- 1
#        for (sv in s_values)
#        {
#            for(k in 1:K)
#            {
#                
#                ik <- test_data[[k]]
#                
#                # fit with training set
#                spatFit <- spatialFit(xyblack_list, xywhite_list, I=I, J=J, x=scaledX[-ik,], 
#                              Response=Response[,-ik], time=time[-ik], Theta=MeanTheta, 
#                              lambda=lam, sv=sv, firstpen=TRUE, centersOnly=TRUE)
#                              
#                # prediction error on test data set     
#                N <- I*J
#                fit <- array(NA,c(N,length(ik)))
#                for(n in 1:N)
#                {
#                    fit[n,] <-  spatFit$reducedTheta[n,spatFit$selected[[n]]]%*%t(scaledX[ik,spatFit$selected[[n]]])
#                
#                }
#                                         
#                prediction_error <- sum((Response[,ik] - fit[,])^2,na.rm=TRUE)
#
#                                                          
#                cv[i,j] <- cv[i,j] + prediction_error 
#                               
#            }   
#            #name <- paste("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatFit_030_k",k,"_i",i,"_j",j,".Rdata", sep="")
#            #save(ik, lam, sv, spatFit,prediction_error,file=name)  
#            
#            print("cv for")
#            print(c(lam,sv))
#            print(cv[i,j])         
#            
#            i <- i+1
#        }
#        j <- j+1
#    }
#
#save(test_data,cv,file="~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/cv_036.Rdata")
#


### ------------------------------------------------------------------------------
## cross-validation  voxelwise fit
#
#library(multicore,lib="~/R/x86_64-pc-linux-gnu-library/2.12/")
#options(cores=6)
#multicore=TRUE
#
#lambda_values <- c(10^(-1),10^(-5),10^(-10))
#s_values <- c(0.1,0.2,0.5)
#
#K <- 5  
#n <- length(time)
#nK <- floor(n/K)
#neworder <- sample(1:n,n)
#
#cv <- matrix(0,length(s_values),length(lambda_values))
#colnames(cv) <- lambda_values
#rownames(cv) <- s_values
#
## Matrix test data sets
#test_data <- c()
#for (k in 1:K)
#{
#    if (k < K)
#    {
#        ik <- neworder[(k-1)*nK + (1:nK)]
#    }
#    else
#    {
#        ik <- neworder[((K-1)*nK+1):n]
#    }
#    test_data[[k]] <- ik
#}
#
## use the same test_data as in...
#load("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/cv_034.Rdata")
#
#j <- 1
#    for (lam in lambda_values)
#    {
#        i <- 1
#        for (sv in s_values)
#        {
#            for(k in 1:K)
#            {
#                
#                ik <- test_data[[k]]
#                                
#                # fit with training set
#                Theta <- Theta2 <- array(NA, c(I*J,Q))
#                Fit <- Fit2 <- array(NA, c(I*J,(T-length(ik))))
#                Selected_array <- vector("list", N)
#                for (pix in 1:length(xysample[,1]))
#                {
#                    single_index <- xysample[pix,]
#                
#                    test <- nnenet2(x=scaledX[-ik,], Response[indx(single_index,I,J),-ik], lambda=lam, sv=sv, firstpen=TRUE)
#                    Theta[indx(single_index,I,J),] <- test$beta
#                    Fit[indx(single_index,I,J),] <- test$fit
#                    Theta2[indx(single_index,I,J),test$selected] <- test$reducedBeta
#                    Fit2[indx(single_index,I,J),] <- test$reducedFit
#                    Selected_array[[indx(single_index,I,J)]] <- test$selected
#                }
#                              
#                # prediction error on test data set     
#                N <- I*J
#                fit <- array(NA,c(N,length(ik)))
#                for(n in 1:N)
#                {
#                    fit[n,] <-  Theta2[n,Selected_array[[n]]]%*%t(scaledX[ik,Selected_array[[n]]])
#                
#                }
#                                         
#                prediction_error <- sum((Response[,ik] - fit[,])^2,na.rm=TRUE)
#
#                                                          
#                cv[i,j] <- cv[i,j] + prediction_error 
#                               
#            }   
#            #name <- paste("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatFit_030_k",k,"_i",i,"_j",j,".Rdata", sep="")
#            #save(ik, lam, sv, spatFit,prediction_error,file=name)  
#            
#            print("cv for")
#            print(c(lam,sv))
#            print(cv[i,j])         
#            
#            i <- i+1
#        }
#        j <- j+1
#    }
#
#print("end of cv")
#print(c(i,j))
#
#save(test_data,cv,file="~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/voxelwise_cv_02.Rdata")
#
#