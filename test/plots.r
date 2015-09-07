setwd("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/")            ## Server

# simulated data
#load("save_simulation_spatial_blocks_withoutnoise.Rdata")
#I <- 25 # numer of rows
#J <- 25 # numer of columns
#load("save_simulation_123Comp.Rdata")
load("save_simulation_123vpComp.Rdata")
load("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/spatFit_365.Rdata")
#load("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/scan03_Theta_voxelwise_lambda0k01_sv0k75.Rdata")

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


# ------------------------------------------------------------------------------
# plot  maps
require("fields")

z_range <- c(-0.001,0.05)
name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/Theta_map_365.pdf"
pdf(name, width = 32, height = 32)
par(mfrow=c(8,8), mar=c(1, 1, 1, 1), cex=1)
#for(k in c(1:3,55:150))
for(k in 1:Q)
{
    Thetak_matrix <- t(matrix(spatFit$SpatialTheta[,k], nrow=J))
    image(Thetak_matrix, axes=FALSE,zlim=z_range,col=gray(50:0/50),main=paste("k=",k))
}
dev.off()


z_range <- c(-0.001,0.2)
name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/reducedTheta_map_365.pdf"
pdf(name, width = 32, height = 32)
par(mfrow=c(8,8), mar=c(1, 1, 1, 1), cex=1)
for(k in 1:Q)
{
    Thetak_matrix <- t(matrix(spatFit$reducedTheta[,k], nrow=J))
    image(Thetak_matrix, axes=FALSE,zlim=z_range,col=gray(50:0/50),main=paste("k=",k))
}
dev.off()



require("fields")


i <- 20
j <- 20

#i <- 55
#j <- 30
name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/conc_scan13.pdf"
pdf(name, width = 4, height = 4)
par(mfrow=c(1,1), mar=c(4, 4, 4, 4), cex=1)

    #image.plot(conc[,,10], axes=FALSE,main="conc")
    image(1:I,1:J,conc[,,10], main="conc",xlab="i",ylab="j")
    abline(h=j)
    abline(v=i)
    points(i,j,cex=2)

dev.off()




#
##############################
## voxelwise
#z_range <- c(-0.001,0.05)
#name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/Theta_voxelwise_scan1.pdf"
#pdf(name, width = 32, height = 32)
#par(mfrow=c(8,8), mar=c(1, 1, 1, 1), cex=1)
##for(k in c(1:3,55:150))
#for(k in 1:Q)
#{
#    Thetak_matrix <- t(matrix(Theta[,k], nrow=J))
#    image(Thetak_matrix, axes=FALSE,zlim=z_range,col=gray(50:0/50),main=paste("k=",k))
#}
#dev.off()
#
#sse <- array(NA, I*J)
#error <- Response - Fit
#for(n in 1:N)
#{
#	 sse[n] <- sum(error[n,]^2)
#}
#
#name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/sse_voxelwise_scan1.pdf"
#pdf(name, width = 4, height = 4)
#par(mfrow=c(1,1), mar=c(1, 1, 1, 1), cex=1)
#
#    sse_matrix <- t(matrix(sse, nrow=J))
#    z_range <- c(0,0.15)
#    image.plot(sse_matrix, axes=FALSE,zlim=z_range,col=gray(0:50/50),main="sse")
#
#dev.off()


#z_range <- c(-0.001,0.05)
#name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/Theta2_map_365.pdf"
#pdf(name, width = 32, height = 32)
#par(mfrow=c(8,8), mar=c(1, 1, 1, 1), cex=1)
##for(k in c(1:3,55:150))
#for(k in 1:Q)
#{
#    Thetak_matrix <- t(matrix(Theta2[,k], nrow=J))
#    image(Thetak_matrix, axes=FALSE,zlim=z_range,col=gray(50:0/50),main=paste("k=",k))
#}
#dev.off()
#
#
#z_range <- c(-0.001,0.05)
#name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/selectedTheta_map_365.pdf"
#pdf(name, width = 32, height = 32)
#par(mfrow=c(8,8), mar=c(1, 1, 1, 1), cex=1)
##for(k in c(1:3,55:150))
#for(k in 1:Q)
#{
#    Thetak_matrix <- t(matrix(reducedTheta[,k], nrow=J))
#    image(Thetak_matrix, axes=FALSE,zlim=z_range,col=gray(50:0/50),main=paste("k=",k))
#}
#dev.off()


#z_range <- c(-0.001,0.5)
#name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/reducedTheta_map_365.pdf"
#pdf(name, width = 16, height = 8)
#par(mfrow=c(2,3), mar=c(1, 1, 1, 1), cex=1)
#for(k in 1:length(spatFit$reduced_k))
#{
#    Thetak_matrix <- t(matrix(spatFit$reducedTheta[,k], nrow=J))
#    image(Thetak_matrix, axes=FALSE,zlim=z_range,col=gray(50:0/50),main=paste("k=",spatFit$reduced_k[k]))
#}
#dev.off()




load("save_simulation_123vpComp.Rdata")
PAR_VP[which(is.na(PAR_VP)==TRUE)] <- 0

name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/Basis_Selection_365.pdf"
pdf(name, width = 20, height = 12)
par(mfrow=c(3,5), mar=c(4, 4, 1, 1), cex=1)
y_range <- c(-0.01,0.4)
for(i in c(5,35,70))
{
  for(j in c(5,20,40,50,70))
  {

      plot(1:Q, spatFit$SpatialTheta[indx(c(i,j),I,J),], ylab=expression(K^{trans}), type="h", ylim=y_range, main=paste("voxel",i,j))

      true_k <- c(1,c(15,38,45)+1)
      points(true_k, c(PAR_VP[i,j],PAR_KTRANS1[i,j],PAR_KTRANS2[i,j],PAR_KTRANS3[i,j]), type="h", col="green")

      points(1:Q,spatFit$reducedTheta[indx(c(i,j),I,J),]/sdX,col="red", type="h")
      points(1:Q, spatFit$SpatialTheta[indx(c(i,j),I,J),], type="h") # nochmal dr?ber

      #points(Selected_array[[indx(c(i,j),I,J)]],spatFit$SpatialTheta[indx(c(i,j),I,J),Selected_array[[indx(c(i,j),I,J)]]],col="red", type="h")
  }
}
dev.off()

#name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/Fit_365.pdf"
#pdf(name, width = 4, height = 4)
#par(mfrow=c(1,1), mar=c(4, 4, 1, 1), cex=1)
#y_range <- c(-0.01,1)
#n <- 2316
#      plot(time, spatFit$reducedFit[n,], ylab="t", type="l", ylim=y_range, main=paste("voxel",n))
#      points(time, Response[n,])
#
#      #points(Selected_array[[indx(c(i,j),I,J)]],spatFit$SpatialTheta[indx(c(i,j),I,J),Selected_array[[indx(c(i,j),I,J)]]],col="red", type="h")
#
#dev.off()
#
#name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/Selected_365.pdf"
#pdf(name, width = 4, height = 4)
#par(mfrow=c(1,1), mar=c(4, 4, 1, 1), cex=1)
#y_range <- c(-0.01,0.1)
#n <- 2316
#      plot(1:Q, spatFit$SpatialTheta[n,], ylab=expression(K^{trans}), type="h", ylim=y_range, main=paste("voxel",i,j))
#
#      points(1:Q,spatFit$reducedTheta[n,]/sdX,col="red", type="h")
#      points(1:Q, spatFit$SpatialTheta[n,], type="h") # nochmal dr?ber
#
#dev.off()
#


sse <- array(NA, I*J)
error <- Response - spatFit$SpatialFit
for(n in 1:N)
{
	 sse[n] <- sum(error[n,]^2)
}

name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/sse_Fit_365.pdf"
pdf(name, width = 4, height = 4)
par(mfrow=c(1,1), mar=c(1, 1, 1, 1), cex=1)

    sse_matrix <- t(matrix(sse, nrow=J))
    z_range <- c(0,0.15)
    image(sse_matrix, axes=FALSE,zlim=z_range,col=gray(0:50/50))

dev.off()


sse <- array(NA, I*J)
error <- Response - spatFit$reducedFit
for(n in 1:N)
{
	 sse[n] <- sum(error[n,]^2)
}



name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/sse_reducedFit_365.pdf"
pdf(name, width = 4, height = 4)
par(mfrow=c(1,1), mar=c(1, 1, 1, 1), cex=1)

    sse_matrix <- t(matrix(sse, nrow=J))
    #z_range <- c(0,0.15)
    z_range <- c(0,0.25)
    image(sse_matrix, axes=FALSE,zlim=z_range,col=gray(0:50/50))
    box()

dev.off()




# number of compartments
number <- array(NA, I*J)
vp_array <- array(NA, I*J)
for(n in 1:N)
{
	 #number[n] <- length(which(spatFit$selected[[n]]!=1))
	 number[n] <- length(which(which(spatFit$reducedTheta[n,]!=0)!=1))

	 vp_array[n] <- length(which(which(spatFit$reducedTheta[n,]!=0)==1))

	 #vp_array[n] <-  length(which(spatFit$selected[[n]]==1))
}

#name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/q_estimates_365.pdf"
#pdf(name, width = 4, height = 4)
#par(mfrow=c(1,1), mar=c(1, 1, 1, 1), cex=1)
#
#    number_matrix <- t(matrix(number, nrow=J))
#    z_range <- c(0,5)
#    image(number_matrix, axes=FALSE,zlim=z_range,col=heat.colors(n=5, alpha = 1))
#
#dev.off()
#
#name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/q_legend.pdf"
#pdf(name, width = 1, height = 4)
#par(mfrow=c(1,1), mar=c(0.1, 0.1, 0.1, 0.1),  cex=1.4)
#
#    image.plot(number_matrix, legend.only = TRUE, axes=FALSE,zlim=z_range,col=heat.colors(n=5, alpha = 1))
#
#dev.off()


name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/q_estimates_365.pdf"
pdf(name, width = 4, height = 4)
par(mfrow=c(1,1), mar=c(1, 1, 1, 1), cex=1)

    number_matrix <- t(matrix(number, nrow=J))
    z_range <- c(0,4)
    image(number_matrix, axes=FALSE,zlim=z_range,col=gray.colors(4, start=0.95, end=0))
    box()

dev.off()

name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/q_legend.pdf"
pdf(name, width = 1, height = 4)
par(mfrow=c(1,1), mar=c(0.1, 0.1, 0.1, 0.1),  cex=1.4)

    image.plot(number_matrix, legend.only = TRUE, axes=FALSE,axis.args = list(at=c(0.5,1.5,2.5,3.5), labels=c("1","2","3","4")),zlim=z_range,col=gray.colors(4, start=0.95, end=0))

dev.off()


name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/vp_estimates_365.pdf"
pdf(name, width = 4, height = 4)
par(mfrow=c(1,1), mar=c(1, 1, 1, 1), cex=1)

    vp_matrix <- t(matrix(vp_array, nrow=J))
    z_range <- c(0,1)
    image(vp_matrix, axes=FALSE,zlim=z_range,col=gray.colors(2, start=1, end=0.3))
    box()

dev.off()

name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/vp_legend.pdf"
pdf(name, width = 1, height = 4)
par(mfrow=c(1,1), mar=c(0.1, 0.1, 0.1, 0.1),  cex=1.4)

    image.plot(vp_matrix, legend.only = TRUE, axes=FALSE, axis.args = list(at=c(0.25,0.75), labels=c("=0",">0")),zlim=z_range,col=gray.colors(2, start=1, end=0.3))

dev.off()







# parameter maps
kep1 <- array(NA, I*J)
kep2 <- array(NA, I*J)
kep3 <- array(NA, I*J)
ktrans1 <- array(NA, I*J)
ktrans2 <- array(NA, I*J)
ktrans3 <- array(NA, I*J)

theta1 <- array(NA, I*J)
theta2 <- array(NA, I*J)
theta3 <- array(NA, I*J)

vp <- array(NA, I*J)

vp_comp <- which(spatFit$reducedTheta[,1]!=0)
vp[vp_comp] <-  spatFit$reducedTheta[vp_comp,1]

for(n in 1:N)
{

    first_comp <- which(spatFit$reducedTheta[n,-1]!=0)[1]
    #second_comp <- which(spatFit$reducedTheta[n,-1]!=0)[2]
    last_comp <- which(spatFit$reducedTheta[n,-1]!=0)[length(which(spatFit$reducedTheta[n,-1]!=0))]


    if(length(which(spatFit$reducedTheta[n,-1]!=0))>=3)
    {
        second_comp <- which(spatFit$reducedTheta[n,-1]!=0)[2]
    }

    #if(last_comp==first_comp||last_comp==second_comp){last_comp <- NA}

    if(last_comp==first_comp){last_comp <- NA}


    kep_values <- exp(seq(-3,3,0.1))
    kep1[n] <- kep_values[first_comp]

    ktrans1[n] <- spatFit$reducedTheta[n,first_comp+1]/sdX[first_comp+1]

    theta1[n] <- spatFit$reducedTheta[n,first_comp+1]


    if(length(which(spatFit$reducedTheta[n,-1]!=0))>=3)
    {
        kep2[n] <- kep_values[second_comp]
        ktrans2[n] <- spatFit$reducedTheta[n,second_comp+1]/sdX[second_comp+1]
        theta2[n] <- spatFit$reducedTheta[n,second_comp+1]

    }


    if(!is.na(last_comp))
    {
        kep3[n] <- kep_values[last_comp]
        ktrans3[n] <- spatFit$reducedTheta[n,last_comp+1]/sdX[last_comp+1]
        theta3[n] <- spatFit$reducedTheta[n,last_comp+1]
    }
}

name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/kep1_365.pdf"
pdf(name, width = 4, height = 4)
par(mfrow=c(1,1), mar=c(1, 1, 1, 1), cex=1)

    kep1_matrix <- t(matrix(kep1, nrow=J))
    z_range <- c(0,0.5)
    image.plot(kep1_matrix, axes=FALSE,zlim=z_range,col=gray(50:0/50),main="kep1")

dev.off()

name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/kep2_365.pdf"
pdf(name, width = 4, height = 4)
par(mfrow=c(1,1), mar=c(1, 1, 1, 1), cex=1)

    kep2_matrix <- t(matrix(kep2, nrow=J))
    z_range <- c(0,8)
    image.plot(kep2_matrix, axes=FALSE,zlim=z_range,col=gray(50:0/50),main="kep2")

dev.off()

name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/kep3_365.pdf"
pdf(name, width = 4, height = 4)
par(mfrow=c(1,1), mar=c(1, 1, 1, 1), cex=1)

    kep3_matrix <- t(matrix(kep3, nrow=J))
    z_range <- c(0,8)
    image.plot(kep3_matrix, axes=FALSE,zlim=z_range,col=gray(50:0/50),main="kep3")

dev.off()




name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/ktrans1_365.pdf"
pdf(name, width = 4, height = 4)
par(mfrow=c(1,1), mar=c(1, 1, 1, 1), cex=1)

    ktrans1_matrix <- t(matrix(ktrans1, nrow=J))
    z_range <- c(0,0.3)
    image.plot(ktrans1_matrix, axes=FALSE,zlim=z_range,col=gray(50:0/50),main="ktrans1")

dev.off()

name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/ktrans2_365.pdf"
pdf(name, width = 4, height = 4)
par(mfrow=c(1,1), mar=c(1, 1, 1, 1), cex=1)

    ktrans2_matrix <- t(matrix(ktrans2, nrow=J))
    z_range <- c(0,4)
    image.plot(ktrans2_matrix, axes=FALSE,zlim=z_range,col=gray(50:0/50),main="ktrans2")

dev.off()

name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/ktrans3_365.pdf"
pdf(name, width = 4, height = 4)
par(mfrow=c(1,1), mar=c(1, 1, 1, 1), cex=1)

    ktrans3_matrix <- t(matrix(ktrans3, nrow=J))
    z_range <- c(0,4)
    image.plot(ktrans3_matrix, axes=FALSE,zlim=z_range,col=gray(50:0/50),main="ktrans3")

dev.off()


z_range <- c(-0.001,0.2)
name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/trueEstiTheta_map_123vpComp.pdf"

#true_k <- c(1,70,102,129)     # 161 basis functions
true_k <- c(1,15,31,47)     # 61 basis functions

pdf(name, width = 16, height = 8)
par(mfrow=c(2,4), mar=c(1, 1, 1, 1), cex=1.8)

image(PAR_VP*sdX[1], axes=FALSE,zlim=z_range,col=gray(50:0/50),main="k=0")
box()
image(PAR_KTRANS1*sdX[15], axes=FALSE,zlim=z_range,col=gray(50:0/50),main="k=14")
box()
image(PAR_KTRANS2*sdX[31], axes=FALSE,zlim=z_range,col=gray(50:0/50),main="k=30")
box()
image(PAR_KTRANS3*sdX[47], axes=FALSE,zlim=z_range,col=gray(50:0/50),main="k=46")
box()

vp_matrix <- t(matrix(vp, nrow=J))
image(vp_matrix, axes=FALSE,zlim=z_range,col=gray(50:0/50))
box()
ktrans1_matrix <- t(matrix(theta1, nrow=J))
image(ktrans1_matrix, axes=FALSE,zlim=z_range,col=gray(50:0/50))
box()
ktrans2_matrix <- t(matrix(theta2, nrow=J))
image(ktrans2_matrix, axes=FALSE,zlim=z_range,col=gray(50:0/50))
box()
ktrans3_matrix <- t(matrix(theta3, nrow=J))
image(ktrans3_matrix, axes=FALSE,zlim=z_range,col=gray(50:0/50))
box()

dev.off()



name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/vp_365.pdf"
pdf(name, width = 4, height = 4)
par(mfrow=c(1,1), mar=c(1, 1, 1, 1), cex=1)

    vp_matrix <- t(matrix(vp, nrow=J))
    z_range <- c(0,0.2)
    image(vp_matrix, axes=FALSE,zlim=z_range,col=gray(50:0/50))

dev.off()

name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/theta1_365.pdf"
pdf(name, width = 4, height = 4)
par(mfrow=c(1,1), mar=c(1, 1, 1, 1), cex=1)

    ktrans1_matrix <- t(matrix(theta1, nrow=J))
    z_range <- c(-0.001,0.2)
    image(ktrans1_matrix, axes=FALSE,zlim=z_range,col=gray(50:0/50))

dev.off()

name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/theta2_365.pdf"
pdf(name, width = 4, height = 4)
par(mfrow=c(1,1), mar=c(1, 1, 1, 1), cex=1)

    ktrans2_matrix <- t(matrix(theta2, nrow=J))
    z_range <- c(-0.001,0.2)
    image(ktrans2_matrix, axes=FALSE,zlim=z_range,col=gray(50:0/50))

dev.off()

name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/theta3_365.pdf"
pdf(name, width = 4, height = 4)
par(mfrow=c(1,1), mar=c(1, 1, 1, 1), cex=1)

    ktrans3_matrix <- t(matrix(theta3, nrow=J))
    z_range <- c(-0.001,0.2)
    image(ktrans3_matrix, axes=FALSE,zlim=z_range,col=gray(50:0/50))

dev.off()





# ------------------------------------------------------------------------------
# plot curves
i <- 50
j <- 51

i <- 52
j <- 7

name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/sse_red_pix_365.pdf"
pdf(name, width = 4, height = 4)
par(mfrow=c(1,1), mar=c(4, 4, 4, 4), cex=1)

    #image.plot(conc[,,10], axes=FALSE,main="conc")
    sse_matrix <- t(matrix(sse, nrow=J))
    z_range <- c(0,0.15)
    #image(1:I,1:J,conc[,,10], main="conc",xlab="i",ylab="j")
    image.plot(1:I,1:J,sse_matrix, main="sse",zlim=z_range,col=gray(0:50/50),xlab="i",ylab="j")
    abline(h=j)
    abline(v=i)
    points(i,j,cex=2)

dev.off()



name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/Fit_365.pdf"
pdf(name, width = 8, height = 12)
par(mfrow=c(3,2), mar=c(4, 4, 1, 1), cex=1)
y_range <- c(-0.01,1)
for(i in c(50,52,70))
{
  for(j in c(7,51))
  {
      plot(time, Response[indx(c(i,j),I,J),], ylab="C(t)", ylim=y_range)
      points(time,spatFit$reducedFit[indx(c(i,j),I,J),],col="black", type="l", lwd=1.8)
  }
}
dev.off()

########################
# Fit, selection...
name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/Fit_Selection_365.pdf"
pdf(name, width = 12, height = 8)
par(mfrow=c(2,3), mar=c(4, 4, 1, 1), cex=1)
y_range <- c(-0.01,0.16)

#i <- 5  # q=1
#j <- 20 # without vp

i <- 5  # q=1
j <- 8 # without vp

plot(1:Q, spatFit$SpatialTheta[indx(c(i,j),I,J),], ylab=expression(theta), xlab=expression(k),type="h", ylim=y_range, main=paste("voxel",i,j),col="gray")
true_k <- c(1,c(15,38,45)+1)
points(true_k, c(PAR_VP[i,j],PAR_KTRANS1[i,j],PAR_KTRANS2[i,j],PAR_KTRANS3[i,j])*sdX[true_k], type="h", col="blue",lwd=1.5)
#points(1:Q,spatFit$reducedTheta[indx(c(i,j),I,J),]/sdX,col="red", type="h")
#points(1:Q, spatFit$SpatialTheta[indx(c(i,j),I,J),]/sdX, type="h") # nochmal dr?ber
points(1:Q,spatFit$reducedTheta[indx(c(i,j),I,J),],col="black", type="h",lwd=1.8)
#points(1:Q, spatFit$SpatialTheta[indx(c(i,j),I,J),], type="h") # nochmal dr?ber
abline(0,0)
legend("topright", # x-y coordinates for location of the legend
		cex=0.8,		# character expansion factor
		legend=c("True","Final estimates","First estimation step"),      # Legend labels
		col=c("blue","black","gray"),   # Color of points or lines
		lty=c(1,1,1),                    # Line type
		lwd=c(1.5,1.8,1))                    # Line width

#i <- 38  # q=2
#j <- 20 # without vp

i <- 38
j <- 8

plot(1:Q, spatFit$SpatialTheta[indx(c(i,j),I,J),], ylab=expression(theta), xlab=expression(k),type="h", ylim=y_range, main=paste("voxel",i,j),col="gray")
true_k <- c(1,c(15,38,45)+1)
points(true_k, c(PAR_VP[i,j],PAR_KTRANS1[i,j],PAR_KTRANS2[i,j],PAR_KTRANS3[i,j])*sdX[true_k], type="h", col="blue",lwd=1.5)
#points(1:Q,spatFit$reducedTheta[indx(c(i,j),I,J),]/sdX,col="red", type="h")
#points(1:Q, spatFit$SpatialTheta[indx(c(i,j),I,J),]/sdX, type="h") # nochmal dr?ber
points(1:Q,spatFit$reducedTheta[indx(c(i,j),I,J),],col="black", type="h",lwd=1.8)
#points(1:Q, spatFit$SpatialTheta[indx(c(i,j),I,J),], type="h") # nochmal dr?ber
abline(0,0)
legend("topright", # x-y coordinates for location of the legend
		cex=0.8,		# character expansion factor
		legend=c("True","Final estimates","First estimation step"),      # Legend labels
		col=c("blue","black","gray"),   # Color of points or lines
		lty=c(1,1,1),                    # Line type
		lwd=c(1.5,1.8,1))                    # Line width

#i <- 70  # q=3
#j <- 20 # without vp

i <- 52
j <- 8

plot(1:Q, spatFit$SpatialTheta[indx(c(i,j),I,J),], ylab=expression(theta), xlab=expression(k),type="h", ylim=y_range, main=paste("voxel",i,j),col="gray")
true_k <- c(1,c(15,38,45)+1)
points(true_k, c(PAR_VP[i,j],PAR_KTRANS1[i,j],PAR_KTRANS2[i,j],PAR_KTRANS3[i,j])*sdX[true_k], type="h", col="blue",lwd=1.5)
#points(1:Q,spatFit$reducedTheta[indx(c(i,j),I,J),]/sdX,col="red", type="h")
#points(1:Q, spatFit$SpatialTheta[indx(c(i,j),I,J),]/sdX, type="h") # nochmal dr?ber
points(1:Q,spatFit$reducedTheta[indx(c(i,j),I,J),],col="black", type="h",lwd=1.8)
#points(1:Q, spatFit$SpatialTheta[indx(c(i,j),I,J),], type="h") # nochmal dr?ber
abline(0,0)
legend("topright", # x-y coordinates for location of the legend
		cex=0.8,		# character expansion factor
		legend=c("True","Final estimates","First estimation step"),      # Legend labels
		col=c("blue","black","gray"),   # Color of points or lines
		lty=c(1,1,1),                    # Line type
		lwd=c(1.5,1.8,1))                    # Line width


y_range <- c(-0.01,0.75)
i <- 5  # q=1
j <- 8 # without vp

plot(time, Response[indx(c(i,j),I,J),], ylab="C(t)", xlab="t",ylim=y_range)
points(time, CONC_SIM[i,j,], col="blue", type="l", lty=1, lwd=1.8)
points(time,spatFit$reducedFit[indx(c(i,j),I,J),],col="black", type="l", lwd=1.8)


#i <- 38  # q=2
#j <- 20 # without vp
i <- 38
j <- 8

plot(time, Response[indx(c(i,j),I,J),], ylab="C(t)", xlab="t",ylim=y_range)
points(time, CONC_SIM[i,j,], col="blue", type="l", lty=1, lwd=1.8)
points(time,spatFit$reducedFit[indx(c(i,j),I,J),],col="black", type="l", lwd=1.8)



#i <- 70  # q=3
#j <- 20 # without vp
i <- 52
j <- 8

plot(time, Response[indx(c(i,j),I,J),], ylab="C(t)", xlab="t", ylim=y_range)
points(time, CONC_SIM[i,j,], col="blue", type="l", lty=1, lwd=1.8)
points(time,spatFit$reducedFit[indx(c(i,j),I,J),],col="black", type="l", lwd=1.8)

dev.off()


########################
# Fit, selection...
name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/Fit_Selection_simple_365.pdf"
pdf(name, width = 12, height = 8)
par(mfrow=c(2,3), mar=c(4, 4, 1, 1), cex=1)
y_range <- c(0,0.16)

#i <- 5  # q=1
#j <- 20 # without vp

i <- 5  # q=1
j <- 8 # without vp

true_k <- c(1,c(15,38,45)+1)
plot(1:Q, spatFit$reducedTheta[indx(c(i,j),I,J),], ylab=expression(theta),xlab=expression(k), type="h", lwd=1.8,ylim=y_range, main=paste("voxel",i,j))
points(true_k, c(PAR_VP[i,j],PAR_KTRANS1[i,j],PAR_KTRANS2[i,j],PAR_KTRANS3[i,j])*sdX[true_k], type="h", col="gray",lwd=1.8)
abline(a=0,b=0)

legend("topright", # x-y coordinates for location of the legend
		cex=0.8,		# character expansion factor
		legend=c("True","Estimates"),     # Legend labels
		col=c("gray","black"),   # Color of points or lines
		lty=c(1,1),                    # Line type
		lwd=c(1.8,1.8))                    # Line width

i <- 38  # q=2
j <- 8 # without vp

true_k <- c(1,c(15,38,45)+1)
plot(1:Q, spatFit$reducedTheta[indx(c(i,j),I,J),], ylab=expression(theta),xlab=expression(k), type="h", lwd=1.8,ylim=y_range, main=paste("voxel",i,j))
points(true_k, c(PAR_VP[i,j],PAR_KTRANS1[i,j],PAR_KTRANS2[i,j],PAR_KTRANS3[i,j])*sdX[true_k], type="h", col="gray",lwd=1.8)
abline(a=0,b=0)

legend("topright", # x-y coordinates for location of the legend
		cex=0.8,		# character expansion factor
		legend=c("True","Estimates"),     # Legend labels
		col=c("gray","black"),   # Color of points or lines
		lty=c(1,1),                    # Line type
		lwd=c(1.8,1.8))                    # Line width

i <- 52
j <- 8

true_k <- c(1,c(15,38,45)+1)
plot(1:Q, spatFit$reducedTheta[indx(c(i,j),I,J),], ylab=expression(theta),xlab=expression(k), type="h", lwd=1.8,ylim=y_range, main=paste("voxel",i,j))
points(true_k, c(PAR_VP[i,j],PAR_KTRANS1[i,j],PAR_KTRANS2[i,j],PAR_KTRANS3[i,j])*sdX[true_k], type="h", col="gray",lwd=1.8)
abline(a=0,b=0)

legend("topright", # x-y coordinates for location of the legend
		cex=0.8,		# character expansion factor
		legend=c("True","Estimates"),     # Legend labels
		col=c("gray","black"),   # Color of points or lines
		lty=c(1,1),                    # Line type
		lwd=c(1.8,1.8))                    # Line width


y_range <- c(-0.01,0.75)
i <- 5  # q=1
j <- 8 # without vp

plot(time, Response[indx(c(i,j),I,J),], ylab="C(t)", xlab="t", ylim=y_range)
points(time, CONC_SIM[i,j,], col="gray", type="l", lty=1, lwd=1.8)
points(time,spatFit$reducedFit[indx(c(i,j),I,J),],col="black", type="l", lwd=1.8)

#i <- 38  # q=2
#j <- 20 # without vp

i <- 38  # q=2
j <- 8 # without vp

#i <- 51
#j <- 51

plot(time, Response[indx(c(i,j),I,J),], ylab="C(t)", xlab="t", ylim=y_range)
points(time, CONC_SIM[i,j,], col="gray", type="l", lty=1, lwd=1.8)
points(time,spatFit$reducedFit[indx(c(i,j),I,J),],col="black", type="l", lwd=1.8)

#i <- 70  # q=3
#j <- 20 # without vp
i <- 52
j <- 8

plot(time, Response[indx(c(i,j),I,J),], ylab="C(t)", xlab="t", ylim=y_range)
points(time, CONC_SIM[i,j,], col="gray", type="l", lty=1, lwd=1.8)
points(time,spatFit$reducedFit[indx(c(i,j),I,J),],col="black", type="l", lwd=1.8)

dev.off()



name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/Residuals_365.pdf"
pdf(name, width = 8, height = 12)
par(mfrow=c(3,2), mar=c(4, 4, 1, 1), cex=1)
y_range <- c(-0.25,0.25)
for(i in c(5,35,70))
{
  for(j in c(5,70))
  {
      plot(time, Response[indx(c(i,j),I,J),]-spatFit$SpatialFit[indx(c(i,j),I,J),], ylab="residuals", type="p", ylim=y_range)
      abline(a=0,b=0)
  }
}
dev.off()

name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/reducedFit_365.pdf"
pdf(name, width = 8, height = 12)
par(mfrow=c(3,2), mar=c(4, 4, 1, 1), cex=1)
y_range <- c(-0.01,1)
for(i in c(5,35,70))
{
  for(j in c(5,70))
  {
      plot(time, Response[indx(c(i,j),I,J),], ylab="C(t)", ylim=y_range)
      points(time,spatFit$reducedFit[indx(c(i,j),I,J),],col="black", type="l", lwd=1.8)
  }
}
dev.off()


name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/reducedResiduals_365.pdf"
pdf(name, width = 8, height = 12)
par(mfrow=c(3,2), mar=c(4, 4, 1, 1), cex=1)
y_range <- c(-0.25,0.25)
for(i in c(5,35,70))
{
  for(j in c(5,70))
  {
      plot(time, Response[indx(c(i,j),I,J),]-spatFit$reducedFit[indx(c(i,j),I,J),], ylab="residuals", type="p", ylim=y_range)
      abline(a=0,b=0)
  }
}
dev.off()


# ------------------------------------------------------------------------------
# plot true parameter maps
load("save_simulation_123vpComp.Rdata")
PAR_VP[which(is.na(PAR_VP)==TRUE)] <- 0

z_range <- c(-0.001,0.2)
name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/trueTheta_map_123vpComp.pdf"

#true_k <- c(1,70,102,129)     # 161 basis functions
true_k <- c(1,15,31,47)     # 61 basis functions

pdf(name, width = 16, height = 4)
par(mfrow=c(1,4), mar=c(1, 1, 1, 1), cex=1)

image(PAR_VP*sdX[1], axes=FALSE,zlim=z_range,col=gray(50:0/50),main="k=1")
image(PAR_KTRANS1*sdX[15], axes=FALSE,zlim=z_range,col=gray(50:0/50),main="k=15")
image(PAR_KTRANS2*sdX[31], axes=FALSE,zlim=z_range,col=gray(50:0/50),main="k=31")
image(PAR_KTRANS3*sdX[47], axes=FALSE,zlim=z_range,col=gray(50:0/50),main="k=47")

dev.off()




# ------------------------------------------------------------------------------
# plot curves crossvalidation
setwd("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/")
load("spatFit_016_k1_i1_j1.Rdata")
name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/Fit_016_k1.pdf"
pdf(name, width = 8, height = 12)
par(mfrow=c(3,2), mar=c(4, 4, 1, 1), cex=1)
y_range <- c(-0.01,1)
for(i in c(5,35,70))
{
  for(j in c(5,70))
  {
      plot(time, Response[indx(c(i,j),I,J),], ylab="C(t)", ylim=y_range)
      points(time[-ik], Response[indx(c(i,j),I,J),-ik], col="red", ylim=y_range)
      points(time[-ik],spatFit$SpatialFit[indx(c(i,j),I,J),],col="black", type="l", lwd=1.8)

      predict <-  spatFit$reducedTheta%*%t(scaledX[ik,spatFit$reduced_k])

      points(time[ik], predict[indx(c(i,j),I,J),], col="blue")
  }
}
dev.off()






# Plot AIF
aif_TK <- Cp(time, model="tofts.kermode")
aif_FH <- Cp(time, model="fritz.hansen")

name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/Aif.pdf"
pdf(name, width = 4.0, height = 4.0)
par(mfrow=c(1,1), cex=1.2)

plot(time, aif_TK, col="black", type="l", lwd=1.5, ylab="C(t)")
points(time, aif_FH, col="grey", type="l", lwd=1.5)

legend("topright", # x-y coordinates for location of the legend
		cex=0.5,		# character expansion factor
		legend=c("Weinmann","FritzHansen"),     # Legend labels
		col=c("black", "grey"),   # Color of points or lines
		lty=c(1,1),                    # Line type
		lwd=c(1.5,1.5),                    # Line width
)

dev.off()

# Plot Predictors / Basis Functions
name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/Predictors.pdf"
pdf(name, width = 4.0, height = 4.0)
#par(mfrow=c(1,1), cex=1.2)
par(mfrow=c(1,1), mar=c(4, 4, 1, 1), cex=1.2)
plot(time, scaledX[,(Q-1)], col="black", type="l", lwd=1, ylab="Predictors", xlab="t")
for(k in seq(from=2,to=Q,by=5))
{
    points(time, scaledX[,k], col="black", type="l", lwd=1)
}
points(time, aif_TK, col="grey", type="l", lwd=1.5)
dev.off()

name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/Predictork2.pdf"
pdf(name, width = 4.0, height = 4.0)
par(mfrow=c(1,1), cex=1.2)


    plot(time, scaledX[,2], col="black", type="l", lwd=1,main="k=2")

#points(time, aif_TK, col="grey", type="l", lwd=1.5)
dev.off()

# kep <- infinity ausprobieren
aif_TK <- Cp(time, model="tofts.kermode")
kep <- 10000
ft <- model.weinmann(time,kep,model="tofts.kermode")

name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/Test_kep_inf.pdf"
pdf(name, width = 4.0, height = 4.0)
par(mfrow=c(1,1), cex=1.2)

plot(time, aif_TK, col="black", type="l", lwd=1.5, ylab="C(t)")
points(time, ft*10^4, col="grey", type="l", lwd=1.5)


dev.off()


# Plot simulated Curves
setwd("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/")
load("save_simulation_123vpComp.Rdata")

#y_range <- range(CONC_SIM_NOISE)
y_range <- c(-0.1,0.7)

name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/sim_curves.pdf"
pdf(name, width = 4.0, height = 4.0)
par(mfrow=c(1,1), mar=c(4, 4, 1, 1), cex=1.2)

plot(time, CONC_SIM[5,5,], col="grey", type="l", lwd=1.95, xlab="t",ylab="C(t)",ylim=y_range)
points(time, CONC_SIM[35,5,], col="black", type="l", lty=1, lwd=1.95)
points(time, CONC_SIM[70,5,], col="blue", type="l", lty=1, lwd=1.95)

points(time, CONC_SIM_NOISE[5,5,], col="grey", type="p", pch=20, lwd=0.5)
points(time, CONC_SIM_NOISE[35,5,], col="black", type="p", pch=20, lwd=0.5)
points(time, CONC_SIM_NOISE[70,5,], col="blue", type="p", pch=20, lwd=0.5)

# # Add a legend to the plot
legend("bottomright", # x-y coordinates for location of the legend
		cex=1,		# character expansion factor
		legend=c("q=1","q=2","q=3"),     # Legend labels
		col=c("grey", "black", "blue"),   # Color of points or lines
		lty=c(1,1,1),                    # Line type
		lwd=c(1.95,1.95,1.95),                    # Line width
)

dev.off()


name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/sim_curves_with_vp.pdf"
pdf(name, width = 4.0, height = 4.0)
par(mfrow=c(1,1), cex=1.2)

plot(time, CONC_SIM[5,75,], col="grey", type="l", lwd=1.5, xlab="t",ylab="C(t)",ylim=y_range)
points(time, CONC_SIM[35,75,], col="black", type="l", lty=1, lwd=1.5)
points(time, CONC_SIM[70,75,], col="blue", type="l", lty=1, lwd=1.5)

points(time, CONC_SIM_NOISE[5,75,], col="grey", type="p", pch=20, lwd=1)
points(time, CONC_SIM_NOISE[35,75,], col="black", type="p", pch=20, lwd=1)
points(time, CONC_SIM_NOISE[70,75,], col="blue", type="p", pch=20, lwd=1)

# # Add a legend to the plot
legend("bottomright", # x-y coordinates for location of the legend
		cex=0.5,		# character expansion factor
		legend=c("q=1","q=2","q=3"),     # Legend labels
		col=c("grey", "black", "blue"),   # Color of points or lines
		lty=c(1,1,1),                    # Line type
		lwd=c(1.5,1.5,1.5),                    # Line width
)

dev.off()








selected <- which(Theta2[n,]>10^(-10),arr.ind=TRUE)
        #print("selected")
        #print(selected)

selected <- c(1,2,3,49, 54, 55, 56, 61)

selected <- c(1,2,3,49, 54, 55, 56)

        selected_1 <- c(Q,selected[-length(selected)]) # right shift
        selected_2 <- c(selected[-1],Q+1) # left shift

        left_bounds <- which(selected!=selected_1+1)
        right_bounds <- which(selected!=selected_2-1)

#        print("bounds")
#        print(left_bounds)
#        print(right_bounds)

        selected_max <- c()
        for(k in 1:length(left_bounds))
        {
            index_max <- which.max(Theta2[n,selected[left_bounds[k]:right_bounds[k]]])   # index within k-th group
            index_max <- left_bounds[k] + index_max -1     # index within selected
            selected_max <- c(selected_max,selected[index_max])
        }


index <- c(50,2)
Selected_array(indx(index,I,J))

test3 <- callSpatEnetSel(index, I=I, J=J, x=x, Response=Response,
  Theta=reducedTheta, lambda=lambda2, sv=sv2,firstpen=TRUE,Selected=Selected_array)


N <- I*J
fit <- array(NA,c(N,T))
for(n in 1:N)
{
    fit[n,] <-  spatFit$reducedTheta[n,spatFit$selected[[n]]]%*%t(scaledX[,spatFit$selected[[n]]])

}



############################
## voxelwise plots
#############################
load("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/sim_Theta_voxelwise_016.Rdata")


z_range <- c(-0.001,0.05)
name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/Theta_map_voxelwise_016.pdf"
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

name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatial/sse_Fit_voxelwise_016.pdf"
pdf(name, width = 4, height = 4)
par(mfrow=c(1,1), mar=c(1, 1, 1, 1), cex=1)

    sse_matrix <- t(matrix(sse, nrow=J))
    z_range <- c(0,0.15)
    image.plot(sse_matrix, axes=FALSE,zlim=z_range,col=gray(0:50/50),main="sse")

dev.off()

################
# AIC
sse <- array(NA, I*J)
q_number <- array(NA, I*J)
error <- Response - spatFit$reducedFit
for(n in 1:N)
{
	 sse[n] <- sum(error[n,]^2)
	 q_number[n] <- length(which(spatFit$reducedTheta[n,]!=0))
}

# AIC = -2log-likeli+2q
true_sigma <- 0.05
AIC <- sum(sse)/(N*T*true_sigma^2) + 2*sum(q_number)/(N*T)



