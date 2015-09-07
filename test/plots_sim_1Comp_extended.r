setwd("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/")            ## Server## breast cancer data
library(oro.nifti)

# simulated data


load("save_simulation_123vpComp_extended.Rdata")
I <- 75 # numer of rows
J <- 75 # numer of columns
N <- I*J # number of pixels
K <- 30  # number of simulated images


conc <-  CONC_SIM_NOISE
mask <-  array(TRUE, c(I,J))


# least squares fit of 1Comp model
require("dcemriS4")
require("fields")
require("minpack.lm")
library(multicore,lib="~/R/x86_64-pc-linux-gnu-library/2.12/")
options(cores=5)
multicore=FALSE

Fit1Comp_list <- c()

for(k in 1:K)
{

   Fit1Comp <- dcemri.lm(conc[,,,k], time, mask, model="extended", aif="tofts.kermode", multicore=TRUE, verbose=TRUE)               
   Fit1Comp_list[[k]] <- Fit1Comp   
}   


name <- paste("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/1Comp/lsFit_sim")
name <- paste(name,"_extended.Rdata",sep="")

save(Fit1Comp_list,file=name)


# least squares fit extended 1Comp model
name <- paste("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/1Comp/lsFit_sim")
name <- paste(name,"_extended.Rdata",sep="")    
load(name)

########################################################################
# plots 1Comp
require("fields")

#name <- paste("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/1Comp/vp_choice_sim")
#name <- paste(name,"_2.pdf",sep="")
#z_range <- c(0,1)
#pdf(name, width = 4, height = 4)
#par(mfrow=c(1,1), mar=c(1, 1, 1, 1), cex=1)
#    image(Fit1Comp$vp[,,1]>0, axes=FALSE,zlim=z_range,col=gray.colors(2, start=1, end=0.3))
#    #image(Fit1Comp$vp[,,1]>10^(-5), axes=FALSE,zlim=z_range,col=gray.colors(2, start=1, end=0.3))
#    box()
#dev.off()
#
#name <- paste("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/1Comp/ktrans_choice_sim")
#name <- paste(name,"_2.pdf",sep="")
#
#z_range <- c(0,4)
#
#pdf(name, width = 4, height = 4)
#par(mfrow=c(1,1), mar=c(1, 1, 1, 1), cex=1)
#    image(Fit1Comp$ktrans[,,1]>0, axes=FALSE,zlim=z_range,col=gray.colors(4, start=0.95, end=0))
#    box()
#dev.off()
#
#
#name <- paste("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/1Comp/vp_map_sim")
#name <- paste(name,"_2.pdf",sep="")
#
#z_range <- range(Fit1Comp$vp[,,1],na.rm=TRUE)
#
#pdf(name, width = 4, height = 4)
#par(mfrow=c(1,1), mar=c(1, 1, 1, 1), cex=1)
#    image(Fit1Comp$vp[,,1], axes=FALSE,zlim=z_range,col=gray(50:0/50),main="vp")
#dev.off()
#
#name <- paste("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/1Comp/ktrans_map_sim")
#name <- paste(name,"_2.pdf",sep="")
#
#z_range <- c(0,2)
#
#pdf(name, width = 4, height = 4)
#par(mfrow=c(1,1), mar=c(1, 1, 1, 1), cex=1)
#    image(Fit1Comp$ktrans[,,1], axes=FALSE,zlim=z_range,col=gray(50:0/50),main="ktrans")
#dev.off()
#
#
#name <- paste("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/1Comp/kep_map_sim")
#name <- paste(name,"_2.pdf",sep="")
#
#z_range <- c(0,2)
#
#pdf(name, width = 4, height = 4)
#par(mfrow=c(1,1), mar=c(1, 1, 1, 1), cex=1)
#    image(Fit1Comp$kep[,,1], axes=FALSE,zlim=z_range,col=gray(50:0/50),main="kep")
#dev.off()
#


  msse <- array(NA, c(I,J))   
  sse <- array(NA, c(I,J,T,K))
  
  
  for(k in 1:K)
  {
     sse[,,,k] <-  Fit1Comp_list[[k]]$sse
  }
  
  #error <- Response - spatFit$SpatialFit
  for(i in 1:I)
  {
    for(j in 1:J)
    {
      	 msse[i,j] <-  sum(sse[i,j,1,])/K
    }
  }




name <- paste("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/1Comp/msse_map_sim")
name <- paste(name,"_extended.pdf",sep="")

z_range <- c(0,0.25)

pdf(name, width = 4, height = 4)
par(mfrow=c(1,1), mar=c(1, 1, 1, 1), cex=1)
    image(msse, axes=FALSE,zlim=z_range,col=gray(0:50/50))
    box()
dev.off()


  #######################################
  # extended 1Comp
  #name <- paste("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/1Comp/lsFit_scan", var_s, sep="")
  #name <- paste(name,"_2.Rdata",sep="")
  #load(name)
  
  # least squares fit extended 1Comp model
  name <- paste("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/1Comp/lsFit_sim")
  name <- paste(name,"_extended.Rdata",sep="")    
  load(name)


#  q_1Comp <- c()
#  sse_1Comp <- c()
#  aic_1Comp <- c()
#  bic_1Comp <- c()

  aic <- array(NA,K)
  bic <- array(NA,K)

  for(k in 1:K)
  {
      sse <- Fit1Comp_list[[k]]$sse[,,1]
      
      true_sigma <- 0.05  # for simulated data
     
      # AIC = -2log-likeli+2q  
      q <- 2*N # fixed model, vp + ktrans1 per voxel
      aic[k] <- sum(sse,na.rm=TRUE)/(true_sigma^2) + 2*q  
      bic[k] <- sum(sse,na.rm=TRUE)/(true_sigma^2) + log(N*T) * q
      
      #q_1Comp <- c(q_1Comp,q)
#      sse_1Comp <- c(sse_1Comp,sum(sse,na.rm=TRUE))
#      aic_1Comp <- c(aic_1Comp,aic)
#      bic_1Comp <- c(bic_1Comp,bic)
#      
#      print(q_1Comp)
#      print(round(sse_1Comp))
#      print(round(aic_1Comp))
#      print(round(bic_1Comp))
  }


########################################################################
# plots true
setwd("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/") 
load("save_simulation_123vpComp.Rdata")
I <- 75 # numer of rows
J <- 75 # numer of columns
N <- I*J # number of pixels

true_sigma <- 0.05
T <- length(time)

PAR_VP[which(is.na(PAR_VP)==TRUE)] <- 0

PAR_KTRANS1[which(PAR_KTRANS1>0)] <- 1
PAR_KTRANS2[which(PAR_KTRANS2>0)] <- 1
PAR_KTRANS3[which(PAR_KTRANS3>0)] <- 1

PAR_Q <- PAR_KTRANS1 + PAR_KTRANS2 + PAR_KTRANS3


name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/1Comp/sse_true.pdf"
z_range <- c(0,0.25)
pdf(name, width = 4, height = 4)
par(mfrow=c(1,1), mar=c(1, 1, 1, 1), cex=1)
    image(array((true_sigma)^2*T,c(I,J)), axes=FALSE,zlim=z_range,col=gray(0:50/50))
    box()
dev.off()


name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/1Comp/vp_true.pdf"
z_range <- c(0,1)
pdf(name, width = 4, height = 4)
par(mfrow=c(1,1), mar=c(1, 1, 1, 1), cex=1)
    image(PAR_VP>0, axes=FALSE,zlim=z_range,col=gray.colors(2, start=1, end=0.3))
    box()
dev.off()

name <- "~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/1Comp/q_true.pdf"
z_range <- c(0,4)
pdf(name, width = 4, height = 4)
par(mfrow=c(1,1), mar=c(1, 1, 1, 1), cex=1)
    image(PAR_Q, axes=FALSE,zlim=z_range,col=gray.colors(4, start=0.95, end=0))
    box()
dev.off()