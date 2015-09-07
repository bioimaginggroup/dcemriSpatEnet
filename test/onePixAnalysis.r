# data
#load("Z:/Research/DCE-MRI/Response.RData")
load("Z://Dateien//Projekte//Jan-Gertheis//Response.RData")

# functions
#source("Z:/Research/DCE-MRI/onePixFunctions.r")
source("Z://Dateien//Projekte//Jan-Gertheis//onePixFunctions.r")

# data visualization

Images <- vector("list",46)
npix <- 41*41
cuts <- cbind(seq(1,npix*45+1,by=npix),seq(npix,npix*46,by=npix))

for (tp in 1:46)
  {
    Images[[tp]] <- matrix(MRIdata$Ct[cuts[tp,1]:cuts[tp,2]],41,41)
  }

xp <- 29
yp <- 30
time <- seq(0,9.1,length=46)

maxInt <- max(unlist(lapply(Images,max)))
minInt <- min(unlist(lapply(Images,min)))
plotimages(tps=1:46, wait=0.5, x=xp, y=yp, time=time, limits=c(minInt,maxInt))

tps <- c(5,13,25,37)
#postscript("Z://Dateien//Projekte//Jan-Gertheis//plotInt1.ps",height=8,width=8,
#horizontal=F)
pdf("Z://Dateien//Projekte//Jan-Gertheis//plotInt1.pdf",height=8,width=8)
par(mfrow=c(2,2))
plotimage(tp=tps[1],x=xp,y=yp,limits=c(minInt,maxInt))
title(paste("time = ", round(time[tps[1]],2), sep=""))
plotimage(tp=tps[2],x=xp,y=yp,limits=c(minInt,maxInt))
title(paste("time = ", round(time[tps[2]],2), sep=""))
plotimage(tp=tps[3],x=xp,y=yp,limits=c(minInt,maxInt))
title(paste("time = ", round(time[tps[3]],2), sep=""))
plotimage(tp=tps[4],x=xp,y=yp,limits=c(minInt,maxInt))
title(paste("time = ", round(time[tps[4]],2), sep=""))
par(mfrow=c(1,1))
dev.off()

# concentration time curve
Response <- as.vector(MRIdata[(MRIdata$x==xp)&(MRIdata$y==yp),1])

# Xmatrix
Matrix <- Cp(time) # First curve: plasma concentration
# basis function for different kep values
for (kep in exp(seq(-2,2,0.01)))
  {    
    ft <- model.weinmann(time,kep)
    Matrix <- cbind(Matrix,ft)
}

#postscript("Z:/Research/DCE-MRI/Paper/plotVox1.ps",height=6,width=6,
#horizontal=F)
plot(Matrix[,1], type="l", ylim=c(0,3.3), xlab="time point",
ylab="response / predictors")
for (i in 2:ncol(Matrix))
lines(Matrix[,i])
lines(Response,col="blue",lwd=2)
title(paste("x =", xp, ", y =", yp))
#dev.off()

# ------------------------------------------------------------------------------

# data for nnent
p <- ncol(Matrix)

sdX <- apply(Matrix,2,sd)
scaledX <- scale(Matrix,center=F,scale=sdX)

# model fitting
nnenet1 <- nnenet(scaledX, Response, 10^(-5), 0.125, firstpen=T)                
        
# plots
postscript("Z:/Research/DCE-MRI/Paper/plotCoef1a.ps",width=10,height=4,
horizontal=F)
plot(seq(-2,2,by=0.01),(nnenet1$beta/sdX)[-1],type="h",xlab=expression(log(k[ep])),
ylab=expression(K^{trans}))
dev.off()

postscript("Z:/Research/DCE-MRI/Paper/plotFit1a.ps",width=5,height=5,
horizontal=F)
plot(Response,nnenet1$fit,xlab="observed concentrations",
ylab="fitted concentrations")
abline(0,1,lty=2)
dev.off()

# non-zero betas
which(nnenet1$beta > 10^(-10))

# ------------------------------------------------------------------------------

# tuning parameter evaluation 
xypix <- cbind(sort(rep(1:41,41)),rep(1:41,41))

pixsample <- 1:1681

image(1:41,1:41,Images[[1]],xlab="x",ylab="y")
points(xypix[pixsample,],lwd=2)

colnames(xypix) <- c("x","y")
write.table(xypix[pixsample,],"Z:/Research/DCE-MRI/CVscores1/xysample.txt")

lambda <- 10^seq(-8,-3,by=0.5)
#lambda <- 10^(-8)
svalues <- c(0.01,0.015,0.02,0.03,0.04,0.055,seq(0.075,0.3,by=0.025))

for (pix in pixsample) # dauert sehr lange
  {
    cat("pix",pix,"\n")
    xp <- xypix[pix,1]
    yp <- xypix[pix,2]
    Response <- as.vector(MRIdata[(MRIdata$x==xp)&(MRIdata$y==yp),1])
    CVscore <- CVnnenet(K=5, x=scaledX, y=Response, lambda=lambda,
    svalues=svalues, printFold=F, firstpen=F)
    write.table(CVscore,
    paste("Z:/Research/DCE-MRI/CVscores1/CVscore",xp,yp,".txt",sep="q"))
  }

# ------------------------------------------------------------------------------

# fitting 
lambda <- 10^seq(-8,-3,by=0.5)
svalues <- c(0.01,0.015,0.02,0.03,0.04,0.055,seq(0.075,0.3,by=0.025))

sdX <- apply(Matrix,2,sd)
scaledX <- scale(Matrix,center=F,scale=sdX)

xypix <- cbind(sort(rep(1:41,41)),rep(1:41,41))
colnames(xypix) <- c("x","y")
xysample <- xypix

Coefs <- matrix(NA,nrow(xysample),ncol(scaledX))
Fit <- matrix(NA,nrow(xysample),nrow(scaledX))

for (pix in 1:nrow(xysample)) # dauert auch etwas
  {
    if(is.element(pix,seq(20,1660,by=20)))
     print(pix)
    xp <- xysample[pix,1]
    yp <- xysample[pix,2]
    CVscore <- read.table(paste("Z:/Research/DCE-MRI/CVscores1/CVscore",
    xp,yp,".txt",sep="q"))
    #whichlam <- which.min(apply(CVscore,2,min))
    #bestlam <- lambda[whichlam]
    #bests <- svalues[which.min(CVscore[,whichlam])]
    bestlam <- 10^(-8)
    bests <- svalues[which.min(CVscore[,1])]
    
    Response <- as.vector(MRIdata[(MRIdata$x==xp)&(MRIdata$y==yp),1])
    model <- nnenet(scaledX, Response, bestlam, bests, firstpen=T)
    Coefs[pix,] <- model$beta
    Fit[pix,] <- model$fit
  }
  
Coefs <- round(t(t(Coefs)/sdX),digits=10)


# analyze results
VarFit <- numeric(nrow(Fit))
Residuals <- matrix(NA,nrow(Fit),ncol(Fit))
RSS <- numeric(nrow(Fit))
TSS <- numeric(nrow(Fit))
autocorr <- numeric(nrow(Fit))
n <- ncol(Fit)
for (pix in 1:nrow(Residuals))
  { 
    xp <- xysample[pix,1]
    yp <- xysample[pix,2]
    Response <- as.vector(MRIdata[(MRIdata$x==xp)&(MRIdata$y==yp),1])
    VarFit[pix] <- sum((Fit[pix,] - mean(Fit[pix,]))^2)
    TSS[pix] <- sum((Response - mean(Response))^2)
    Residuals[pix,] <- Response - Fit[pix,]
    RSS[pix] <- sum(Residuals[pix,]^2)
    autocorr[pix] <- cor(Residuals[pix,-1],Residuals[pix,-n])
  }

# estimated variances
postscript("Z:/Research/DCE-MRI/Paper/plotVar1.ps",height=6,width=6,
horizontal=F)
image(1:41,1:41,matrix(RSS/ncol(Fit),41,41,byrow=T),xlab="x",ylab="y",
zlim=c(0,0.0015),main="Estimated Variances")
dev.off()

# autocorrelations
postscript("Z:/Research/DCE-MRI/Paper/plotAutoCorr1a.ps",height=6,width=6,
horizontal=F)
hist(autocorr,main="Histogram of Observed Auto-Correlations",
xlab="observed values", ylab="frequency")
dev.off()

postscript("Z:/Research/DCE-MRI/Paper/plotAutoCorr1b.ps",height=6,width=6,
horizontal=F)
image(1:41,1:41,matrix(autocorr,41,41,byrow=T),xlab="x",ylab="y", zlim=c(-1,1),
main="Spatial Pattern of Observed Auto-Correlations")
dev.off()

# selected coefficients
Selected <- vector("list", nrow(Coefs))
for (pix in 1:length(Selected))
  {
    Selected[[pix]] <- which(Coefs[pix,-1] > 0)
  }
  
# numbner of selected compartments
NBlocks <- numeric(length(Selected))
for (pix in 1:length(NBlocks))  
  {
    kepp <- Selected[[pix]]
    if (length(kepp)>1)
      NBlocks[pix] <- sum((kepp[-1] - kepp[-length(kepp)]) > 1) + 1
    else
      NBlocks[pix] <- 1
  } 
 
table(NBlocks) 

NBlocks3 <- NBlocks
NBlocks3[NBlocks > 3] <- 3

image(1:41,1:41,matrix(NBlocks3,41,41,byrow=T),xlab="x",ylab="y",#zlim=c(1,6),
main="Number of Selected Compartments")

# compare with one-compartment model
oCCoefs <- matrix(NA,nrow(xysample),3)
oCRSS <- numeric(nrow(xysample))

for (pix in 1:nrow(xysample))
  {
    if(is.element(pix,seq(20,1660,by=20)))
     print(pix)
    xp <- xysample[pix,1]
    yp <- xysample[pix,2]

    Response <- as.vector(MRIdata[(MRIdata$x==xp)&(MRIdata$y==yp),1])
    oCmodel <- optim(par = c(0,0,0), fn = oneCompError, CT = Response, 
    time = time, method= "L-BFGS-B", lower = 0)

    oCCoefs[pix,] <- oCmodel$par
    oCRSS[pix] <- oCmodel$value
  }

par(mfrow=c(2,2))
hist((RSS - oCRSS)/TSS, main="RSS(m)/TSS - RSS(1)/TSS", xlab="observed values")
hist((RSS - oCRSS)/ncol(Fit), main="Var(m) - Var(1)", xlab="observed values")

summary((RSS - oCRSS)/TSS)
summary((RSS - oCRSS)/ncol(Fit))

image(1:41,1:41,matrix((RSS-oCRSS)/TSS,41,41,byrow=T),xlab="x",ylab="y",
main="RSS(m)/TSS - RSS(1)/TSS")
image(1:41,1:41,matrix((RSS-oCRSS)/ncol(Fit),41,41,byrow=T),xlab="x",ylab="y",
main="Var(m) - Var(1)")


hist(RSS/TSS, main="RSS(m)/TSS", xlab="observed values", breaks=20)
hist(RSS/ncol(Fit), main="Var(m)", xlab="observed values", xlim=c(0,0.002), breaks=20)
hist(oCRSS/TSS, main="RSS(1)/TSS", xlab="observed values", breaks=20)
hist(oCRSS/ncol(Fit), main="Var(1)", xlab="observed values", xlim=c(0,0.002), breaks=20)

par(mfrow=c(1,1))

xps <- c(10,18,29,33) 
yps <- c(5,27,30,23)

postscript("Z:/Research/DCE-MRI/Paper/plotCurves1.ps",height=8,width=8,
horizontal=F)
par(mfrow=c(2,2))
for (i in 1:4)
  {
    xp <- xps[i]
    yp <- yps[i]
    Response <- as.vector(MRIdata[(MRIdata$x==xp)&(MRIdata$y==yp),1])
    oCC <- oCCoefs[(xysample[,1]==xp)&(xysample[,2]==yp),]
    oCFit <- oneComp(time, oCC[1], oCC[2], oCC[3])

    plot(time,Response,type="p",ylim=c(0,0.5),xlab="time",ylab="concentration",
    pch=1,cex=0.7)
    lines(time,Fit[(xysample[,1]==xp) & (xysample[,2]==yp),],col=2,lwd=2)
    lines(time,oCFit,lwd=2,lty=2)
    title(paste("voxel = {", xp, ",", yp, "}", sep=""))
  }
dev.off()
