setwd("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/")            ## Server## breast cancer data
library(oro.nifti)


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

    # least squares fit of 1Comp model
    require("dcemriS4")
    require("fields")
    require("minpack.lm")
    library(multicore,lib="~/R/x86_64-pc-linux-gnu-library/2.12/")
    options(cores=4)
    multicore=TRUE
    
    Fit1Comp <- dcemri.lm(conc, time, mask, model="extended", aif="tofts.kermode", multicore=TRUE, verbose=TRUE)
    
    name <- paste("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/1Comp/lsFit_scan", var_s, sep="")
    name <- paste(name,"_3.Rdata",sep="")
    
    save(Fit1Comp,file=name)
   

    # least squares fit extended 1Comp model
    name <- paste("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/1Comp/lsFit_scan", var_s, sep="")
    name <- paste(name,"_3.Rdata",sep="")    
    load(name)
    
    require("fields")
    
    
    
    name <- paste("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/1Comp/vp_map_scan", var_s, sep="")
    name <- paste(name,"_3.pdf",sep="")
    
    z_range <- range(Fit1Comp$vp[,,1],na.rm=TRUE)
    
    pdf(name, width = 4, height = 4)
    par(mfrow=c(1,1), mar=c(1, 1, 1, 1), cex=1)
        image(Fit1Comp$vp[,,1], axes=FALSE,zlim=z_range,col=gray(50:0/50),main="vp")
    dev.off()
    
    
    name <- paste("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/1Comp/ktrans_map_scan", var_s, sep="")
    name <- paste(name,"_3.pdf",sep="")
    
    z_range <- c(0,2)
    
    pdf(name, width = 4, height = 4)
    par(mfrow=c(1,1), mar=c(1, 1, 1, 1), cex=1)
        image(Fit1Comp$ktrans[,,1], axes=FALSE,zlim=z_range,col=gray(50:0/50),main="ktrans")
    dev.off()
    
    
    name <- paste("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/1Comp/kep_map_scan", var_s, sep="")
    name <- paste(name,"_3.pdf",sep="")
    
    z_range <- c(0,2)
    
    pdf(name, width = 4, height = 4)
    par(mfrow=c(1,1), mar=c(1, 1, 1, 1), cex=1)
        image(Fit1Comp$kep[,,1], axes=FALSE,zlim=z_range,col=gray(50:0/50),main="kep")
    dev.off()
    
    
    name <- paste("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/1Comp/sse_map_scan", var_s, sep="")
    name <- paste(name,"_3.pdf",sep="")
    
    z_range <- c(0,0.15)
    
    pdf(name, width = 4, height = 4)
    par(mfrow=c(1,1), mar=c(1, 1, 1, 1), cex=1)
        image.plot(Fit1Comp$sse[,,1], axes=FALSE,zlim=z_range,col=gray(0:50/50),main="sse")
    dev.off()

}



