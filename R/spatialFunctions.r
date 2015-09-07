# get number from two-dimensional index
indx <- function(index,I,J)
{ 
  help <-  (index[1]-1)*J + index[2]
  return(help)
}


# function that returns list of indices of all neighbors of Pixel p = (x_p, y_p)
# four nearest neighbors within image of dimension P1 x P2
getNeighbors <- function(p, P1, P2, Response)
{
    if(length(p)!=2){print("wrong dimension of pixel coordinates")
                      print(p)}
    
    xp <- p[1]
    yp <- p[2]
    
    neighbors <- list()
    count <- 1
    if(xp-1 > 0 && !is.na(Response[indx(c(xp-1,yp),I,J),1]))
    {
        neighbors[[count]] <- c(xp-1,yp)
        count <- count+1
    }
    if(xp+1 <= P1 && !is.na(Response[indx(c(xp+1,yp),I,J),1]))
    {
        neighbors[[count]] <- c(xp+1,yp)
        count <- count+1
    }
    if(yp-1 >0 && !is.na(Response[indx(c(xp,yp-1),I,J),1]))
    {
        neighbors[[count]] <- c(xp,yp-1)
        count <- count+1
    }
    if(yp+1 <= P2 && !is.na(Response[indx(c(xp,yp+1),I,J),1]))
    {
        neighbors[[count]] <- c(xp,yp+1)
        count <- count+1
    }
    
    return(neighbors)    
}

# theta: List of I*J pixel for Q parameters
getXi <- function(theta, neighbors, I, J, firstpen)
{
    n <- length(neighbors)   # number of neighbors
    Xi <- c()
    
    for(k in 1:n)
      {
          Xi <- c(Xi,theta[indx(neighbors[[k]],I,J),])    
      }
                
    if(!firstpen)
    {   
        # parameters corresponding to Cp are not penalized
        p <- length(Xi)/n
        for(k in 1:n)
        {
            Xi[(k-1)*p+1] <- 0
        }   
    }    
    return(Xi)
}

# q: length of theta, number of basis functions
getD <- function(q,neighbors, firstpen)
{
    n <- length(neighbors)
    if(firstpen)
    {
        I <- diag(1,nrow=q)
    }
    else
    {
        I <- diag(c(0,rep(1,q-1)))
    }
        
    D <- c()
    for(k in 1: n)
    {
      D <- rbind(D,I)
    }    
    return(D)
}


# function to be called for spatial fit of one voxel
# x: Matrix of basis functions scaled such that entries have unit variance
# y: Response, Concentration (for one pixel)
# lambda: strength of penalization (Ridge term)
# sv: penalization parameter for Lasso term
# firstpen: indicates, if first term (vp-term) is penalized
# Xi: vector of parameter values from neighboring pixels
# D: design matrix extension for pseudo-observations
spatEnet <- function(x, y, lambda, sv, firstpen, Xi, D)
{
    x0 <- x
  
    T <- length(y)   # T (length of time)
    q <- ncol(x)     # q (number of basis functions)
      
    # extent observation vector with pseudo-observations
    y <- c(y,sqrt(lambda)*Xi)
          
    # extent design matrix    
    x <- rbind(x,sqrt(lambda)*D)

    dvec <- t(x)%*%y
    Dmat <- t(x)%*%x  
    if(firstpen)
    {
      Amat <- t(rbind(diag(1,q),rep(-1,q)))
    }
    else
    {
      Amat <- t(rbind(diag(1,q),c(0,rep(-1,q-1))))
    }
    
    bvec <- c(rep(0,q),-sv)
    
    nnen <- solve.QP(Dmat=Dmat, dvec=dvec, Amat=Amat, bvec=bvec)
    
    if(warning()!="")
    {
         print(warning())
    }
    #nnen <- (1+lambda)*nnen$solution
    nnen <- nnen$solution
    return(list("beta"=nnen,"fit"=x0%*%nnen))
}

# function to get all information needed to call spatEnet parallelized
# index: pixel index
# I, J: dimension of original Response-Array
# Response: list of Responses
# Theta: list of parameters per pixel
callSpatEnet <- function(index, I, J, x, Response, Theta, lambda, sv, firstpen=TRUE)
{

   if(is.na(Response[indx(index,I,J),1]))
   {
        #print(index)
        #print("is.na")
        return(list("beta"=NA,"fit"=NA))
   }

   y <- Response[indx(index,I,J),]
   
   neighborList <- getNeighbors(index,I,J,Response)
     
   Xi <- getXi(Theta, neighborList, I, J, firstpen)
   
   x <- as.matrix(x)
   
   D <- getD(ncol(x),neighborList, firstpen)
     
   theta_new <- spatEnet(x, y, lambda, sv, firstpen, Xi, D)
   
   return(theta_new)
    
}

callSpatEnetSel <- function(index, I, J, x, Response, Theta, lambda, sv, firstpen=TRUE, Selected)
{
   if(is.na(Response[indx(index,I,J),1]))
   {             
        return(list("beta"=NA,"fit"=NA))        
   }
    
   y <- Response[indx(index,I,J),]
   
   neighborList <- getNeighbors(index,I,J,Response)
      
   n <- indx(index,I,J)
   
   Theta_n <- as.matrix(Theta[,Selected[[n]]])   
   x_n <- as.matrix(x[,Selected[[n]]])
          
   Xi <- getXi(Theta_n, neighborList, I, J, firstpen)
   D <- getD(ncol(x_n),neighborList, firstpen)
        
   theta_new <- spatEnet(x_n, y, lambda, sv, firstpen, Xi, D)
      
   return(theta_new)
    
}

# function for (reduced) spatial fit
# iteratively update of pixels (checkerboard pattern)
# first with many basis functions
# refit with reduced basis
# Theta: starting value for parameters
# return: Theta2, Fit2 (full spatial fit), reducedTheta, reducedFit
# centersOnly: TRUE -> interval centers as reduced basis
#              FALSE -> all 'non-vanishing' basis functions 
spatialFit <- function(xyblack_list, xywhite_list, I, J, x, Response, time,
                        Theta, lambda, sv, firstpen, centersOnly)
{
    
    print("in spatialFit")
    print(c(lambda,sv))
    
    T <- length(time)
    Q <- ncol(x)
    N <- I*J
    Theta2 <- Theta
    Fit2 <- array(NA, c(I*J,T))

    delta <- 100
    count <- 0

    max_it <- 2500
    min_it <- 100

#    #while(delta > 3)
#    #while(delta > (N*0.0005) & count < max_it )  
#    while(delta > (N*0.0002) & count < max_it | count < min_it)  
#    while(delta > (N*sv/1000) & count < max_it | count < min_it) 
    while(delta > 0.0001 & count < max_it | count < min_it)
    {   
        Theta_old <- Theta2

        # parallel update black pixels
        test3 <- mclapply(xyblack_list, FUN = callSpatEnet, I=I, J=J, x=x, Response=Response,
                        Theta=Theta2, lambda=lambda, sv=sv,firstpen=firstpen)
        # update Theta2-values
        for(i in 1:length(xyblack_list))
        {
             Theta2[indx(xyblack_list[[i]],I,J),] <- test3[[i]]$beta
             Fit2[indx(xyblack_list[[i]],I,J),] <- test3[[i]]$fit
        }

        # parallel update white pixels
        test3 <- mclapply(xywhite_list, FUN = callSpatEnet, I=I, J=J, x=x, Response=Response,
                        Theta=Theta2, lambda=lambda, sv=sv,firstpen=firstpen)

        # update Theta2-values
        for(i in 1:length(xywhite_list))
        {
             Theta2[indx(xywhite_list[[i]],I,J),] <- test3[[i]]$beta
             Fit2[indx(xywhite_list[[i]],I,J),] <- test3[[i]]$fit
        }
        
        Theta_delta <- Theta_old - Theta2
        delta <- sum(abs(Theta_delta),na.rm=TRUE) / sum(abs(Theta_old),na.rm=TRUE)
        count <- count + 1
        
        #print(count)
        #print(delta)
                      
    }
    
    #name <- paste("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatFit_scan03_20.Rdata")
    #save(Theta2,Fit2,delta,file=name)  
          
    
 #   # TODO: wieder raus
 #   load("~/Dateien/Projekte/Jan-Gertheis/Code_spatEnet/results/spatFit_331.Rdata")
#    Theta2 <-  spatFit$SpatialTheta
#    Fit2 <-    spatFit$SpatialFit
#    sv <- 0.2
#    #x <- scaledX
    
    # --------------------------------------------------------------------------
    # choose reduced, non-vanishing basis functions per voxel
    # --------------------------------------------------------------------------
    Selected_array <- vector("list", N)
        
    for(n in 1:N)               
    {        
        #selected <- which(Theta2[n,]>10^(-5),arr.ind=TRUE) 
        #selected <- which(Theta2[n,]>sv/Q,arr.ind=TRUE)
        #selected <- which(Theta2[n,]>sum(Theta2[n,])/Q,arr.ind=TRUE)
        
        
        if(is.na(Theta2[n,1]))
        {
            Selected_array[[n]] <- NA
        }
        
        if(!is.na(Theta2[n,1]))
        {
        
            #selected <- which(Theta2[n,]>max(Theta2[n,])/50,arr.ind=TRUE)
            selected <- which(Theta2[n,]>max(Theta2[n,])/10,arr.ind=TRUE)
            
            if(length(selected)==0)
            {
                selected <- which.max(Theta2[n,])
            }     
            #print("selected")
            #print(selected)
            
            # treat first basis function seperately
            first <- c()
            selected_rest <- selected
            if(selected[1]==1)
            {
                first <- 1
                selected_rest <- selected[-1]
            }
             
            selected_1 <- c(Q,selected_rest[-length(selected_rest)]) # right shift
            selected_2 <- c(selected_rest[-1],0) # left shift
             
            left_bounds <- c(first,length(first)+which(selected_rest!=selected_1+1))
            right_bounds <- c(first,length(first)+which(selected_rest!=selected_2-1))
            
            #print("bounds")
            #print(left_bounds)
            #print(right_bounds)
            #print(n)
            #print(selected)
             
            selected_max <- c() 
            for(k in 1:length(left_bounds))
            {
                index_max <- which.max(Theta2[n,selected[left_bounds[k]:right_bounds[k]]])   # index within k-th group
                index_max <- left_bounds[k] + index_max -1     # index within selected
                selected_max <- c(selected_max,selected[index_max])                
            }
                            
            Selected_array[[n]] <- selected_max  
        }    
    } 

    # --------------------------------------------------------------------------
    # reestimate with reduced basis
    # --------------------------------------------------------------------------
    reducedTheta <-  Theta2
    
    for(n in 1:N)
    {
        # nonzero only for selected coefficients
        reducedTheta[n,-Selected_array[[n]]]  <- 0
    }    
                
    reducedFit <- Fit2
        
    min_it <- 80   
    max_it <- 200 
    delta <- 100
    count2 <- 0
    s <- 0
        
    lambda2 <- 10^(-10)  # relax smoothing
    sv2 <- 100    # unsrestricted reestimation
    
    #while(delta > (3/Q) | s < min_it )
    while(delta > 0.0001 & count2 < max_it | count2 < min_it)
    {    
        Theta_old <- reducedTheta
        
        # alle schwarzen Pixel parallel updaten
        
        test3 <- mclapply(xyblack_list, FUN = callSpatEnetSel, I=I, J=J, x=x, Response=Response,
                        Theta=reducedTheta, lambda=lambda2, sv=sv2,firstpen=TRUE,Selected=Selected_array)
        # aktualisiere Werte Theta2
        for(i in 1:length(xyblack_list))
        {
             reducedTheta[indx(xyblack_list[[i]],I,J),Selected_array[[indx(xyblack_list[[i]],I,J)]]] <- test3[[i]]$beta
             reducedFit[indx(xyblack_list[[i]],I,J),] <- test3[[i]]$fit             
        }
    
        # alle weissen Pixel parallel updaten
        test3 <- mclapply(xywhite_list, FUN = callSpatEnetSel, I=I, J=J, x=x, Response=Response,
                        Theta=reducedTheta, lambda=lambda2, sv=sv2,firstpen=TRUE,Selected=Selected_array)
    
        # aktualisiere Werte Theta2
        for(i in 1:length(xywhite_list))
        {
             reducedTheta[indx(xywhite_list[[i]],I,J),Selected_array[[indx(xywhite_list[[i]],I,J)]]] <- test3[[i]]$beta
             reducedFit[indx(xywhite_list[[i]],I,J),] <- test3[[i]]$fit            
        }  
        
        Theta_delta <- Theta_old - reducedTheta
        delta <- sum(abs(Theta_delta),na.rm=TRUE) / sum(abs(Theta_old),na.rm=TRUE)
        
#        if(delta==0)
#        {
#            print("delta is zero")
#            print(count2)
#            #print(Theta_old[5000:5600])
#            #print(reducedTheta[5000:5600])
#            warnings()
#        }
        
        count2 <- count2 + 1        
        s <- s+1                
    }

    return(list("SpatialTheta"=Theta2,"SpatialFit"=Fit2, "reducedTheta"=reducedTheta,"reducedFit"=reducedFit,"selected"=Selected_array,"numberIter"=c(count,count2),"delta"=delta ))
}