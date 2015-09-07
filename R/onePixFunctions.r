# visualization

plotimage <- function(tp,x=0,y=0,limits=NULL)
  {
    image(1:41,1:41,Images[[tp]],zlim=limits,xlab="x",ylab="y")
    abline(h=y)
    abline(v=x)
    points(x,y,cex=2)
  }

plotimages <- function(tps, wait=0.1, x=0, y=0, time=tps, limits=NULL)
  {
    for (tp in tps)
      {
        plotimage(tp,x,y,limits=limits)
        title(paste("time = ", round(time[tp],2)))
        Sys.sleep(wait)
      }
  }


# models

Cp <- function(time, model=fritz.hansen)
{
        if(model=="fritz.hansen")
        {
            D <- 1
            a1 <- 2.4
            a2 <- 0.62
            m1 <- 3.01
            m2 <- 0.016
            return(D * (a1 * exp(-m1 * time) + a2 * exp(-m2 * time)) )   
        }
        if(model=="tofts.kermode")
        {
            D <- 0.1
            a1 <- 3.99
            a2 <- 4.78
            m1 <- 0.144
            m2 <- 0.0111
            return(D * (a1 * exp(-m1 * time) + a2 * exp(-m2 * time))   ) 
        }
        else
        {
            print(model)
            print(" is not a valid model")
            return()
        }
        
        
}

model.weinmann <- function(time, th3, model=fritz.hansen)
{

       if(model=="fritz.hansen")
       {
           D <- 1
           a1 <- 2.4
           a2 <- 0.62
           m1 <- 3.01
           m2 <- 0.016
           erg <- D  * ((a1/(m1 - th3)) * (exp(-(time *
               th3)) - exp(-(time * m1))) + (a2/(m2 - th3)) *
               (exp(-(time * th3)) - exp(-(time * m2))))
           erg[time <= 0] <- 0
           return(erg)
       }   
       if(model=="tofts.kermode")
       {
           D <- 0.1
           a1 <- 3.99
           a2 <- 4.78
           m1 <- 0.144
           m2 <- 0.0111
           erg <- D  * ((a1/(m1 - th3)) * (exp(-(time *
               th3)) - exp(-(time * m1))) + (a2/(m2 - th3)) *
               (exp(-(time * th3)) - exp(-(time * m2))))
           erg[time <= 0] <- 0
           return(erg)
       }   
       else
        {
            print(model)
            print(" is not a valid model")
            return()
        }
}

oneComp <- function(time, vp, Ktrans, kep)
  {
    vp*Cp(time) + Ktrans*model.weinmann(time,kep)
  }

oneCompError <- function(theta, CT, time)
  {
    vp <- theta[1]
    Ktrans <- theta[2]
    kep <- theta[3]
    fit <- oneComp(time, vp, Ktrans, kep)
    return(sum((CT - fit)^2))
  }


# non-negative elastic net
library(quadprog)

# fit
nnenet <- function(x, y, lambda, sv, firstpen=T)
  {
    
    x <- as.matrix(x)
    
    n <- length(y)
    p <- ncol(x)
    Xy <- t(x)%*%y
    XX <- t(x)%*%x
    dvec <- Xy
    if (firstpen)
      {
        Dmat <- XX + lambda*diag(1,p)
        Amat <- t(rbind(diag(1,p),rep(-1,p)))
      }
    else
      {
        Dmat <- XX + lambda*diag(c(0,rep(1,p-1)))
        Amat <- t(rbind(diag(1,p),c(0,rep(-1,p-1))))
      }
    bvec <- c(rep(0,p),-sv)

    nnen <- solve.QP(Dmat=Dmat, dvec=dvec, Amat=Amat, bvec=bvec)
    nnen <- (1+lambda)*nnen$solution
    return(list("beta"=nnen,"fit"=x%*%nnen))
  }
  
# two-step voxelwise fit
nnenet2 <- function(x, y, lambda, sv, firstpen=T)
  {
    
    first <- nnenet(x,y,lambda,sv,firstpen=T)
        
    fit <- first$fit
    beta <- first$beta
    
    #print("after first fit")
    #print(beta)
    #print(fit)
    
    Q <- ncol(x)
    
    selected <- which(beta>max(beta)/50)
            
    if(length(selected)==0)
    {
        selected <- which.max(beta)
    }     
        
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
    #print(selected)
                
    selected_max <- c() 
    for(k in 1:length(left_bounds))
    {
        index_max <- which.max(beta[selected[left_bounds[k]:right_bounds[k]]])   # index within k-th group
        index_max <- left_bounds[k] + index_max -1     # index within selected
        selected_max <- c(selected_max,selected[index_max])                
    }
                    
    selected_final <- selected_max  
    
    #print("selected_final")
    #print(selected_final)
    
    # refit relaxed
    #lambda <- 10^(-10)
    lambda <- 10^(-10)
    sv <- 100
    
    second <-  nnenet(x[,selected_final],y,lambda,sv,firstpen=T)
    
    reducedFit <- second$fit
    reducedBeta <- second$beta
        
    return(list("beta"=beta,"fit"=fit,"reducedBeta"=reducedBeta,"reducedFit"=reducedFit,"selected"=selected_final))
  }  


 



## spatial fit
### x: Matrix of basis functions scaled such that entries have unit variance
### y: Response, Concentration
### lambda: strength of penalization (Ridge term)
### sv: penalization parameter for Lasso term
### firstpen: indicates, if first term (vp-term) is penalized
### Xi: vector of parameter values from neighboring pixels
#spatEnet <- function(x, y, lambda, sv, firstpen=T, Xi)
#  {
#    x0 <- x
#  
#    T <- length(y)   # T (length of time)
#    q <- ncol(x)     # q (number of basis functions)
#  
#    # extent observation vector with pseudo-observations
#    y <- c(y,sqrt(lambda)*Xi)
#      
#    # extent design matrix    
#    d1 <- rbind(c(-1,rep(0,q-1)))
#    D2 <- rbind(d1,d1,d1,d1)
#    D <- D2
#    for(k in 1:(q-1))
#    {
#        D <- rbind(D,D2[,c((q-k+1):q,1:(q-k))] )
#    }
#    # X = (Z, sqrt(lambda)D)'
#    x <- rbind(x,sqrt(lambda)*D)
#    
#    # dimensions of extended problem with pseudo-observations  
#    n <- length(y)
#    p <- ncol(x)
#    
#    dvec <- Xy <- t(x)%*%y
#    Dmat <- XX <- t(x)%*%x
#    Amat <- t(rbind(diag(1,q),rep(-1,q)))
#    bvec <- c(rep(0,q),-sv)
#
#    nnen <- solve.QP(Dmat=Dmat, dvec=dvec, Amat=Amat, bvec=bvec)
#    nnen <- (1+lambda)*nnen$solution
#    return(list("beta"=nnen,"fit"=x0%*%nnen))
#  }
#


# cross-validation
CVnnenet <- function(K=5, x, y, lambda, svalues, neworder=NULL, printFold=T,
firstpen = T)
  {
    n <- nrow(x)
    nK <- floor(n/K)
    if (length(neworder)==0)
      {
        neworder <- sample(1:n,n)
      }
      
    cv <- matrix(0,length(svalues),length(lambda))
    colnames(cv) <- lambda
    rownames(cv) <- svalues
    for (k in 1:K)
      {
        if (printFold)
          cat("Fold",k,"\n")
        if (k < K)
          {
            ik <- neworder[(k-1)*nK + (1:nK)]
          }
        else
          {
            ik <- neworder[((K-1)*nK+1):n]
          }
        j <- 1
        for (lam in lambda)
          {
            i <- 1
            for (sv in svalues)
              {
                result <- nnenet(y=y[-ik], x=x[-ik,], lam, sv, firstpen=firstpen)
                cv[i,j] <- cv[i,j] + sum((y[ik] - x[ik,]%*%result$beta)^2)
                i <- i+1
              }
            j <- j+1
          }
      }
    return(cv)
  }
