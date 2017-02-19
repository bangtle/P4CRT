#' @title power2L
#' @description This function returns the power for a two-level cluster randomized trials designs.
#' @param n Number of units in level 1 (individual level).
#' @param J Number of units in level 2 (cluster level)
#' @param delta Standardized effect size.
#' @param rhoL2 Intraclass correlation (ICC) at level 2.
#' @param numcovL1 Number of covariate(s) at level 1.
#' @param numcovL2 Number of covariate(s) at level 2.
#' @param R2L1 R-squared at level 1.
#' @param R2L2 R-squared at level 2.
#' @param modL1 Existence of moderator at level 1.
#' @param modL2 Existence of moderator at level 2.
#' @param alpha Probability of type I error.
#' @seealso \code{\link{PowerUpR}}
#' @examples 
#' power2L(n=100,J=40,delta=0.2,rhoL2=0.3,numcovL1=1,numcovL2=0,R2L1=0.2,R2L2=0,modL1=TRUE,modL2=FALSE,alpha=0.05)



power2L <- function(n,J,delta,rhoL2,numcovL1=0,numcovL2=0,R2L1=0,R2L2=0,modL1=FALSE,modL2=FALSE,alpha=0.05){

  if ((n<=0)|!(n%%1==0)) stop("n must be a positive integer!")
  else if((J<=0)|!(J%%1==0)) stop("J must be a positive integer!")
  else if ((modL1==TRUE)&(J<=4)) stop("L2 moderator exists, so J must be greater than 4!")
  else if (delta<0) stop("delta cannot be negative!")
  else if ((rhoL2<0)|(rhoL2>1)) stop("rhoL2 must be between 0 and 1!")
  else if ((R2L1<0)|(R2L1>=1)) stop("L1 R-squared must be positive and less than 1!")
  else if ((R2L2<0)|(R2L2>=1)) stop("L2 R-squared must be positive and less than 1!")
  else if ((numcovL1<0)|!(numcovL1%%1==0)) stop("number of level-1 covariate(s) must be a nonnegative integer!")
  else if ((numcovL1==0) & !(R2L1==0)) stop("numcovL1 = 0, so R2L1 must also be zero!")
  else if ((numcovL2<0)|!(numcovL2%%1==0)) stop("number of level-2 covariate(s) must be a nonnegative integer!")
  else if ((numcovL2==0) & !(R2L2==0)) stop("numcovL2 = 0, so R2L2 must also be zero!")
  else if ((alpha<0)|(alpha>1)) stop("alpha must be between 0 and 1!")

  else if ((numcovL1==0)&(numcovL2==0)&(modL1==FALSE)&(modL2==FALSE)){   # model (3.1.1)
    denom.df <- J-2
    lambda <- delta^2/(4*(rhoL2+(1-rhoL2)/n)/J)
  }
  else if (!(numcovL1==0)&(numcovL2==0)&(modL1==FALSE)&(modL2==FALSE)){ # model (3.1.2)
    denom.df <- J-2
    lambda <- delta^2/(4*(rhoL2+(1-R2L1)*(1-rhoL2)/n)/J)
  }
  else if ((numcovL1==0)&!(numcovL2==0)&(modL1==FALSE)&(modL2==FALSE)){ # model (3.1.3)
    denom.df <- J-2-numcovL2
    lambda <-  delta^2/(4*((1-R2L2)*rhoL2+(1-rhoL2)/n)/J)
  }
  else if (!(numcovL1==0)&!(numcovL2==0)&(modL1==FALSE)&(modL2==FALSE)){ # model (3.1.4)
    denom.df <- J-2-numcovL2
    lambda <-  delta^2/(4*((1-R2L2)*rhoL2+(1-R2L1)*(1-rhoL2)/n)/J)
  }
  else if ((numcovL1==0)&(numcovL2==0)&(modL1==FALSE)&(modL2==TRUE)){ # model (3.2.1)
    denom.df <- J-4
    lambda <- delta^2/(16*((1-R2L2)*rhoL2+(1-rhoL2)/n)/J)
  }
  else if ((numcovL1==0)&(!numcovL2==0)&(modL1==FALSE)&(modL2==TRUE)){ # model (3.2.2)
    denom.df <- J-4-numcovL2
    lambda <-  delta^2/(16*((1-R2L2)*rhoL2+(1-rhoL2)/n)/J)
  }
  else if ((numcovL1==0)&(numcovL2==0)&(modL1==TRUE)&(modL2==FALSE)){ # model (3.3.1)
    denom.df <- n*J-J-2
    lambda <-  delta^2/((16*(1-R2L1)*(1-rhoL2))/(n*J))
  }
  else if (!(numcovL1==0)&(numcovL2==0)&(modL1==TRUE)&(modL2==FALSE)){ # model (3.3.2)
    denom.df <- n*J-J-2-numcovL1
    lambda <-  delta^2/((16*(1-R2L1)*(1-rhoL2))/(n*J))
  }

  else stop("Unrecognized model. Please specify a valid model!")

  power <- 1-pf(qf(1-alpha,1,denom.df),1,denom.df,lambda)

  result <- data.frame(n,J,delta,rhoL2,numcovL1,numcovL2,R2L1,R2L2,modL1,modL2,alpha,power)
  return(result)
}

