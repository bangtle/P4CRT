#' @title power3L
#' @description This function returns the power for a three-level cluster randomized trials designs.
#' @param n Number of subjects in level 1 (individual level).
#' @param J Number of units in level 2 (subcluster level).
#' @param K Number of units in level 3 (cluster level).
#' @param delta Standardized effect size.
#' @param rhoL2 Intraclass correlation (ICC) at level 2.
#' @param rhoL3 Intraclass correlation (ICC) at level 3.
#' @param numcovL1 Number of covariate(s) at level 1. Default value is \code{0}.
#' @param numcovL2 Number of covariate(s) at level 2. Default value is \code{0}.
#' @param numcovL3 Number of covariate(s) at level 3. Default value is \code{0}.
#' @param R2L1 R-squared at level 1. Default value is \code{0}.
#' @param R2L2 R-squared at level 2. Default value is \code{0}.
#' @param R2L3 R-squared at level 3. Default value is \code{0}.
#' @param modL1 Existence of moderator at level 1. Default value is \code{FALSE}.
#' @param modL2 Existence of moderator at level 2. Default value is \code{FALSE}.
#' @param modL3 Existence of moderator at level 3. Default value is \code{FALSE}.
#' @param alpha Probability of type I error. Default value is \code{0.05}.
#' @seealso \code{\link{PowerUpR}}
#' @examples
#' # All arguments are specified
#' power3L(n=100,J=40,K=40,delta=0.2,rhoL2=0.3,rhoL3=0.2,numcovL1=1,numcovL2=1,numcovL3=1,R2L1=0.2,R2L2=0.2,modL1=FALSE,modL2=FALSE,modL3=FALSE,alpha=0.05)
#'
#'# Unspecified arguments take their default values
#' power3L(n=100,J=40,K=40,delta=0.2,rhoL2=0.3,rhoL3=0.2)


power3L <- function(n,J,K,delta,rhoL2,rhoL3,numcovL1=0,numcovL2=0,numcovL3=0,R2L1=0,R2L2=0,R2L3=0,modL1=FALSE,modL2=FALSE,modL3=FALSE,alpha=0.05){

  if ((n<=0)|!(n%%1==0)) stop("n must be a positive integer!")
  else if((J<=0)|!(J%%1==0)) stop("J must be a positive integer!")
  else if ((K<=0)|!(K%%1==0)) stop("K must be a positive integer!")
  #else if ((modL1==TRUE)&(J<=4)) stop("L2 moderator exists, so J must be greater than 4!")
  else if (delta<0) stop("delta cannot be negative!")
  else if ((rhoL2<0)|(rhoL2>1)) stop("rhoL2 must be between 0 and 1!")
  else if ((rhoL3<0)|(rhoL3>1)) stop("rhoL3 must be between 0 and 1!")
  else if ((R2L1<0)|(R2L1>=1)) stop("L1 R-squared must be positive and less than 1!")
  else if ((R2L2<0)|(R2L2>=1)) stop("L2 R-squared must be positive and less than 1!")
  else if ((R2L3<0)|(R2L3>=1)) stop("L3 R-squared must be positive and less than 1!")
  else if ((numcovL1<0)|!(numcovL1%%1==0)) stop("number of level-1 covariate(s) must be a nonnegative integer!")
  else if ((numcovL1==0) & !(R2L1==0)) stop("numcovL1 = 0, so R2L1 must also be zero!")
  else if ((numcovL2<0)|!(numcovL2%%1==0)) stop("number of level-2 covariate(s) must be a nonnegative integer!")
  else if ((numcovL2==0) & !(R2L2==0)) stop("numcovL2 = 0, so R2L2 must also be zero!")
  else if ((numcovL3<0)|!(numcovL3%%1==0)) stop("number of level-3 covariate(s) must be a nonnegative integer!")
  else if ((numcovL3==0) & !(R2L3==0)) stop("numcovL3 = 0, so R2L3 must also be zero!")
  else if ((alpha<0)|(alpha>1)) stop("alpha must be between 0 and 1!")

  else if ((numcovL1==0)&(numcovL2==0)&(numcovL3==0)&(modL1==FALSE)&(modL2==FALSE)&(modL3==FALSE)){ # model (4.1.1)
    denom.df <- K-2
    lambda <- delta^2/(4*(rhoL3+(rhoL2+(1-rhoL3-rhoL2)/n)/J)/K)
  }
  else if (!(numcovL1==0)&(numcovL2==0)&(numcovL3==0)&(modL1==FALSE)&(modL2==FALSE)&(modL3==FALSE)){ # model (4.1.2)
    denom.df <- K-2
    lambda <- delta^2/(4*(rhoL3+(rhoL2+(1-R2L1)*(1-rhoL3-rhoL2)/n)/J)/K)
  }
  else if ((numcovL1==0)&!(numcovL2==0)&(numcovL3==0)&(modL1==FALSE)&(modL2==FALSE)&(modL3==FALSE)){ # model (4.1.3)
    denom.df <- K-2
    lambda <- delta^2/(4*(rhoL3+((1-R2L2)*rhoL2+(1-rhoL3-rhoL2)/n)/J)/K)
  }
  else if ((numcovL1==0)&(numcovL2==0)&!(numcovL3==0)&(modL1==FALSE)&(modL2==FALSE)&(modL3==FALSE)){ # model (4.1.4)
    denom.df <- K-2-numcovL3
    lambda <- delta^2/(4*((1-R2L3)*rhoL3+(rhoL2+(1-rhoL3-rhoL2)/n)/J)/K)
  }
  else if (!(numcovL1==0)&!(numcovL2==0)&!(numcovL3==0)&(modL1==FALSE)&(modL2==FALSE)&(modL3==FALSE)){ # model (4.1.5)
      denom.df <- K-2-numcovL3
      lambda <- delta^2/(4*((1-R2L3)*rhoL3+((1-R2L2)*rhoL2+(1-R2L1)*(1-rhoL3-rhoL2)/n)/J)/K)
  }
  else if ((numcovL1==0)&(numcovL2==0)&(numcovL3==0)&(modL1==FALSE)&(modL2==FALSE)&(modL3==TRUE)){ # model (4.2.1)
    denom.df <- K-4
    lambda <- delta^2/(16*((1-R2L3)*rhoL3+(rhoL2+(1-rhoL3-rhoL2)/n)/J)/K)
  }
  else if ((numcovL1==0)&(numcovL2==0)&!(numcovL3==0)&(modL1==FALSE)&(modL2==FALSE)&(modL3==TRUE)){ # model (4.2.2)
    denom.df <- K-4-numcovL3
    lambda <- delta^2/(16*((1-R2L3)*rhoL3+(rhoL2+(1-rhoL3-rhoL2)/n)/J)/K)
  }
  else if ((numcovL1==0)&!(numcovL2==0)&(numcovL3==0)&(modL1==FALSE)&(modL2==TRUE)&(modL3==FALSE)){ # model (4.3.1)
    denom.df <- J*K-J-2
    lambda <- delta^2/(16*(((1-R2L2)*rhoL2+(1-rhoL3-rhoL2)/n)/J)/K)
  }
  else if ((numcovL1==0)&!(numcovL2==0)&(numcovL3==0)&(modL1==FALSE)&(modL2==FALSE)&(modL3==FALSE)){ # model (4.3.2)
    denom.df <- J*K-J-2-numcovL2
    lambda <- delta^2/(16*(((1-R2L2)*rhoL2+(1-rhoL3-rhoL2)/n)/J)/K)
  }
  else if ((numcovL1==0)&(numcovL2==0)&(numcovL3==0)&(modL1==TRUE)&(modL2==FALSE)&(modL3==FALSE)){ # model (4.4.1)
    denom.df <- n*J*K-J*K-K-2
    lambda <- delta^2/(16*(1-R2L1)*(1-rhoL3-rhoL2)/(n*J*K))
  }
  else if (!(numcovL1==0)&(numcovL2==0)&(numcovL3==0)&(modL1==TRUE)&(modL2==FALSE)&(modL3==FALSE)){ # model (4.4.2)
    denom.df <- n*J*K-J*K-K-2-numcovL1
    lambda <- delta^2/(16*(1-R2L1)*(1-rhoL3-rhoL2)/(n*J*K))
  }
  else stop("Unrecognized model. Please specify a valid model!")

  power <- 1-pf(qf(1-alpha,1,denom.df),1,denom.df,lambda)

  result <- data.frame(n,J,K,delta,rhoL2,rhoL3,numcovL1,numcovL2,numcovL3,R2L1,R2L2,R2L3,modL1,modL2,modL3,alpha,power)
  return(result)
}

