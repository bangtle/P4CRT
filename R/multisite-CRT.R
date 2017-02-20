#' @title powerMSCRT
#' @description This function returns the power for a multi-site cluster randomized trials designs with treatment at level 2.
#' @param n Number of subjects in level 1 (individual level).
#' @param J Number of units in level 2 (cluster level).
#' @param K Number of units in level 3 (site level).
#' @param delta Standardized effect size.
#' @param sigmad Variance in treatment effects across sites.
#' @param rhoL2 Intraclass correlation (ICC) at level 2.
#' @param rhoL3 Intraclass correlation (ICC) at level 3.
#' @param numcovL1 Number of covariate(s) at level 1.
#' @param numcovL2 Number of covariate(s) at level 2.
#' @param R2L1 R-squared at level 1.
#' @param R2L2 R-squared at level 2.
#' @param randeff Existence of random site effects.
#' @param alpha Probability of type I error.
#' @seealso \code{\link{PowerUpR}}
#' @examples
#' powerMSCRT(n=100,J=40,K=40,delta=0.2,rhoL2=0.3,rhoL3=0.2,numcovL1=0,numcovL2=0,R2L1=0,R2L2=0,randeff=TRUE,alpha=0.05)
#'

powerMSCRT <- function(n,J,K,delta,sigmad,rhoL2,rhoL3,numcovL1=0,numcovL2=0,R2L1=0,R2L2=0,randeff=TRUE,alpha=0.05)
{
  if ((n<=0)|!(n%%1==0)) stop("n must be a positive integer!")
  else if((J<=0)|!(J%%1==0)) stop("J must be a positive integer!")
  else if((K<=0)|!(K%%1==0)) stop("K must be a positive integer!")
  else if (delta<0) stop("delta cannot be negative!")
  else if (sigmad<0) stop("sigmad cannot be negative!")
  else if ((rhoL2<0)|(rhoL2>1)) stop("rhoL2 must be between 0 and 1!")
  else if ((rhoL3<0)|(rhoL3>1)) stop("rhoL3 must be between 0 and 1!")
  else if ((R2L1<0)|(R2L1>=1)) stop("L1 R-squared must be positive and less than 1!")
  else if ((R2L2<0)|(R2L2>=1)) stop("L2 R-squared must be positive and less than 1!")
  else if ((numcovL1<0)|!(numcovL1%%1==0)) stop("number of level-1 covariate(s) must be a nonnegative integer!")
  else if ((R2L2<0)|(R2L2>=1)) stop("L2 R-squared must be positive and less than 1!")
  else if ((numcovL2<0)|!(numcovL2%%1==0)) stop("number of level-2 covariate(s) must be a nonnegative integer!")
  else if ((numcovL2==0) & !(R2L2==0)) stop("numcovL2 = 0, so R2L2 must also be zero!")
  else if (!(randeff==TRUE)|!(randeff==T)&!(randeff==FALSE)|!(randeff==F)) stop("randeff must be either TRUE (or T) or FALSE (or F) only!")
  else if ((alpha<0)|(alpha>1)) stop("alpha must be between 0 and 1!")

  else if ((numcovL1==0)&(numcovL2=0)&(randeff==TRUE)){ # model (6.1.1)
    num.df <- 1
    denom.df <- K-1
    lambda <- K*delta^2/(sigmad^2+4*(rhoL2+(1-rhoL3-rhoL2)/n)/J)
  }
  else if (!(numcovL1==0)&(numcovL2=0)&(randeff==TRUE)){ # model (6.1.2)
    num.df <- 1
    denom.df <- K-1
    lambda <- K*delta^2/(sigmad^2+4*(rhoL2+(1-R2L1)*(1-rhoL3-rhoL2)/n)/J)
  }
  else if ((numcovL1==0)&!(numcovL2=0)&(randeff==TRUE)){ # model (6.1.3)
    num.df <- 1
    denom.df <- K-1-numcovL2
    lambda <- K*delta^2/(sigmad^2+4*((1-R2L2)*rhoL2+(1-rhoL3-rhoL2)/n)/J)
  }
  else if (!(numcovL1==0)&!(numcovL2=0)&(randeff==TRUE)){ # model (6.1.4)
    num.df <- 1
    denom.df <- K-1-numcovL2
    lambda <- K*delta^2/(sigmad^2+4*((1-R2L2)*rhoL2+(1-R2L1)*(1-rhoL3-rhoL2)/n)/J)
  }
  else if ((numcovL1==0)&(numcovL2=0)&(randeff==FALSE)){ # model (6.2.1)
    num.df <- 1
    denom.df <- K-1
    lambda <- K*J*delta^2/(4*(rhoL2+(1-rhoL3-rhoL2)/n)/J)
  }
  else if (!(numcovL1==0)&(numcovL2=0)&(randeff==FALSE)){ # model (6.2.2)
    num.df <- 1
    denom.df <- K-1
    lambda <- K*delta^2/(4*(rhoL2+(1-R2L1)*(1-rhoL3-rhoL2)/n)/J)
  }
  else if ((numcovL1==0)&!(numcovL2=0)&(randeff==FALSE)){ # model (6.2.3)
    num.df <- 1
    denom.df <- K-1-numcovL2
    lambda <- K*J*delta^2/(4*((1-R2L2)*rhoL2+(1-rhoL3-rhoL2)/n)/J)
  }
  else if (!(numcovL1==0)&!(numcovL2=0)&(randeff==FALSE)){ # model (6.2.4)
    num.df <- 1
    denom.df <- K-1-numcovL2
    lambda <- K*J*delta^2/(4*((1-R2L2)*rhoL2+(1-R2L1)*(1-rhoL3-rhoL2)/n)/J)
      }
  else stop("Unrecognized model. Please specify a valid model!")

  power <- 1-pf(qf(1-alpha,num.df,denom.df),num.df,denom.df,lambda)

  result <- data.frame(n,J,delta,sigmad,rhoL2,numcovL1,R2L1,randeff,alpha,power)
  return(result)
}
