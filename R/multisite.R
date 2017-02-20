#' @title powerMS
#' @description This function returns the power for multi-site (blocked) trials designs with treatment at level 2 (site level). It addresses the following models:
#' \itemize{
#' \item Main treatment effects (assuming random site effects) with no additional covariates.
#' \item Main treatment effects (assuming random site effects) with covariate(s) at level 1.
#' \item Main treatment effects (assuming fixed site effects) with no additional covariates.
#' \item Main treatment effects (assuming fixed site effects) with covariate(s) at level 1.
#' }
#' @param n Number of subjects in level 1 (individual level).
#' @param J Number of units in level 2 (site level).
#' @param delta Standardized effect size.
#' @param sigmad Variance in treatment effects across sites.
#' @param rhoL2 Intraclass correlation (ICC) at level 2.
#' @param numcovL1 Number of covariate(s) at level 1. Default value is \code{0}.
#' @param R2L1 R-squared at level 1. Default value is \code{0}.
#' @param randeff Random site effects. Default value is \code{FALSE}.
#' @param alpha Probability of type I error. Default value is \code{0.05}.
#' @seealso \code{\link{PowerUpR}}
#' @examples
#' # All arguments are specified
#' powerMS(n=100,J=40,delta=0.2,sigmad=0.2,rhoL2=0.3,numcovL1=1,R2L1=0.2,randeff=FALSE,alpha=0.05)
#'
#' # Unspecified arguments take their default values
#' powerMS(n=100,J=40,delta=0.2,sigmad=0.2,rhoL2=0.3)

powerMS <- function(n,J,delta,sigmad,rhoL2,numcovL1=0,R2L1=0,randeff=TRUE,alpha=0.05)
{
  if ((n<=0)|!(n%%1==0)) stop("n must be a positive integer!")
  else if((J<=0)|!(J%%1==0)) stop("J must be a positive integer!")
  else if (delta<0) stop("delta cannot be negative!")
  else if (sigmad<0) stop("sigmad cannot be negative!")
  else if ((rhoL2<0)|(rhoL2>1)) stop("rhoL2 must be between 0 and 1!")
  else if ((R2L1<0)|(R2L1>=1)) stop("L1 R-squared must be positive and less than 1!")
  else if ((numcovL1<0)|!(numcovL1%%1==0)) stop("number of level-1 covariate(s) must be a nonnegative integer!")
  else if ((numcovL1==0) & !(R2L1==0)) stop("numcovL1 = 0, so R2L1 must also be zero!")
  else if ((alpha<0)|(alpha>1)) stop("alpha must be between 0 and 1!")

  else if ((numcovL1==0)&(randeff==FALSE)){ # model (5.1.1)
    num.df <- 1
    denom.df <- J-1
    lambda <- J*delta^2/(sigmad^2+4*(1-rhoL2)/n)
  }
  else if (!(numcovL1==0)&!(randeff==FALSE)){ # model (5.1.2)
    num.df <- 1
    denom.df <- J-1
    lambda <-  J*delta^2/(sigmad^2+4*(1-R2L1)*(1-rhoL2)/n)
  }
  else if ((numcovL1==0)&(randeff==TRUE)){ # model (5.2.1)
    num.df <- J-1
    denom.df <- J*(n-2)
    lambda <-  J*delta^2/(4*(1-rhoL2)/n)
  }
  else if (!(numcovL1==0)&(randeff==TRUE)){ # model (5.2.2)
    num.df <- J-1
    denom.df <- J*(n-2)-numcovL1
    lambda <- J*delta^2/(4*(1-R2L2)*(1-rhoL2)/n)
  }
  else stop("Unrecognized model. Please specify a valid model!")

  power <- 1-pf(qf(1-alpha,num.df,denom.df),num.df,denom.df,lambda)

  result <- data.frame(n,J,delta,sigmad,rhoL2,numcovL1,R2L1,randeff,alpha,power)
  return(result)
}
