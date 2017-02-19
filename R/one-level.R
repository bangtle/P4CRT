#' @title power1L
#' @description This function returns the power for a single-level randomized trials designs.
#' @param n Number of subjects in level 1 (individual level).
#' @param delta Standardized effect size.
#' @param numcovL1 Number of covariate(s) at level 1.
#' @param R2L1 R-squared at level 1.
#' @param alpha Probability of type I error.
#' @seealso \code{\link{PowerUpR}}
#' @return NULL
#' @examples
#' power1L(n=100,delta=0.2,numcovL1=1,R2L1=0.2,alpha=0.05)




power1L <- function(n,delta,numcovL1=0,R2L1=0,alpha=0.05){

  if((n<=(2+numcovL1))|!(n%%1==0)) stop("Sample n must be an integer greater than number of L1 covariate(s) + 2!")
  else if (delta<0) stop("delta cannot be negative!")
  else if ((R2L1<0)|(R2L1>=1)) stop("The R-squared must be positive and less than 1!")
  else if ((alpha<0)|(alpha>1)) stop("alpha must be between 0 and 1!")
  else if ((numcovL1<0)|!(numcovL1%%1==0)) stop("number of level-1 covariates must be a nonnegative integer!")
  else if ((numcovL1==0) & !(R2L1==0)) stop("numcovL1 = 0, so R2L1 must also be zero!")
  else {
    denom.df <- n - 2 - numcovL1
    lambda <- n*delta^2/(4*(1-R2L1))
    power <- 1-pf(qf(1-alpha,1,denom.df),1,denom.df,lambda)
  }
  result <- data.frame(n,delta,numcovL1,R2L1,alpha,power)
  return(result)
}
