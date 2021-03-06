% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multisite.R
\name{powerMS}
\alias{powerMS}
\title{powerMS}
\usage{
powerMS(n, J, delta, sigmad, rhoL2, numcovL1 = 0, R2L1 = 0,
  randeff = TRUE, alpha = 0.05)
}
\arguments{
\item{n}{Number of subjects in level 1 (individual level).}

\item{J}{Number of units in level 2 (site level).}

\item{delta}{Standardized effect size.}

\item{sigmad}{Variance in treatment effects across sites.}

\item{rhoL2}{Intraclass correlation (ICC) at level 2.}

\item{numcovL1}{Number of covariate(s) at level 1. Default value is \code{0}.}

\item{R2L1}{R-squared at level 1. Default value is \code{0}.}

\item{randeff}{Random site effects. Default value is \code{FALSE}.}

\item{alpha}{Probability of type I error. Default value is \code{0.05}.}
}
\description{
This function returns the power for multi-site (blocked) trials designs with treatment at level 2 (site level). It addresses the following models:
\itemize{
\item Main treatment effects (assuming random site effects) with no additional covariates.
\item Main treatment effects (assuming random site effects) with covariate(s) at level 1.
\item Main treatment effects (assuming fixed site effects) with no additional covariates.
\item Main treatment effects (assuming fixed site effects) with covariate(s) at level 1.
}
}
\examples{
# All arguments are specified
powerMS(n=100,J=40,delta=0.2,sigmad=0.2,rhoL2=0.3,numcovL1=1,R2L1=0.2,randeff=FALSE,alpha=0.05)

# Unspecified arguments take their default values
powerMS(n=100,J=40,delta=0.2,sigmad=0.2,rhoL2=0.3)
}
\seealso{
\code{\link{PowerUpR}}
}
