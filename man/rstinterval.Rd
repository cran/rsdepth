\name{rstinterval}
\alias{rstinterval}
\alias{rstinterval.default}
\title{Bivariate RS Rings}
\description{
Computes the Ray Shooting rings of a bivariate data set.
}
\usage{
rstinterval(pt, beta=c(0.90), sampleSize=c(250), M=c(50),clr=FALSE, ...)
}
\arguments{
  \item{pt}{The data as a matrix, data frame or list. If it is a matrix or data frame, then each row is viewed as one multivariate observation. If it is a list, all components must be numerical vectors of equal length (coordinates of observations). Only 2-dimensional data is acceptable in this version.}
  \item{beta}{beta is a parameter between 0 and 1 determines the accuracy of the interval. Set to 0.90 by default.}
  \item{sampleSize}{Size of the sample data set.}
  \item{M}{Size of test data set.}
  \item{clr}{Clear the canvas before use or not. Boolean and set to FALSE by default}
  \item{...}{Reserved for future use.}
}
\details{
This function creates a two dimension generalization of confidence intervals of data. A bag that contains beta fraction of data points is constructed.
}
\value{Should not return anything}

\references{N. Mustafa, S. Ray, and M. Shabbir, Statistical Data Depth of Pointsets in the Plane,\emph{in prep.}.
}
\author{Nabil Mustafa, Saurabh Ray, and Mudassir Shabbir.}

\seealso{\code{\link{rsdepth}} for depth function}
\examples{
##  RS median of a two-dimensional data set
set.seed(617)
zz <- matrix(rnorm(600), nc = 2)
rstinterval(zz)
}
\keyword{bivariate}
\keyword{multivariate}
\keyword{nonparametric}
\keyword{robust}
