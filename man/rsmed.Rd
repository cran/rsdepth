\name{rsmed}
\alias{rsmed}
\alias{rsmed.default}
\title{Bivariate RS median}
\description{
Computes the Ray Shooting median of a bivariate data set.
}
\usage{
rsmed(pt, eps=c(0),...)
}
\arguments{
  \item{pt}{The data as a matrix, data frame or list. If it is a matrix or data frame, then each row is viewed as one multivariate observation. If it is a list, all components must be numerical vectors of equal length (coordinates of observations). Only 2-dimensional data is acceptable in this version.}
  \item{eps}{eps is an optional parameter used for approximating a median in case of large data sets. It takes valye of a real between 0 and 1 and is by default set to 0 which means no approximation is used if eps is not given.}
  \item{...}{Reserved for future use.}
}
\details{
Finds out the an arbitrary point among the Ray Shooting median set of given point set. Current uses the brute-force algorithm on all O(n^4) possible points 
in the arrangment of all possible lines in complete graph on \code{pt}. For each point O(n log n) is used to find out depth so overall complexity of this algorithm is O(n^5 log n).
When approximation parameter is provided then algorithm tries to approximate by finding median of a uniform sample subset of pt of size 1/eps^2*log(1/eps). If this constant is more than the size of original set then eps value is ignored and exact median is calculate on original point set.
}
\value{A point in two dimension is returned as a single row two column vector}

\references{N. Mustafa, S. Ray, and M. Shabbir, Statistical Data Depth of Pointsets in the Plane,\emph{in prep.}.
}
\author{Nabil Mustafa, Saurabh Ray, and Mudassir Shabbir.}

\seealso{\code{\link{rsdepth}} for depth function}
\examples{
##  RS median of a two-dimensional data set
set.seed(617)
zz <- matrix(rnorm(120), nc = 2)
rsmed(zz,eps=0.2)

}
\keyword{bivariate}
\keyword{multivariate}
\keyword{nonparametric}
\keyword{robust}
