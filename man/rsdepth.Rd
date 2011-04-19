\name{rsdepth}
\alias{rsdepth}
\alias{rsdepth.default}
\concept{RS depth}
\title{RS Depth calculation}
\description{
Computes the Ray Shooting depth of a point with respect to a bivariate data set.
}
\usage{
rsdepth(pt,q, ...)
}
\arguments{
  \item{q}{Numerical vector whose depth is to be calculated. 
  Data needs to be 2-dimensional.}
  \item{pt}{The data as a matrix, data frame or list. If it is a matrix or data frame, then each row is viewed as one bivariate observation. If it is a list, all components must be numerical vectors of equal length (coordinates of observations).}
  \item{...}{For future use.}
}
\details{ In dimension 2, calculates ray shooting depth of a given point with respect to the point set. Time complexity of the simple algorithms implemented is O(n log n).

}
\value{
Returns the exact depth of bivariate point \code{q} with respect to data set \code{pt}, the number of line segments interested by a ray from \code{q}, minimum over all rays.
}
\references{N. Mustafa, S. Ray, and M. Shabbir, Statistical Data Depth of Pointsets in the Plane,\emph{in prep.}.
}
\author{Nabil Mustafa, Saurabh Ray, and Mudassir Shabbir.}

\seealso{\code{\link{rsmed}}}
\examples{## calculation of RS depth
z = matrix(rnorm(24),nc=2)
x = matrix(rnorm(2),nc=2)
rsdepth(z, x)

}
\keyword{bivariate}
\keyword{multivariate}
\keyword{nonparametric}
\keyword{robust}
