\name{centroid}
\alias{centroid}
\alias{centroid.default}
\concept{Centroid of a convex polygon}
\title{Centroid of a convex polygon}
\description{
Computes Centroid of a convex polygon in plane.
}
\usage{
centroid(x, y=NULL,...)
}
\arguments{
  \item{x}{The data as a matrix, data frame or list. If it is a matrix or data frame, then each row is viewed as one bivariate observation. If it is a list, all components must be numerical vectors of equal length (coordinates of observations).}
  \item{y}{The data as a matrix, data frame or list. If it is a matrix or data frame, then each row is viewed as one bivariate observation. If it is a list, all components must be numerical vectors of equal length (coordinates of observations).} 
  \item{...}{For future use.}
}
\details{ In dimension 2, calculates centroid of a convex polygon.

}
\value{
Returns with respect to data set, the centroid point in plane.
}
\references{N. Mustafa, S. Ray, and M. Shabbir, Statistical Data Depth of Pointsets in the Plane,\emph{in prep.}.
}
\author{Nabil Mustafa, Saurabh Ray, and Mudassir Shabbir.}

\seealso{\code{\link{inflate}}}
\examples{## calculation of centroid of a random pointset
z = matrix(rnorm(24),nc=2)
x = centroid(z)
}
\keyword{bivariateCentroid}
