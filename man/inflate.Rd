\name{inflate}
\alias{inflate}
\alias{inflate.default}
\concept{inflates a convex polygon}
\title{inflates a convex polygon}
\description{
Inflates a convex polygon
}
\usage{
inflate(x, y=NULL, factor=2, ...)
}
\arguments{
  \item{x}{The data as a matrix, data frame or list. If it is a matrix or data frame, then each row is viewed as one bivariate observation. If it is a list, all components must be numerical vectors of equal length (coordinates of observations).}
  \item{y}{The data as a matrix, data frame or list. If it is a matrix or data frame, then each row is viewed as one bivariate observation. If it is a list, all components must be numerical vectors of equal length (coordinates of observations).} 
  \item{factor}{An integer by default set to 2.} 
  \item{...}{For future use.}
}
\details{ 

In dimension 2, inflates a convex polygon

}
\value{
Returns nothing.
}
\references{N. Mustafa, S. Ray, and M. Shabbir, Statistical Data Depth of Pointsets in the Plane,\emph{in prep.}.
}
\author{Nabil Mustafa, Saurabh Ray, and Mudassir Shabbir.}

\seealso{\code{\link{convexhull}}}
\examples{## calculation of centroid of a random pointset
z = matrix(rnorm(24),nc=2)
x = convexhull(z)
y= inflate(x)
}
\keyword{bivariateConvexhull}
