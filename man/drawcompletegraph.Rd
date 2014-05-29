\name{drawcompletegraph}
\alias{drawcompletegraph}
\alias{drawcompletegraph.default}
\concept{Draws Complete Graph of a pointset}
\title{Draws Complete Graph of a pointset}
\description{
Draws Complete Graph of a pointset in plane.
}
\usage{
drawcompletegraph(x, y=NULL,startcanvas=TRUE,...)
}
\arguments{
  \item{x}{The data as a matrix, data frame or list. If it is a matrix or data frame, then each row is viewed as one bivariate observation. If it is a list, all components must be numerical vectors of equal length (coordinates of observations).}
  \item{y}{The data as a matrix, data frame or list. If it is a matrix or data frame, then each row is viewed as one bivariate observation. If it is a list, all components must be numerical vectors of equal length (coordinates of observations).} 
  \item{startcanvas}{A boolean value to let the function whether there is already a plot that we want to use or create a new canvas. Be default set to TRUE.} 
  \item{...}{For future use.}
}
\details{ 

In dimension 2, draws complete graph on a pointset.

}
\value{
Returns nothing.
}
\references{N. Mustafa, S. Ray, and M. Shabbir, Statistical Data Depth of Pointsets in the Plane,\emph{in prep.}.
}
\author{Nabil Mustafa, Saurabh Ray, and Mudassir Shabbir.}

\seealso{\code{\link{inflate}}}
\examples{## calculation of centroid of a random pointset
z = matrix(rnorm(24),nc=2)
x = drawcompletegraph(z)
}
\keyword{bivariateConvexhull}
