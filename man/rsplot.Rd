\name{rsplot}
\alias{rsplot}
\alias{rsplot.default}
\concept{RS Bag}
\title{Ray Shooting depth Bag}
\description{
Computes the Ray Shooting depth ISO of a point with respect to a bivariate data set.
}
\usage{
rsplot(x, y=NULL, factorsecondbag=2,mring=T,...)
}
\arguments{
  \item{x}{The data as a matrix, data frame or list. If it is a matrix or data frame, then each row is viewed as one bivariate observation. If it is a list, all components must be numerical vectors of equal length (coordinates of observations).}
  \item{y}{The data as a matrix, data frame or list. If it is a matrix or data frame, then each row is viewed as one bivariate observation. If it is a list, all components must be numerical vectors of equal length (coordinates of observations).}
  \item{factorsecondbag}{Factor for second bag set to 2 by default.}
  \item{mring}{Boolean value set to TRUE by default.}
  \item{...}{For future use.}
}
\details{ In dimension 2, calculates ray shooting depth of a given point with respect to the point set. Time complexity of the simple algorithms implemented is O(n log n). ISO

}
\value{
Returns with respect to data set \code{pt}, the number of line segments interested by a ray from , minimum over all rays. ISO
}
\references{N. Mustafa, S. Ray, and M. Shabbir, Statistical Data Depth of Pointsets in the Plane,\emph{in prep.}.
}
\author{Nabil Mustafa, Saurabh Ray, and Mudassir Shabbir.}

\seealso{\code{\link{rsdepth}}}
\examples{## calculation of RS depth
z = matrix(rnorm(24),nc=2)
x = rsplot(z)
}
\keyword{bivariateISO}
