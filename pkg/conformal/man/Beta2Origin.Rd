\name{Beta2Origin}
\alias{Beta2Origin}
\alias{Origin2Beta}
\alias{RBeta2Origin}
\alias{ROrigin2Beta}
\alias{matrixBeta2Origin}
\alias{matrixOrigin2Beta}
\title{Moebius transformation}
\description{
Moebius transformations to map the unit circle to the unit circle and a point
\code{beta} to the origin, or vice versa.
}
\usage{
Beta2Origin(z, beta)
Origin2Beta(z, beta)
RBeta2Origin(x, x0)
ROrigin2Beta(x, x0)
matrixBeta2Origin(x, x0)
matrixOrigin2Beta(x, x0)
}
\arguments{
  \item{z}{a complex vector of points to transform.}
  \item{beta}{a complex vector of length 1. The (pre)image of the origin.}
  \item{x}{real vector or matrix representation of \code{z}.}
  \item{x0}{real vector representation of \code{beta}.}
}
\details{
The Moebius transform \eqn{z -> (beta+z)/(1+Conj(beta)z)} maps the unit disk to 
itself and the point \code{beta} to the origin. This is implemented in \code{Beta2Origin};
the reverse transformation is \code{Origin2Beta}. In the versions \code{RBeta2Origin} 
and \code{ROrigin2Beta}, the input and output use the real representation of complex values (see \code{\link{C2R}}).
In the versions \code{matrixBeta2Origin} and \code{matrixOrigin2Beta}, the input 
and output use the two-column matrix representation.
}
\value{
  For \code{Beta2Origin} and \code{Origin2Beta} the transformed complex vector. 
  For \code{RBeta2Origin} and \code{ROrigin2Beta} the real representation of the 
  transformed complex vector. For \code{matrixBeta2Origin} and \code{matrixOrigin2Beta}, 
  the two-column matrix representation of the transformed complex vector.
}
\author{Nick Ellis, nick.ellis@csiro.au}
\seealso{\code{\link{C2R}}}
\examples{
CircleGrid <- function(nradial) { # create a random pattern to fill a disk
  r <- seq(0,1,len=nradial)
  nth <- pmax(1,ceiling(2*pi*r*(nradial-1)))
  th0 <- runif(nth,-pi/nth,pi/nth)
  list(x=unlist(lapply(1:nradial,function(i,r,nth,th0) r[i]*cos(th0[i]+seq(0,2*pi,len=nth[i]+1)[-nth[i]-1]),r=r,nth=nth,th0=th0)),
       y=unlist(lapply(1:nradial,function(i,r,nth,th0) r[i]*sin(th0[i]+seq(0,2*pi,len=nth[i]+1)[-nth[i]-1]),r=r,nth=nth,th0=th0)))
}
cDisk <- R2C(listxy2R(CircleGrid(30)))
centre <- R2C(c(0.5,-0.1))
#centre <- R2C(listxy2R(locator(1)))
cPred <- Origin2Beta(R2C(disk), centre)
cDisk <- R2C(disk)
cols <- hsv(h=Arg(cDisk)/2/pi+0.5,s=pmin(1,Mod(cDisk)))
par(mfrow=c(1,2))
plot(cDisk, col=cols, asp=1, pch=16, main="Colour wheel")
plot(cDisk, type='n', asp=1, main="Moebius transformation")
points(cPred, col=cols,pch=16)
points(centre, pch=4)
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
