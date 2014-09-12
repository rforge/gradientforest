\name{R2C}
\alias{R2C}
\alias{listxy2C}
\alias{matrix2C}
\title{Real/complex manipulation}
\description{
Convert objects to a complex vector
}
\usage{
R2C(ww)
listxy2C(x)
matrix2C(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ww}{A real vector representing alternating real/imaginary pairs or a complex vector}
  \item{x}{A list having two components or a matrix having two columns representing the real and imaginary parts of a complex vector}
}
\value{
  a complex vector
}
\author{Nick Ellis, nick.ellis@csiro.au}
\seealso{\code{\link{C2R}}, \code{\link{listxy2R}}, \code{\link{matrix2R}}}
\examples{
x <- 1:10
y <- 10:1
mat <- cbind(x,y)
lis <- list(x,y)
xy <- as.vector(t(mat))
z1 <- R2C(xy)
z2 <- listxy2C(lis)
z3 <- matrix2C(mat)
all(z1==z2)
all(z1==z3)
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
