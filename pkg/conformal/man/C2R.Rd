\name{C2R}
\alias{C2R}
\alias{listxy2R}
\alias{matrix2R}
\title{Real/complex manipulation}
\description{
Convert objects to a real vector representation of a complex vector
}
\usage{
C2R(z)
listxy2R(lxy)
matrix2R(x)
}
\arguments{
  \item{z}{A complex vector}
  \item{lxy}{A list having two components representing real and imaginary parts}
  \item{x}{A matrix having two columns representing real and imaginary parts}
}
\value{
  A real vector whose odd elements are the real part and whose even elements are 
  the imaginary part of the complex vector
}
\author{Nick Ellis, nick.ellis@csiro.au}
\seealso{\code{\link{R2C}}, \code{\link{listxy2C}}, \code{\link{matrix2C}}}
\examples{
re <- 1:10
im <- 10:1
z <- complex(real=re,imag=im)
lis <- list(re,im)
mat <- cbind(re,im)
x1 <- C2R(z)
x2 <- listxy2R(lis)
x3 <- matrix2R(mat)
all(x1==x2)
all(x1==x3)
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
