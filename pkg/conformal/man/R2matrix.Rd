\name{R2matrix}
\alias{R2matrix}
\alias{C2matrix}
\alias{R2listxy}
\alias{C2listxy}
\title{Real/complex manipulation}
\description{
Convert objects to matrix or list representation of a complex vector
}
\usage{
R2matrix(x)
C2matrix(x)
R2listxy(z)
C2listxy(z)
}
\arguments{
  \item{z,x}{A complex vector or real representation of one}
}
\value{
  For \code{R2matrix} or \code{R2matrix}, a matrix having two columns representing the real and imaginary parts of the complex vector. 
  For \code{R2listxy} or \code{C2listxy}, a list having two components representing the real and imaginary parts of the complex vector. 
}
\author{Nick Ellis, nick.ellis@csiro.au}
\seealso{\code{\link{C2R}}, \code{\link{listxy2R}}, \code{\link{matrix2R}}}
\examples{
re <- 1:10
im <- 10:1
z <- complex(real=re,imag=im)
x <- C2R(z)
mat1 <- C2matrix(z)
mat2 <- R2matrix(x)
lis1 <- C2listxy(z)
lis2 <- R2listxy(x)
all.equal(mat1,mat2)
all.equal(lis1,lis2)
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
