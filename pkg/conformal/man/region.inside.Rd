\name{region.inside}
\alias{region.inside}
\title{Test for points inside a polygon}
\description{
Given a set of vertices of a polygon, determine which test points lie inside
the polygon.
}
\usage{
region.inside(test, vertices)
}
\arguments{
  \item{test}{a two-column matrix of points to be tested.}
  \item{vertices}{a two-column matrix of vertices of the polygon. The last vertex
  should be the same as the first.}
}
\details{
  The algorithm works by checking whether the angle from the test point to the
  polygon changes by 2pi over one circuit.
}
\value{
  A logical vector taking value \code{TRUE} where the point lies inside the 
  polygon and \code{FALSE} otherwise.
}
\author{Peter Perkins (peter@caliban.ucsd.edu) on S-news: http://math.yorku.ca/Who/Faculty/Monette/S-news/0474.html}
\examples{
test <- as.matrix(expand.grid(x=0:10,y=0:10))
poly <- jitter(rbind(c(1,3), c(7,5), c(3,9)))[c(1:3,1),]
inside <- region.inside(test, poly)
plot(test, col=ifelse(inside,"red","black"), pch=16)
lines(poly)
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
