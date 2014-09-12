Origin2Beta <- function(z,beta) (beta+z)/(1+Conj(beta)*z) # map disk to disk sending 0 to beta
Beta2Origin <- function(z,beta) (z-beta)/(1-Conj(beta)*z) # map disk to disk sending beta to 0
ROrigin2Beta <- function(x,x0) C2R(Origin2Beta(R2C(x),R2C(x0)))
RBeta2Origin <- function(x,x0) C2R(Origin2Beta(R2C(x),R2C(x0)))
R2matrix <- function(x) matrix(x,,2,T)
matrix2R <- function(x) as.vector(t(x))
matrix2C <- function(x) R2C(matrix2R(x))
C2matrix <- function(x) R2matrix(C2R(x))
matrixOrigin2Beta <- function(x,x0) R2matrix(C2R(Origin2Beta(R2C(matrix2R(x)),R2C(x0))))
matrixBeta2Origin <- function(x,x0) R2matrix(C2R(Beta2Origin(R2C(matrix2R(x)),R2C(x0))))
C2R <- function(z) as.vector(rbind(Re(z),Im(z)))
listxy2R <- function(lxy) as.vector(do.call("rbind",lxy))
listxy2C <- function(x) R2C(listxy2R(x))
R2C <- function(ww) if(is.complex(ww)) ww else complex(real=matrix(ww,,2,T)[,1],imag=matrix(ww,,2,T)[,2])
C2listxy <- function(z) list(x=Re(z),y=Im(z))
R2listxy <- function(z) C2listxy(R2C(z))


region.inside <- 
function(test, vertices)
{
#
# Originally from the Snews list - I have lost track of the original author
#
# determine whether test points are inside a simply closed polygon.
# this routine uses EUCLIDEAN geometry, not SPHERICAL.
#
# test:         matrix (or 2-vector) of test points
# vertices:     matrix of polygon vertices, last vertex same as first
#
# col1 = x-coords, col2 = y-coords
#
	n <- nrow(vertices) - 1
	if(!is.matrix(test)) test <- matrix(test, nrow = 1, byrow = T)	#
# translate to coords centered at each test point
	x <-  - outer(test[, 1], vertices[1:n, 1], "-")
	y <-  - outer(test[, 2], vertices[1:n, 2], "-")	#
#
# tan(theta) = (v1 x v2) / (v1 * v2)
	i2 <- c(2:n, 1)
	theta <- atan2(x * y[, i2] - y * x[, i2], x * x[, i2] + y * y[, i2])	#
#
# sum angles about each test point
	theta <- abs(theta %*% rep(1, n))
	ifelse(abs(theta - 2 * pi) < 1e-006, T, F)
}



conformalFit <- function(xypoly,wc=c(mean(xypoly$x),mean(xypoly$y)),nptsq=12) {
  n <- length(xypoly$x)
  dx <- diff(xypoly$x[c(n,1:n,1)])
  dy <- diff(xypoly$y[c(n,1:n,1)])
  Cos <- (dx[-1]*dx[-n-1]+dy[-1]*dy[-n-1])/sqrt(dx[-1]^2+dy[-1]^2)/sqrt(dx[-n-1]^2+dy[-n-1]^2)
  Sin <- (dy[-1]*dx[-n-1]-dx[-1]*dy[-n-1])/sqrt(dx[-1]^2+dy[-1]^2)/sqrt(dx[-n-1]^2+dy[-n-1]^2)
  betam <- -atan2(Sin,Cos)/pi
  if (all(betam>=0)) stop("the polygon must be entered in anti-clockwise order")
  if (any(betam>0)) warning("the polygon may have been input in clockwise order")
  w <- as.vector(do.call("rbind",xypoly))
  #z <- C2R(exp(complex(modulus=1,argument=seq(0,2*pi,len=n+1)[-1])))
  wc <- wc
  fit <- .C("fit_",n=as.integer(n),w=as.double(w),wc=as.double(wc),
            betam=as.double(betam),nptsq=as.integer(nptsq),
            tol=double(1),errest=double(1),c=double(2),z=double(2*n),
            qwork=double(nptsq*(2*n+3)))
  fit
}

conformalPredict <- function(fit,x,polygon2disk=TRUE) {
  npred <- length(x)/2
  if (polygon2disk) {
    pred <- .C("predictr_",n=as.integer(fit$n),c=as.double(fit$c),
               z=as.double(fit$z),wc=as.double(fit$wc),w=as.double(fit$w),
               betam=as.double(fit$betam),nptsq=as.integer(fit$nptsq),
               qwork=as.double(fit$qwork),ww=as.double(x),
               npred=as.integer(npred),zz=double(2*npred))
    pred$zz
  } else {
    pred <- .C("predict_",n=as.integer(fit$n),c=as.double(fit$c),
               z=as.double(fit$z),wc=as.double(fit$wc),w=as.double(fit$w),
               betam=as.double(fit$betam),nptsq=as.integer(fit$nptsq),
               qwork=as.double(fit$qwork),zz=as.double(x),
               npred=as.integer(npred),ww=double(2*npred))
    pred$ww
  }
}


