`predict.combinedGradientForest` <-
function (object, newdata, extrap=TRUE, ...)
{
    if (!inherits(object,"combinedGradientForest"))
      stop(paste("'object' must be a combinedGradientForest object"))
    linfun <- function(xold,yold,xnew)
        yold[1] + (xnew-xold[1])*diff(yold)/diff(xold)
    if (missing(newdata)) {
      n.gf <- length(object$nspec)
      n.preds <- sapply(object$gf.names, length)
      use.preds <- n.gf == n.preds
      use.preds.names <- names(n.preds[use.preds])
      if(!all(use.preds))
        warning("predict.combinedGradientForest called with missing newdata. By default, fitting data X is used as newdata, but not all gradientForest objects in the combinedGradientForest contain the same sets of predictors. Using the common subset of predictors from X.")
      newdata <- lapply(object$X, function(X) X[use.preds.names])
      newdata <- do.call("rbind",newdata)
      if(ncol(newdata) == 0)
        stop("predict.combinedGradientForest called with missing newdata, but no predictor appears in all gradientForest objects that make up the combinedGradientForest object. No meaningful prediction can be done, please provide newdata, eg. newdata = gf.combined$X[[1]].")
    }
    if(!inherits(newdata,"data.frame"))
        stop("newdata must be a data.frame")
    newnames <- names(newdata)
    if(!all(ok <- newnames %in% names(object$gf.names))) {
        badnames <- paste(newnames[!ok], collapse=", ")
        stop(paste("the following predictors are not in the gradientForest:\n\t",badnames,sep=""))
    }
    for (varX in newnames) {
        ci <- cumimp(object, varX, ...)
        xold <- range(ci$x)
        yold <- range(ci$y)
        xnew <- range(newdata[,varX],na.rm=T)
        if (extrap)
          ynew <- linfun(xold, yold, xnew)
        else 
          ynew <- yold
        if (xnew[1] < xold[1]) {
            ci$x <- c(xnew[1],ci$x)
            ci$y <- c(ynew[1],ci$y)
        }
        if (xnew[2] > xold[2]) {
            ci$x <- c(ci$x,xnew[2])
            ci$y <- c(ci$y,ynew[2])
        }
        f <- approxfun(ci, rule = 2)  
        newdata[,varX] <- f(newdata[,varX])     
    }
    class(newdata) <- c("predict.gradientForest", "data.frame")
    newdata
}
