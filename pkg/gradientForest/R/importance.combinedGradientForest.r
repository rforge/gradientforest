`importance.combinedGradientForest` <-
function (x, type=c("Weighted","Raw","Species")[1], sort=TRUE, ...)
{
    if (!inherits(x,"combinedGradientForest"))
      stop(paste("'x' must be a combinedGradientForest object"))

    imp <- x$imp.rsq[,names(x$imp.rsq) != "predictor"]
    weighted <- rowSums(imp, na.rm=TRUE)/ncol(imp)
    if (sort)
      o <- order(-weighted)
    else o <- 1:length(weighted)
    nam <- x$imp.rsq$predictor
    res <- switch(pmatch(type,c("Weighted","Raw","Species")),
      weighted[o],
      rowSums(sweep(imp,2,x$rsq,"/"), na.rm=TRUE)[o]/ncol(imp),
      if (sort) sort(x$rsq,decreasing=T) else x$rsq
    )
    if (is.null(res))
      stop(paste('Unmatched type "',type,'". Expecting one of "Weighted", "Raw" or "Species"',sep=""))
    else res
}
