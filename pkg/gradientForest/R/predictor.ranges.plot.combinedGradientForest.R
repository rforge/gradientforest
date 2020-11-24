predictor.ranges.plot.combinedGradientForest <-
function (obj, ...)
{
    o <- order(-importance(obj))

    X_r <- do.call("rbind", lapply(names(obj$X), function(gf, obj){
      allpreds <- names(obj$gf.names)
      preds <- allpreds[sapply(obj$gf.names, function(pred) {gf %in% pred})]
      out <- data.frame(gf.name = gf, stack(obj$X[[gf]][preds]))
      names(out) <- c("gf.name","value","predictor")
      return(out)
    }, obj = obj)
    )
    browser()
    X_r$predictor <- ordered(X_r$predictor, levels=names(sort(-importance(obj))))
    print(bwplot(gf.name ~ value|predictor, X_r, scales=list(x=list(relation="free",cex=0.5)),
      par.strip.text=list(cex=0.6),xlab="Predictor value", as.table=T, ...))
}
