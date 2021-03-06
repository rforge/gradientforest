predictor.density.plot <-
function (obj, ...)
{
    bind.varXY <- function(obj) {do.call("rbind", lapply(names(obj), function(predictor) data.frame(predictor=predictor,value=obj[[predictor]]$x,density=obj[[predictor]]$y))) }
    dens <- do.call("rbind", lapply(names(obj$dens), function(gear) cbind(gf.name=gear,bind.varXY(obj$dens[[gear]]))))
    o <- order(-importance(obj))
    dens$predictor <- ordered(dens$predictor, levels=names(sort(-importance(obj))))
    dens$gf.name <- factor(dens$gf.name)
    spl <- trellis.par.get("superpose.line")
    n <- length(levels(dens$gf.name))
    spl$lwd[1] <- 2
    trellis.par.set("superpose.line",spl)
    print(xyplot(
      density~value|predictor,
      data=dens,
      groups=gf.name,
      type='l',
      scales=list(relation="free",cex=0.5),
      par.strip.text=list(cex=0.5),
      ylab="Density",
      xlab="Predictor value",
      as.table=T,
      key=list(
        space="top",
        columns=3,
        text=list(levels(dens$gf.name)),
        lines=list(type='l',col=spl$col[1:n],
        lwd=spl$lwd[1:n])
      )
    ))
 }
