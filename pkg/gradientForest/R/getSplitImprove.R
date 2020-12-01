`getSplitImprove` <-function(fit, X) {
#   return a data-frame: var name, rsq, var number, split value, improvement
    trees <- lapply(1:fit$ntree, function(k) try(getTree(fit, k),silent=TRUE)) #Nick Ellis 10/12/2009
    ok <- sapply(trees, class) != "try-error"
    if(!is.vector(ok)){
      ok <- apply(ok, 2, all)
    }
    tmp <- do.call("rbind", lapply((1:fit$ntree)[ok], function(k) cbind(tree =
                                                                          k, trees[[k]])))
    tmp <- tmp[tmp[,"status"]==-3 & zapsmall(tmp[,"improve"]) > 0,c("split var","split point","improve")]
    colnames(tmp) <- c("var_n","split","improve")
    rownames(tmp)<-NULL     #S.J. Smith 11/05/2009
    res <- cbind(data.frame(var=factor(names(X)[tmp[,"var_n"]], levels = names(X)),rsq=rep(fit$rsq[fit$ntree],nrow(tmp))),tmp)
    ok <- zapsmall(res[,"improve"]) > 0
    res[ok,] 
  }
