`print.combinedGradientForest` <-
function(x,...)
{
  cat("\nCall:\n")
  cat(paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = " ")
  cat("\ngradientForest objects:\n")
  cat(paste(names(x$nspec), sep = " "), "\n", sep = " ")
  cat("\nNumber of Species:\n")
  print(x$nspec)
  cat("\nPredictors:\n")
  cat(paste(names(x$gf.names), sep = " "), "\n", sep = ", ", fill=T)
  invisible(x)
}
