rfNews <- function() {
    newsfile <- file.path(system.file(package="extendedForest"), "NEWS")
    file.show(newsfile)
}
