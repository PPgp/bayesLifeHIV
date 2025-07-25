
.onLoad <- function (lib, pkg) {
    library.dynam("bayesLifeHIV", pkg, lib)
    unlockBinding(".e0options", asNamespace("bayesLife"))
    if ("bayesLife" %in% loadedNamespaces()) {
        assign(".e0options", bayesLifeHIVenv$e0hivoptions, 
               envir = asNamespace("bayesLife"))
  }
}

.onUnload <- function (libpath) {
  library.dynam.unload("bayesLifeHIV", libpath)
}
