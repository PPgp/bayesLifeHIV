.onLoad <- function (lib, pkg) {
    library.dynam("bayesLifeHIV", pkg, lib)
    using.bayesLifeHIV()
}

.onUnload <- function (libpath) {
  library.dynam.unload("bayesLifeHIV", libpath)
}

.onAttach <- function(lib, pkg)
{
    # unlock .e0options variable allowing its modification
    #unlockBinding(".e0hivoptions", asNamespace("bayesLifeHIV"))
    unlockBinding(".e0options", asNamespace("bayesLife"))
    invisible()
}