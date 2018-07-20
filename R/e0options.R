e0hiv.options.default <- function() {
    opts <- bayesLife:::.e0options
    opts$mcmc <- e0hiv.mcmc.options.default()
    opts
}

e0hiv.mcmc.options.default <- function() {
    opts <- bayesLife:::e0.mcmc.options.default()
    within(opts, {
        betanonART <- structure(list(ini = -0.5), npar = 1)
        world.parameters <- c(world.parameters, betanonART = 1)
        outliers <- c(-10, 10)
        estimation.function <- "e0hiv.mcmc.sampling"
        meta.ini.fun <- "e0hiv.meta.ini"
        dlcurves.function <- "e0hiv.get.dlcurves"
        include.hiv.countries <- TRUE
    })
}

using.bayesLifeHIV <- function() {
    # overwrite bayesLife options with the HIV ones
    opt.default <- e0hiv.options.default()
    for (item in names(opt.default))
        bayesLife:::e0.options(item, opt.default[[item]])
}

#.e0hivoptions <- e0hiv.options.default()
using.bayesLifeHIV()
