if(getRversion() >= "2.15.1") utils::globalVariables("loess_sd")
data(loess_sd, envir=environment(), package = "bayesLifeHIV")

e0hiv.get.dlcurves <- function(x, mcmc.list, country.code, burnin, 
                               nr.curves = 2000, predictive.distr = FALSE) {
    add.errors <- isTRUE(mcmc.list[[1]]$meta$constant.variance) && predictive.distr
    dlc <- bayesLife:::e0.get.dlcurves(x, mcmc.list, country.code, burnin, 
                                       nr.curves = nr.curves, predictive.distr = add.errors)
    
    if(!predictive.distr || add.errors) return(dlc)
    # obtain errors using the hiv loess curves
    is.hiv <- FALSE
	allerrs <- c()
    nr.curves.from.mc <- if (!is.null(nr.curves)) ceiling(max(nr.curves, 2000)/length(mcmc.list))
    						else NULL
    if(!is.null(country.code)) {
        country.obj <- get.country.object(country.code, mcmc.list[[1]]$meta)
        is.hiv <- mcmc.list[[1]]$meta$regions$hiv.est[country.obj$index]
    }
    loessSD <- loess.lookup.hiv(x, is.hiv = is.hiv)
    for (mcmc in mcmc.list) {
    	th.burnin <- bayesTFR:::get.thinned.burnin(mcmc,burnin)
    	thincurves.mc <- bayesTFR:::get.thinning.index(nr.curves.from.mc, 
            all.points = mcmc$length - th.burnin)
		omegas <- bayesLife:::load.e0.parameter.traces(mcmc, par.names = 'omega', 
		                                    burnin = th.burnin, 
		                                    thinning.index = thincurves.mc$index)

		errors <- matrix(NA, nrow = nrow(omegas), ncol = length(x))
		n <- ncol(errors)
		for(i in 1:nrow(errors))
			errors[i,] <- rnorm(n, mean = 0, sd = omegas[i]*loessSD)
        allerrs <- rbind(allerrs, errors)
    }
	return (dlc + allerrs)
}
