if(getRversion() >= "2.15.1") utils::globalVariables("loess.sd")


loess.lookup.hiv <- function(look, is.hiv) {
    find.look <- function(value, hiv) {
        look.in <- if(hiv) loess.sd$hiv else loess.sd
        idx <- cut(value, look.in$x, labels=FALSE, include.lowest = TRUE)
        look.in$y[idx]
    }
    mapply(find.look, look, is.hiv)
}

compute.hiv.residuals <- function(sim.dir, burnin = 1000) {
    # Compute residuals (non-hiv and hiv) from an existing simulation.
    # The simulation should have run with constant.variance = TRUE.
    mc <- get.e0.mcmc(sim.dir)
    meta <- mc$meta
    cs.par.names <- c('Triangle.c', 'k.c', 'z.c')
    nT <- dim(meta$e0.matrix)[1]-1
    nC <- get.nrest.countries(meta)
    residuals <- obsdata <- matrix(NA, nrow=nC, ncol=nT)
    beta <- get.e0.parameter.traces(mc$mcmc.list, "betanonART", burnin = burnin)
    dlfunc <- function(t) {
        mean(apply(cbind(cs.pars, beta), 1, 
            function(pars) g.dl6(pars[1:6], l = e0.obs[t], 
                                p1 = meta$dl.p1, p2 = meta$dl.p2) + pars[7]*dlt[t]))
    }
    # Compute residuals
    for(cntry in 1:nC) {
        country.obj <- get.country.object(cntry, meta, index = TRUE)
        cs.pars <- get.e0.parameter.traces.cs(mc$mcmc.list, country.obj, 
                                              cs.par.names, burnin = burnin)
        e0.obs <- meta$e0.matrix[, cntry]
        dlt <- meta$dlt.nart[, cntry]
        dl <- sapply(1:nT, dlfunc)
        residuals[cntry, ] <- abs(e0.obs[-length(e0.obs)] + dl - e0.obs[-1])
        obsdata[cntry, ] <- e0.obs[-1]
    }
    res.hiv <- residuals[meta$regions$hiv.est,]
    obs.hiv <- obsdata[meta$regions$hiv.est,]
    res.nohiv <- residuals[!meta$regions$hiv.est,]
    obs.nohiv <- obsdata[!meta$regions$hiv.est,]
    resdf <- data.frame(x = as.vector(obs.hiv), y = as.vector(res.hiv), is.hiv = TRUE)
    resdf <- na.omit(rbind(resdf, 
                           data.frame(x = as.vector(obs.nohiv), 
                                      y = as.vector(res.nohiv), 
                                      is.hiv = FALSE)))
    return(resdf)
}

compute.hiv.loess <- function(sim.dir = NULL, burnin = 1000, residuals = NULL, 
                          merge.hiv.tails = TRUE) {
    # Compute the loess curves (non-hiv and possibly hiv) from either 
    # an existing simulation (sim.dir) or from given residuals.
    # If sim.dir is used, the simulation should have run with constant.variance=TRUE.
    # If residuals are given, it should have columns x and y and optionally logical is.hiv.
    if(is.null(residuals))
        residuals <- compute.hiv.residuals(sim.dir, burnin = burnin)

    resdf.hiv <- loess.sd.hiv <- NULL
    if("is.hiv" %in% colnames(residuals)) {
        resdf.nohiv <- residuals[residuals$is.hiv == FALSE,]
        resdf.hiv <- residuals[residuals$is.hiv == TRUE,]
    } else 
        resdf.nohiv <- residuals
    loess.sd <- bayesLife:::do.compute.loess(resdf.nohiv)
    if(!is.null(resdf.hiv)) {
        loess.sd.hiv <- bayesLife:::do.compute.loess(resdf.hiv)
        if(merge.hiv.tails) {
            # Find intersections
            splfun <- splinefun(loess.sd$x, loess.sd$y)
            splhivfun <- splinefun(loess.sd.hiv$x, loess.sd.hiv$y)
            cross1 <- optimise(f=function(x) (splfun(x) - splhivfun(x))^2, c(30, 50))$minimum
            cross2 <- optimise(f=function(x) (splfun(x) - splhivfun(x))^2, c(70, 100))$minimum
            # Replace parts of hiv curves (to the left from cross1 and to the right
            # from cross2) with the non-hiv equivalents
            idx <- which(loess.sd.hiv$x < cross1 | loess.sd.hiv$x > cross2)
            loess.sd.hiv$x <- loess.sd.hiv$x[-idx]
            loess.sd.hiv$y <- loess.sd.hiv$y[-idx]
        
            idx <- which(loess.sd$x < cross1)
            loess.sd.hiv$x <- c(loess.sd$x[idx], loess.sd.hiv$x)
            loess.sd.hiv$y <- c(loess.sd$y[idx], loess.sd.hiv$y)
        
            idx <- which(loess.sd$x > cross2)
            loess.sd.hiv$x <- c(loess.sd.hiv$x, loess.sd$x[idx])
            loess.sd.hiv$y <- c(loess.sd.hiv$y, loess.sd$y[idx])
        }
        loess.sd$hiv <- loess.sd.hiv
    }
    return(loess.sd)
}
