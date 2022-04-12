do.e0.proj.hiv <- function(x, beta, l.start, kap, ndart.traj,  
                           n.proj = 11, p1 = 9, p2 = 9, const.var = FALSE){
    proj <- c(l.start, rep(NA, n.proj))
    for(a in 2:(n.proj + 1)){
        proj[a] <- (proj[a-1] + g.dl6(x, proj[a-1], p1 = p1, p2 = p2) + 
                    beta*ndart.traj[a-1] +
                    rnorm(1, mean=0, sd=(kap*if(const.var) 1 else loess.lookup.hiv(proj[a-1], TRUE)))
                    )
    }
    return(proj[-1])
}

generate.e0hiv.trajectory <- function(..., traj, pred.env) {
    ccode <- as.character(pred.env$country.obj$code)
    if(pred.env$country.obj$index %in% pred.env$hiv.country.idx) 
        return(do.e0.proj.hiv(..., beta = pred.env$var.beta[traj,], 
                              ndart.traj = pred.env$ndart.trajs[ccode, traj,]))
    # non-HIV projections
    return(bayesLife:::generate.e0.trajectory(...))
}
    
    
e0hiv.predict <- function(mcmc.set = NULL, end.year = 2100, 
                       sim.dir = file.path(getwd(), 'bayesLifeHIV.output'),
                       replace.output = FALSE, predict.jmale = TRUE, 
                       nr.traj = NULL, thin = NULL, burnin = 10000, 
                       use.diagnostics = FALSE, hiv.countries = NULL, 
                       my.art.file = NULL, my.hivtraj.file = NULL, 
                       scale.hivtraj = FALSE, scale.hivtraj.tofile = NULL,
                       save.as.ascii = 1000, 
                       start.year = NULL, output.dir = NULL, 
                       low.memory = TRUE, ignore.last.observed = FALSE,
                       seed = NULL, verbose = TRUE, ...){
	if(!is.null(mcmc.set)) {
		if (class(mcmc.set) != 'bayesLife.mcmc.set') {
			stop('Wrong type of mcmc.set. Must be of type bayesLife.mcmc.set.')
		}
	} else {                
		mcmc.set <- get.e0.mcmc(sim.dir, low.memory=low.memory, verbose=verbose)
	}
	if(!is.null(seed)) set.seed(seed)
	# Get argument settings from existing convergence diagnostics
	if(use.diagnostics) {
	    diagpars <- bayesLife:::get.nr.traj.burnin.from.diagnostics(mcmc.set$meta$output.dir, verbose = verbose)
	    nr.traj <- diagpars$nr.traj
	    burnin <- diagpars$burnin
	}
    if(!is.null(hiv.countries)) # convert to index
        hiv.countries <- which(mcmc.set$meta$regions$country_code %in% hiv.countries)
	pred <- make.e0hiv.prediction(mcmc.set, end.year = end.year, replace.output = replace.output,  
					nr.traj = nr.traj, thin = thin, burnin = burnin, hiv.countries = hiv.countries, 
					save.as.ascii = save.as.ascii, start.year = start.year,
					output.dir = output.dir, my.art.file = my.art.file, 
					my.hivtraj.file = my.hivtraj.file, scale.hivtraj = scale.hivtraj,
					scale.hivtraj.tofile = scale.hivtraj.tofile, ignore.last.observed = ignore.last.observed,
					verbose = verbose)
	if(predict.jmale && mcmc.set$meta$sex == 'F')
		pred <- e0.jmale.predict(pred, ..., save.as.ascii=save.as.ascii, verbose=verbose)
	invisible(pred)
}

e0hiv.prediction.setup <- function(mcmc.set, ...) {
    convert.to.double <- function(dt) {
        fcols <- setdiff(colnames(dt)[sapply(dt, class) != "numeric"], "country_code")
        if(length(fcols) > 0)
            dt[,(fcols):= lapply(.SD, as.double), .SDcols = fcols]
        dt
    }
    setup <- bayesLife:::e0.prediction.setup(mcmc.set = mcmc.set, ...)
    hiv.country.codes <- c()
    
    ##load beta
    var.beta.names <- c('betanonART')
    var.beta <- get.e0.parameter.traces(setup$load.mcmc.set$mcmc.list, var.beta.names, burnin = 0)
    # load hiv trajectories and convert to a 3d array
    hiv.env <- new.env()
    art <- if(is.null(setup$my.art.file)) read.e0.data.file("ARTcoverage.txt") else fread(setup$my.art.file)
    if("include_code" %in% colnames(art))
        art <- art[include_code == 1, ]
    if(is.null(setup$my.hivtraj.file)) {
        data("HIVprevTrajectories", envir = hiv.env)
        hiv.traj <- data.table(hiv.env$HIVprevTrajectories)
    } else {
        hiv.traj <- fread(setup$my.hivtraj.file)
        if(setup$scale.hivtraj)  {
            scale.to <- NULL
            if(!is.null(setup$scale.hivtraj.tofile))
                scale.to <- fread(setup$scale.hivtraj.tofile)
            hiv.traj <- data.table(scale.hiv.trajectories(
                            data.frame(hiv.traj, check.names = FALSE),
                            scale.to = data.frame(scale.to, check.names = FALSE)))
        }
    }
    # delete some columns
    for(col in c("include_code", "name", "country_name")) {
        if(col %in% colnames(hiv.traj)) hiv.traj[[col]] <- NULL
        if(col %in% colnames(art)) art[[col]] <- NULL
    }
    # merge art and hiv and compute nonart
    artl <- melt(convert.to.double(copy(art)), id.vars = "country_code", variable.name = "period", value.name = "art")
    hivl <- melt(convert.to.double(copy(hiv.traj)), id.vars = c("country_code", "Trajectory"), variable.name = "period", value.name = "hiv")
    hiv.art <- merge(hivl, artl, by = c("country_code" , "period"))
    hiv.art[, nonart := hiv * (100 - art)/100]
    hiv.art[, year := as.integer(substr(period, 1,4))]
    if(setup$year.step > 1)
        hiv.art[, year := year + 3]
    # keep only projected time periods
    hiv.art <- hiv.art[year %in% setup$proj.middleyears,]
    # check if all years available
    uyears <- unique(hiv.art$year)
    missing <- setdiff(setup$proj.middleyears[-1], uyears)
    if(length(missing) > 0)
        stop("HIV trajectories or ART data missing for years ", 
             paste(missing, collapse = ", "))
    hiv.country.idx <- if(is.null(setup$hiv.countries)) which(setup$meta$regions$hiv.pred) else setup$hiv.countries
    hiv.country.codes <- setup$meta$regions$country_code[hiv.country.idx]
    if(any(! hiv.country.codes %in% hiv.art$country_code))
        stop("HIV trajectories or ART data missing for countries ", 
             paste(setdiff(hiv.country.codes, unique(hiv.art$country_code)), collapse = ", "))
    # keep only hiv countries
    hiv.art <- hiv.art[country_code %in% hiv.country.codes,]
    # convert back to wide format
    nonartl <- hiv.art[ , .(country_code, Trajectory, year, nonart)]
    nonart <- dcast(nonartl, country_code + Trajectory ~ year, value.var = "nonart") # wide format

    # get a 3d array of the deltas
    mid.years.char <- as.character(sort(uyears))
    nr.hiv.traj <- length(unique(nonart$Trajectory))
    mid.years.minus1 <-  mid.years.char[-length(mid.years.char)]
    nondart.trajs <- array(NA, dim=c(length(hiv.country.codes), nr.hiv.traj,
                                     length(mid.years.minus1)),
                           dimnames = list(hiv.country.codes, NULL, mid.years.minus1))
    for(cntry in hiv.country.codes) {
        tr <- as.matrix(nonart[country_code == cntry, mid.years.char, with = FALSE])
        nondart.trajs[as.character(cntry),, mid.years.minus1] <- tr[, 2:ncol(tr)] - tr[, 1:(ncol(tr)-1)]
    }
    setup$pred.env$ndart.trajs <- nondart.trajs
    setup$pred.env$hiv.traj.idx <- sample(1:nr.hiv.traj, setup$nr_simu, replace=TRUE)
    setup$pred.env$var.beta <- var.beta
    setup$hiv.country.codes <- hiv.country.codes
    return(setup)
}

make.e0hiv.prediction <- function(mcmc.set, pred.options = NULL, ...){
	# if 'countries' is given, it is an index. The same for hiv.countries.
    if(!is.null(pred.options))
        e0pred.options(pred.options)
    setup <- e0hiv.prediction.setup(mcmc.set, ...)
    pred <- bayesLife:::run.e0.projection.for.all.countries(setup, traj.fun = "generate.e0hiv.trajectory")
    pred$hiv.country.codes <- setup$hiv.country.codes
    bayesLife:::write.to.disk.prediction(pred, setup)
	invisible(pred)
}
