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
                       use.diagnostics = FALSE, hiv.countries = NULL, save.as.ascii = 1000, 
                       start.year = NULL, output.dir = NULL, low.memory = TRUE, 
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
	pred <- make.e0hiv.prediction(mcmc.set, end.year=end.year, replace.output=replace.output,  
					nr.traj=nr.traj, thin=thin, burnin=burnin, hiv.countries = hiv.countries, 
					save.as.ascii=save.as.ascii, start.year=start.year,
					output.dir=output.dir, verbose=verbose)
	if(predict.jmale && mcmc.set$meta$sex == 'F')
		pred <- e0.jmale.predict(pred, ..., save.as.ascii=save.as.ascii, verbose=verbose)
	invisible(pred)
}

e0hiv.prediction.setup <- function(mcmc.set, ...) {
    setup <- bayesLife:::e0.prediction.setup(mcmc.set = mcmc.set, ...)
    hiv.country.codes <- c()
    
    ##load beta
    var.beta.names <- c('betanonART')
    var.beta <- get.e0.parameter.traces(setup$load.mcmc.set$mcmc.list, var.beta.names, burnin = 0)
    # load hiv trajectories and convert to a 3d array
    hiv.env <- new.env()
    art <- read.e0.data.file("ARTcoverage.txt")
    if("include_code" %in% colnames(art))
        art <- art[art$include_code == 1,]
    rownames(art) <- art$country_code
    art <- art[,-which(colnames(art) %in% c("country_code", "name", "country_name", "include_code"))]
    colnames(art) <- as.integer(substr(colnames(art), 1,4))+3
    data("HIVprevTrajectories", envir = hiv.env)
    trajs <- hiv.env$HIVprevTrajectories
    elim.cols <- which(colnames(trajs) %in% c("country_code", "Trajectory"))
    years <- colnames(trajs)[-elim.cols]
    mid.years <- substr(years, 1, 4)
    mid.years <- as.integer(mid.years) + 3
    mid.years.remove <- which(! mid.years %in% setup$proj.middleyears)
    if(length(mid.years.remove) > 0) {
        mid.years <- mid.years[-mid.years.remove]
        elim.cols <- c(elim.cols, which(colnames(trajs) %in% years[mid.years.remove]))
    }
    if(length(mid.years) < length(setup$proj.middleyears))
        stop("HIV trajectories missing for years ", 
             paste(setdiff(setup$proj.middleyears[-1], mid.years), collapse = ", "))
    colnames(trajs)[-elim.cols] <- mid.years
    mid.years.char <- as.character(sort(mid.years))
    hiv.country.idx <- if(is.null(get0("hiv.countries"))) which(setup$meta$regions$hiv.pred) else hiv.countries
    hiv.country.codes <- setup$meta$regions$country_code[hiv.country.idx]
    if(any(! hiv.country.codes %in% trajs$country_code))
        stop("HIV trajectories missing for countries ", 
             paste(setdiff(unique(trajs$country_code), hiv.country.codes), collapse = ", "))
    nr.hiv.traj <- length(unique(trajs$Trajectory))
    mid.years.minus1 <-  mid.years.char[-length(mid.years.char)]
    nondart.trajs <- array(NA, dim=c(length(hiv.country.codes), nr.hiv.traj,
                                     length(mid.years.minus1)),
                           dimnames = list(hiv.country.codes, NULL, mid.years.minus1))
    for(cntry in hiv.country.codes) {
        tr <- (as.matrix(trajs[trajs$country_code == cntry, mid.years.char]) * 
                   as.matrix(100 - art[rep(as.character(cntry), nr.hiv.traj), 
                                       mid.years.char])/100)
        nondart.trajs[as.character(cntry),, mid.years.minus1] <- tr[, 2:ncol(tr)] - tr[, 1:(ncol(tr)-1)]
    }
    setup$pred.env$ndart.trajs <- nondart.trajs
    setup$pred.env$hiv.traj.idx <- sample(1:nr.hiv.traj, setup$nr_simu, replace=TRUE)
    setup$pred.env$var.beta <- var.beta
    return(setup)
}

make.e0hiv.prediction <- function(mcmc.set, pred.options = NULL, ...){
	# if 'countries' is given, it is an index. The same for hiv.countries.
    hiv.pred.opts <- e0hivpred.options()
    hiv.mcmc.opts <- e0hivmcmc.options()
    if(!is.null(pred.options))
        hiv.pred.opts <- e0hivpred.options(pred.options)
    # overwrite bayesLife .e0options
    bayesLife:::e0mcmc.options(hiv.mcmc.opts)
    bayesLife:::e0pred.options(hiv.pred.opts)
    
    setup <- e0hiv.prediction.setup(mcmc.set, ...)
    pred <- bayesLife:::run.e0.projection.for.all.countries(setup, traj.fun = "generate.e0hiv.trajectory")
    pred$hiv.country.codes <- setup$hiv.country.codes
    bayesLife:::write.to.disk.prediction(pred, setup)
	invisible(pred)
}
