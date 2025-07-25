library(bayesLifeHIV)

start.test <- function(name, wpp.year = NULL) cat('\n<=== Starting test of', name, if(!is.null(wpp.year)) paste0('(WPP ', wpp.year, ')') else '', '====\n')
test.ok <- function(name) cat('\n==== Test of', name, 'OK.===>\n')

test.options <- function(){
    test.name <- 'e0 HIV options'
    start.test(test.name)
    # HIV options should be available automatically
    mcopts5y <- e0mcmc.options()
    stopifnot("betanonART" %in% names(mcopts5y))
    stopifnot(mcopts5y$include.hiv.countries == TRUE)
    mcopts1y <- e0mcmc.options(annual = TRUE)
    stopifnot(all(mcopts1y$outliers == c(-5, 5)))

    # switch to bayesLife options
    bayesLife::using.bayesLife()
    stopifnot(! "betanonART" %in% names(e0mcmc.options()))
    
    # back to HIV options
    bayesLifeHIV::using.bayesLifeHIV()
    stopifnot("betanonART" %in% names(e0mcmc.options()))
    test.ok(test.name)
}

test.simulate5y <- function(wpp.year = 2019) {
	sim.dir <- tempfile()
    # run MCMC
    test.name <- 'estimating MCMC (5-year)'
	start.test(test.name, wpp.year)
    m <- run.e0hiv.mcmc(nr.chains = 1, iter = 20, thin = 1, output.dir = sim.dir, 
                     wpp.year = wpp.year, mcmc.options = list(buffer.size = 5))
    stopifnot(m$mcmc.list[[1]]$finished.iter == 20)
	stopifnot(get.total.iterations(m$mcmc.list, 0) == 20)
	stopifnot(m$meta$mcmc.options$buffer.size == 5)
	stopifnot(nrow(hiv.countries.est(m$meta)) == 58)
	stopifnot(nrow(hiv.countries.pred(m$meta)) == 21)
	ncountries <- nrow(get.countries.table(m))
	if(wpp.year == 2019)
	    stopifnot(ncountries == 201)
	if(wpp.year == 2024)
	    stopifnot(ncountries == 237)
	test.ok(test.name)
	
	# run prediction
	test.name <- 'running projections (5-year)'
	start.test(test.name, wpp.year)
	pred <- e0hiv.predict(m, burnin=0, verbose=FALSE)
	spred <- summary(pred)
	stopifnot(spred$nr.traj == 20)
	smpred <- summary(get.e0.jmale.prediction(pred))
	stopifnot(smpred$nr.traj == 20)
	test.ok(test.name)
	
	# obtaining DL curves
	e0 <- seq(40, 90, length=100)
	dl <- e0.country.dlcurves(e0, m, "Congo", burnin=10, predictive.distr = TRUE) 
	stopifnot(all(dim(dl)==c(10,100)))
	
	dlw <- e0.world.dlcurves(e0, m, burnin=10)
	stopifnot(all(dim(dlw)==c(10,100)))
	
	# check that the right function is called
	fnc <- bayesLife:::get.from.options("dlcurves.function", m$mcmc.list[[1]]$meta$mcmc.options)
	stopifnot(fnc == "e0hiv.get.dlcurves")
	
	# is betanonART included in traces?
	traces <- get.e0.parameter.traces(m$mcmc.list, burnin=0, thin=1) 
	stopifnot("betanonART" %in% colnames(traces))
	
	unlink(sim.dir, recursive=TRUE)
}


test.simulate1y <- function(wpp.year = 2022) {
    sim.dir <- tempfile()
    
    test.name <- 'running MCMC with annual data'
    start.test(test.name, wpp.year)
    m <- run.e0hiv.mcmc(iter = 5, nr.chains = 2, thin = 1, output.dir = sim.dir, 
                     annual = TRUE, present.year = 2020, wpp.year = wpp.year)
    stopifnot(get.total.iterations(m$mcmc.list, 0) == 10)
    stopifnot(all(1953:2020 %in% rownames(m$meta$e0.matrix)))
    test.ok(test.name)
    
    test.name <- 'running annual projections'
    start.test(test.name, wpp.year)
    pred <- e0hiv.predict(m, burnin=1, verbose = FALSE)
    spred <- summary(pred)
    tbl <- e0.trajectories.table(pred, "Japan")
    stopifnot(spred$nr.traj == 8)
    stopifnot(all(2021:2100 %in% spred$projection.years))

    mpred <- get.e0.jmale.prediction(pred)
    smpred <- summary(mpred)
    mtbl <- e0.trajectories.table(mpred, "Zambia")
    
    stopifnot(all(2021:2100 %in% smpred$projection.years))
    stopifnot('1956' %in% rownames(mtbl))
    
    test.ok(test.name)
    
    unlink(sim.dir, recursive=TRUE)
}

