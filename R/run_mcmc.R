if(getRversion() >= "2.15.1") utils::globalVariables(c("loess_sd", ".e0hivoptions"))
data(loess_sd, envir = environment(), package = "bayesLifeHIV")

e0hivmcmc.options <- function(...) {
    e0hiv.options("mcmc", ...)
}

e0hivpred.options <- function(...) {
    e0hiv.options("pred", ...)
}


e0hiv.options <- function(what, ...) {
    # this code was adapted from mclust.options 
    current <- .e0hivoptions
    if (nargs() == 1) return(current[[what]])
    args <- list(...)
    if (length(args) == 1 && is.null(names(args))) {
        arg <- args[[1]]
        switch(mode(arg), 
               list = args <- arg, 
               character = return(current[[what]][[arg]]), 
               stop("Invalid argument: ", dQuote(arg))
        )
    }
    if (length(args) == 0) return(current[[what]])
    if (is.null(names(args))) stop("Options must be given by name")
    current[[what]] <- modifyList(current[[what]], args)
    .e0hivoptions <<- current
    invisible(current[[what]])
}

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



run.e0hiv.mcmc <- function(sex = c("Female", "Male"), nr.chains = 3, iter = 160000, 
							output.dir = file.path(getwd(), 'bayesLifeHIV.output'), 
                            thin = 10, replace.output = FALSE,
                            start.year = 1873, present.year = 2015, wpp.year = 2017,
							mcmc.options = NULL, ...) {
    #defoptions <- e0hiv.mcmc.options.default()
    old.opts <- e0hivmcmc.options()
    hiv.opts <- old.opts
    if(!is.null(mcmc.options))
        hiv.opts <- e0hivmcmc.options(mcmc.options)
        #hiv.opts <- modifyList(hiv.opts, mcmc.options)

    res <- run.e0.mcmc(sex = sex, nr.chains = nr.chains, iter = iter, output.dir = output.dir,
                thin = thin, replace.output = replace.output, start.year = start.year,
                present.year = present.year, wpp.year = wpp.year, mcmc.options = hiv.opts
                )
    e0hivmcmc.options(old.opts)
    invisible(res)
}

read.e0.data.file <- function(file)
    read.delim(file = file.path(find.package("bayesLifeHIV"), 
                                "data", file), comment.char = "#", 
               check.names = FALSE)

e0hiv.meta.ini <- function(meta) {
	nT <- nrow(meta$e0.matrix)
	loessSD <- meta$loessSD
	dlt.nart <- NULL
	# read HIV and ART files
	hiv <- read.e0.data.file("HIVprevalence.txt")
	if("include_code" %in% colnames(hiv))
        hiv <- hiv[hiv$include_code == 1,]
	art <- read.e0.data.file("ARTcoverage.txt")
	if("include_code" %in% colnames(art))
	    art <- art[art$include_code == 1,]
	# align the two datasets
	epi.cntries <- intersect(hiv$country_code, art$country_code)
	hiv <- hiv[hiv$country_code %in% epi.cntries,]
	rownames(hiv) <- hiv$country_code
	hiv <- hiv[,-which(colnames(hiv) %in% c("country_code", "name", "country_name", "include_code"))]
	hiv <- hiv[as.character(epi.cntries),]
	art <- art[art$country_code %in% epi.cntries,]
	rownames(art) <- art$country_code
	art <- art[,-which(colnames(art) %in% c("country_code", "name", "country_name", "include_code"))]
	art <- art[as.character(epi.cntries),]
	# align columns
	comcols <- intersect(colnames(hiv), colnames(art))
	hiv <- hiv[,comcols]
	art <- art[,comcols]
	# create nonART dataset
	nonart <- hiv * (100 - art)/100
	# put middle years as colnames
	colnames(nonart) <- as.integer(substr(colnames(nonart), 1,4))+3
	# match years with e0 matrix
	nonart.all <- matrix(0, ncol = nT, nrow = nrow(nonart), 
	                     dimnames = list(rownames(nonart), rownames(meta$e0.matrix)))
	nonart <- nonart[, colnames(nonart) %in% rownames(meta$e0.matrix)]
	nonart.all[, colnames(nonart)] <- as.matrix(nonart)
	meta$regions$hiv.est <- meta$regions$country_code %in% epi.cntries
	# add countries not included in nonart, i.e. all non-hiv countries
	missing.countries <- meta$regions$country_code[!meta$regions$hiv.est]
	nonart <- rbind(nonart.all, 
	                matrix(0, nrow=length(missing.countries), ncol=ncol(nonart.all),
	                       dimnames=list(NULL, colnames(nonart.all))))
	# put into the same order as e0.matrix (the codes must match, i.e. no missing values)
	cntries.order <- match(meta$regions$country_code, c(epi.cntries, missing.countries))
	
	nonart <- nonart[cntries.order,]
	rownames(nonart) <- c(epi.cntries, missing.countries)[cntries.order]
	nonart <- t(nonart)
	dlt.nart <- meta$d.ct # to have a matrix of the right shape
	
	for(i in 2:nT) {
	    isna0 <- is.na(meta$e0.matrix[i-1,])
	    if(all(isna0)) next
	    dlt.nart[i-1, ] <- nonart[i, ] - nonart[i-1, ]
		dlt.nart[i-1, is.na(meta$d.ct[i-1,])] <- NA # outliers
		if (!meta$constant.variance)
			loessSD[i-1,!isna0] <- loess.lookup.hiv(meta$e0.matrix[i-1, !isna0], 
			                                             meta$regions$hiv.est[!isna0])
	}
	suppl <- meta$suppl.data
	loessSD.suppl <- suppl$loessSD
	if(!is.null(suppl$e0.matrix)) {
		# add first time point of the observed data to get the last increment of the supplemental data
		data.suppl <- rbind(suppl$e0.matrix, meta$e0.matrix[1, suppl$index.to.all.countries])
		nT <- nrow(data.suppl)
		for(i in 2:nT) {
			isna0 <- is.na(data.suppl[i-1,])
			if(all(isna0)) next
			loessSD.suppl[i-1,!isna0]<- if(meta$constant.variance) 1 else loess.lookup.hiv(data.suppl[i-1,!isna0], 
				                        meta$regions$hiv.est[suppl$index.to.all.countries][!isna0])
		}
		suppl$loessSD <- loessSD.suppl
	}
	meta$dlt.nart <- dlt.nart
	meta$loessSD <- loessSD
	meta$suppl.data$loessSD <- loessSD.suppl
	return(meta)
}

.e0hivoptions <- e0hiv.options.default()