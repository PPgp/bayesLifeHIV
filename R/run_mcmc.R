if(getRversion() >= "2.15.1") utils::globalVariables(c("loess_sd"))
data(loess_sd, envir = environment(), package = "bayesLifeHIV")

run.e0hiv.mcmc <- function(sex = c("Female", "Male"), nr.chains = 3, iter = 160000, 
							output.dir = file.path(getwd(), 'bayesLifeHIV.output'), 
                            thin = 10, replace.output = FALSE,
                            start.year = 1873, present.year = 2015, wpp.year = 2017,
							my.hiv.file = NULL, my.art.file = NULL,
							mcmc.options = NULL, ...) {
    old.opts <- e0mcmc.options()
    if(is.null(mcmc.options)) mcmc.options <- list()
    mcmc.options$my.hiv.file <- my.hiv.file
    mcmc.options$my.art.file <- my.art.file
    
    res <- run.e0.mcmc(sex = sex, nr.chains = nr.chains, iter = iter, output.dir = output.dir,
                thin = thin, replace.output = replace.output, start.year = start.year,
                present.year = present.year, wpp.year = wpp.year, mcmc.options = mcmc.options, ...
                )
    
    e0mcmc.options(old.opts)
    invisible(res)
}

read.e0.data.file <- function(file)
    fread(file.path(find.package("bayesLifeHIV"), "data", file))

e0hiv.meta.ini <- function(meta) {
    convert.to.double <- function(dt) {
        fcols <- setdiff(colnames(dt)[sapply(dt, class) != "numeric"], "country_code")
        if(length(fcols) > 0)
            dt[,(fcols):= lapply(.SD, as.double), .SDcols = fcols]
        dt
    }
    
	nT <- nrow(meta$e0.matrix)
	loessSD <- meta$loessSD
	dlt.nart <- NULL
	opts <- meta$mcmc.options
	
	# read HIV and ART files and cleanup
	hiv <- if(is.null(opts$my.hiv.file)) read.e0.data.file("HIVprevalence.txt") else fread(opts$my.hiv.file)
	art <- if(is.null(opts$my.art.file)) read.e0.data.file("ARTcoverage.txt") else fread(opts$my.art.file)
	if("include_code" %in% colnames(hiv))
        hiv <- hiv[include_code == 1,]
	if("include_code" %in% colnames(art))
	    art <- art[include_code == 1,]
	# delete some columns
	for(col in c("include_code", "name")) {
	    if(col %in% colnames(hiv)) hiv[[col]] <- NULL
	    if(col %in% colnames(art)) art[[col]] <- NULL
	}
	
	# align the two datasets by converting to a long format and compute nonart
	hivl <- melt(convert.to.double(copy(hiv)), id.vars = "country_code", variable.name = "period", value.name = "hiv")
	artl <- melt(convert.to.double(copy(art)), id.vars = "country_code", variable.name = "period", value.name = "art")
    hiv.art <- merge(hivl, artl, by = c("country_code" , "period"))
    hiv.art[, nonart := hiv * (100 - art)/100]
    hiv.art[, year := as.integer(substr(period, 1,4))+3]

    # match years with e0 matrix and convert back to a wide format
    hiv.art <- hiv.art[year %in% rownames(meta$e0.matrix), ]
    nonartl <- hiv.art[ , .(country_code, year, nonart)]
    nonart <- dcast(nonartl, country_code ~ year, value.var = "nonart") # wide format
    rownames(nonart) <- nonart$country_code
    nonart[, country_code := NULL]
    
    # set countries with epidemics for estimation
	epi.cntries <- unique(hiv.art$country_code)
	meta$regions$hiv.est <- meta$regionsDT$hiv.est <- meta$regions$country_code %in% epi.cntries

	# add non-hiv countries
	missing.countries <- meta$regionsDT[hiv.est == FALSE, country_code]
	nonart.all <- matrix(0, ncol = nT, nrow = nrow(nonart), 
	                     dimnames = list(rownames(nonart), rownames(meta$e0.matrix)))
	nonart <- nonart[, colnames(nonart) %in% rownames(meta$e0.matrix), with = FALSE]
	nonart.all[, colnames(nonart)] <- as.matrix(nonart)
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

init.nodes.e0hiv <- function() {
    library(bayesLifeHIV)
}