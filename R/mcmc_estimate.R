e0hiv.mcmc.sampling <- function(mcmc, thin = 1, start.iter = 2, verbose = FALSE, verbose.iter = 10) {
	if (!is.null(mcmc$rng.state)) .Random.seed <- mcmc$rng.state
	niter <- mcmc$iter
	if (start.iter > niter) return(mcmc)
	mcmc$thin <- thin
	
	# Create an environment for the mcmc stuff in order to avoid 
	# copying of the whole mcmc list
	mcenv <- as.environment(mcmc)
	meta <- as.environment(mcenv$meta)
	
	ctrlenv <- bayesLife:::create.ctrl.env(mcenv, meta)
	ctrlenv$DLdata <- get.DLdata.for.hiv.estimation(meta, 1:ctrlenv$C)
	
	ctrlenv <- within(ctrlenv, {
	    X <- as.vector(t(meta$dlt.nart))
	    idxX <- !is.na(X) # can be NA where outliers
	    X <- X[idxX]
	    loess.vector <- as.vector(t(meta$loessSD))[idxX]
	    dct.vector <- as.vector(t(meta$d.ct))[idxX]
	    sdbeta <- 0.5*sd(dct.vector)/sd(X)
	    sigbetainv <- 1/sdbeta[1]^2
	    newDL <- meta$d.ct
	    newDL[] <- NA
	})

	for(iter in start.iter:niter) {
		if(verbose.iter > 0 && (iter %% verbose.iter == 0))
			cat('\nIteration:', iter, '--', date())
		bayesLife:::unblock.gtk('bDem.e0mcmc')
        
		update <- bayesLife:::update.mcmc.parameters(mcenv, ctrlenv, meta$mcmc.options)
		mcenv <- update$mc
		ctrlenv <- update$ctrl
		
		ctrlenv <- within(ctrlenv, {
		    newDL[] <- NA
		    for(country in 1:C) 
                newDL[colnames(DLdata[[country]])[DLdata[[country]]['post1950',]==1],country] <- dlf[[country]][DLdata[[country]]['post1950',]==1]
		    # Update beta - Gibbs Sampler
		    ##########################################
    	    XX <- tcrossprod(t(X)%*%diag(1/sqrt(mcenv$omega^2*loess.vector)))

            newDL.vector <- as.vector(t(newDL))[idxX]
		    sxy <- t(X)%*%(diag(1/(mcenv$omega^2*loess.vector))%*%(dct.vector-newDL.vector))

		    mcenv$betanonART <- rnorm(1,mean=sxy/(XX+sigbetainv),sd=sqrt(1/(XX+sigbetainv)))
        
            for(country in 1:C){
			    new <- DLdata[[country]]['observed.dct',DLdata[[country]]['post1950',]==1]-mcenv$betanonART*meta$dlt.nart[colnames(DLdata[[country]])[DLdata[[country]]['post1950',]==1],country]
			    DLdata[[country]]['dct',DLdata[[country]]['post1950',]==1] <- as.numeric(new)
			}	
		})
		# write samples simu/thin to disk
		mcenv <- bayesLife:::store.sample.to.disk(iter, niter, mcenv, verbose = verbose)

	}
	resmc <- as.list(mcenv)
	class(resmc) <- class(mcmc)
	return(resmc)
}


get.DLdata.for.hiv.estimation <- function(meta, countries) {
	DLdata <- list()
    T.suppl.end <- if(!is.null(meta$suppl.data$e0.matrix)) nrow(meta$suppl.data$e0.matrix) else 0
    for(country in countries) {
    	idx <- which(!is.na(meta$d.ct[, country]))
    	DLdata[[country]] <- matrix(NA, nrow=5, ncol=length(idx), 
    							dimnames=list(c('e0', 'dct', 'loess', 'observed.dct', 'post1950'), 
    										rownames(meta$d.ct[idx,])))
    	DLdata[[country]][1,] <- meta$e0.matrix[idx,country]
    	DLdata[[country]][2,] <- DLdata[[country]][4,] <- meta$d.ct[idx,country]
    	DLdata[[country]][3,] <- meta$loessSD[idx,country]
    	DLdata[[country]][5,] <- 1

    }
    if(T.suppl.end > 0) {
    	for(country in 1:ncol(meta$suppl.data$e0.matrix)) {
    		cidx <- meta$suppl.data$index.to.all.countries[country]
    		if (!is.element(cidx, countries)) next
    		idx <- which(!is.na(meta$suppl.data$d.ct[, country]))
    		if(length(idx) <= 0) next
    		start.col <- ncol(DLdata[[cidx]]) + 1
    		DLdata[[cidx]] <- cbind(DLdata[[cidx]], 
    								matrix(NA, nrow=5, ncol=length(idx),
    										dimnames=list(rownames(DLdata[[cidx]]), 
    											rownames(meta$suppl.data$d.ct[idx,]))))
    		end.col <- ncol(DLdata[[cidx]])
    		DLdata[[cidx]][1,start.col:end.col] <- meta$suppl.data$e0.matrix[idx,country]
       		DLdata[[cidx]][2,start.col:end.col] <- DLdata[[cidx]][4,start.col:end.col] <- meta$suppl.data$d.ct[idx,country]
          	DLdata[[cidx]][3,start.col:end.col] <- meta$suppl.data$loessSD[idx,country]
          	DLdata[[cidx]][5,start.col:end.col] <- 0
  		}
    }
	return(DLdata=DLdata)	
}

