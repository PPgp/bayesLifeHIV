scale.hiv.trajectories <- function(trajectories = NULL, scale.to = NULL,
                                   logit.adjust = 0.001) {
    # Scale given trajectories to a data frame given by scale.to.
    # scale.to should have a column country_code. Other columns should 
    # be named by time periods, e.g. 2025-2030, 2030-2035 etc.
    # Argument trajectories is a data frame with columns country_code, Trajectory and 
    # time periods. The country codes and time periods are matched 
    # between the two datasets. There must be the same number of trajectories 
    # for all countries.
    # If trajectories is not given, dataset HIVprevTrajectories is used. 
    # Default for scale.to is dataset HIVprevalence.
    # The scaling is done using an adjusted logit.
    
    rep.row <- function(x,n) matrix(rep(x,each=n), nrow=n)

    expit <- function(x, adjust) { # inverse logit with adjustment
        a <- 1 - 2 * adjust
        (exp(x)- 1)/(2*a*(1+exp(x))) + 1/2
    }
    env <- new.env()
    if(is.null(scale.to)) {
        data("HIVprevalence", envir = env)
        scale.to <- env$HIVprevalence
    }
    if(is.null(trajectories)) {
        data("HIVprevTrajectories", envir = env)
        trajectories <- env$HIVprevTrajectories
    }
    if(! "country_code" %in% colnames(trajectories))
        stop("Column country_code is missing in dataset trajectories.")
    if(! "country_code" %in% colnames(scale.to))
        stop("Column country_code is missing in dataset scale.to.")
    
    cntries <- sort(intersect(scale.to$country_code, trajectories$country_code))
    cntries.char <- as.character(cntries)
    cols <- intersect(colnames(scale.to), colnames(trajectories))
    dont.include <- c("include_code", "name", "Trajectory")
    if(any(cols %in% dont.include))
        cols <- cols[-which(cols %in% dont.include)]
    if(length(cols) <= 1)
        stop("Datasets trajectories and scale.to do not have time periods in common. Check column names.")
    
    unhiv <- scale.to[scale.to$country_code %in% cntries, cols]
    unhivmtx <- unhiv[, -which(colnames(unhiv) == "country_code")]
    rownames(unhivmtx) <- unhiv$country_code
    # reorder so that countries are sorted
    unhivmtx <- unhivmtx[cntries.char, ]
    
    if(! "Trajectory" %in% colnames(trajectories))
        stop("Column Trajectory is missing in dataset trajectories.")
    
    trajs <- trajectories[order(trajectories$country_code, trajectories$Trajectory),]
    trajs <- trajs[trajs$country_code %in% cntries, ]
    trajs.red <- trajs[, cols]
    trajs.red <- trajs.red[,-which(colnames(trajs.red) == "country_code")]
    spl.trajs <- split(trajs.red, trajs$country_code)
    # compute medians for each country
    trajs.med.spl <- lapply(spl.trajs, function(m) apply(m, 2, median))
    trajs.med <- do.call(rbind, trajs.med.spl)
    # reorder so that countries are sorted
    trajs.med.mtx <- trajs.med[cntries.char,]
    
    # Put the median dataset an the scale.to dataset into the same shape 
    # as trajectories by repeating rows nr.traj times.
    nr.trajs <- length(unique(trajectories$Trajectory))
    trajs.med.big <- apply(trajs.med.mtx, 2, rep.row, nr.trajs)
    trajs.mtx <- as.matrix(trajs.red)
    unhiv.big <- apply(unhivmtx, 2, rep.row, nr.trajs)
    
    x <- (logit(trajs.mtx/100, adjust=logit.adjust) + 
              logit(unhiv.big/100, adjust=logit.adjust) - 
              logit(trajs.med.big/100, adjust=logit.adjust))
    scaled <- 100*pmax(expit(x, adjust=logit.adjust), 0)
    scaled <- cbind(data.frame(country_code = trajs$country_code,
                    Trajectory = trajs$Trajectory), 
                    scaled)
    return(scaled)
}
