\name{compute.hiv.loess}
\alias{compute.hiv.loess}
\alias{compute.hiv.residuals}

\title{
Computing Residuals and Fitting Local Polynomial Regression
}
\description{
Functions for computing residuals (separately for HIV and non-HIV data) from the observed life expectancy and MCMC estimation, and fitting a local polynomial regression. 
}
\usage{
compute.hiv.residuals(sim.dir, burnin = 1000)

compute.hiv.loess(sim.dir = NULL, burnin = 1000, residuals = NULL,
    merge.hiv.tails = TRUE)
}

\arguments{
  \item{sim.dir}{Directory with the MCMC estimation. In \code{compute.hiv.loess}, it is only used if \code{residuals} is not given.}
  \item{burnin}{Number of (unthinned) iterations to be discarded. In \code{compute.hiv.loess}, it is only used if \code{residuals} is not given.}
  \item{residuals}{Residuals can be computed outside of the \code{compute.hiv.loess} function and passed as the \code{residuals} argument. If not given, the \code{compute.hiv.loess} function calls \code{compute.hiv.residuals} internally, in which case \code{sim.dir} must be given.}
  \item{merge.hiv.tails}{Logical determining if an intersection should be found between the non-hiv loess and the the hiv loess.}
}

\details{
The Bayesian hierarchical model for life expectancy uses a lowess curve as a multiplier of the variance. The dataset is stored in the package as the \code{loess_sd} dataset. It has an HIV and non-HIV part. These functions can be used to re-compute this \code{loess_sd} dataset. In such a case, the simulation should be run with the argument \code{constant.variance = TRUE} (in \code{\link{run.e0hiv.mcmc}}).

The residuals are computed for each country as the absolute differences between the observed life expectancy increases and the mean of the estimated double logistic function at the corresponding life expectancy level.
}
\value{
\code{compute.hiv.residuals} returns a data frame with columns \sQuote{x} (life expectancy levels), \sQuote{y} (absolute residuals) and \sQuote{is.hiv} (\code{TRUE} for countries considered in the estimation as HIV/AIDS epidemic).

\code{compute.hiv.loess} also returns a data frame with columns \sQuote{x},  \sQuote{y} and \sQuote{hiv}, where \sQuote{x} is the same as before (with added a minimum and maximum), \sQuote{y} is the local polynomial fit with constant tails for non-HIV countries, and \sQuote{hiv} is the local polynomial fit for HIV countries. If there are no such countries, this column is missing.
}

\author{
Hana Sevcikova
}

\seealso{
\code{\link{run.e0hiv.mcmc}}, \code{\link{hiv.countries.est}}
}

\examples{
\dontrun{
sim.dir <- tempfile()
m <- run.e0hiv.mcmc(nr.chains = 1, iter = 100, thin = 1, 
  output.dir = sim.dir, verbose = TRUE, constant.variance = TRUE)
resid <- compute.hiv.residuals(sim.dir, burnin = 30)
lws <- compute.hiv.loess(residuals = resid)

# plot residuals and loess
plot(resid$x, resid$y, ylim = c(0, 4))
lines(lws$x, lws$y, col = "red")

unlink(sim.dir, recursive=TRUE)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ loess }% use one of  RShowDoc("KEYWORDS")

