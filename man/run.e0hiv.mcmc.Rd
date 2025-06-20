\name{run.e0hiv.mcmc}
\alias{run.e0hiv.mcmc}

\title{
Running Bayesian Hierarchical Model for Life Expectancy with HIV/AIDS
}
\description{
Runs Markov Chains Monte Carlo (MCMCs) for simulating the life expectancy for all countries of the world while taking into account HIV/AIDS epidemics, as described in Godwin and Raftery (2017).
}
\usage{
run.e0hiv.mcmc(sex = c("Female", "Male"), nr.chains = 3, iter = 160000, 
    output.dir = file.path(getwd(), "bayesLifeHIV.output"), 
    thin = 10, replace.output = FALSE, 
    start.year = 1873, present.year = 2020, wpp.year = 2019, 
    my.hiv.file = NULL, my.art.file = NULL, 
    mcmc.options = NULL, \dots)
}

\arguments{
  \item{sex}{Sex for which to run the simulation.}
  \item{nr.chains}{Number of MCMC chains to run.}
  \item{iter}{Number of iterations to run in each chain. In addition to a single value, it can have the value \sQuote{auto} for an automatic assessment of the convergence. In such a case, the function runs for the number of iterations given in the global option \code{auto.conf} list (see \code{\link[bayesLife]{e0mcmc.options}}), then checks if the MCMCs converged (using the \code{auto.conf} settings). If it did not converge, the procedure is repeated until convergence is reached or the number of repetition exceeded \code{auto.conf$max.loops}.}
  \item{output.dir}{Directory which the simulation output should be written into.}
  \item{thin}{Thinning interval between consecutive observations to be stored on disk.}
  \item{replace.output}{If \code{TRUE}, existing outputs in \code{output.dir} will be replaced by results of this simulation.}
  \item{start.year}{Start year for using historical data.}
  \item{present.year}{End year for using historical data.}
  \item{wpp.year}{Year for which WPP data is used. The functions loads a package called \pkg{wpp}\eqn{x} where \eqn{x} is the \code{wpp.year} and uses the \code{\link[wpp2019]{e0}*} datasets.}
  \item{my.hiv.file}{File name containing user-specified historical HIV prevalence for countries with past epidemics. It should have the same format as the default \code{\link{HIVprevalence}} dataset. Only countries are considered as epidemic in the estimation, if their \dQuote{include_code} column is 1.}
  \item{my.art.file}{File name containing user-specified historical values of ART coverage for countries with past epidemics. It has the same format as \code{my.hiv.file}. If not given, the default dataset, \code{\link{ARTcoverage}}, is used.}
  \item{mcmc.options}{List of options that overwrites global MCMC options as defined in \code{\link{e0mcmc.options}}. Type \code{using.bayesLifeHIV(); e0mcmc.options()} to view default values.}
  \item{\dots}{Additional parameters to be passed to \code{\link[bayesLife]{run.e0.mcmc}}.}
}
\details{
The function implements the Godwin and Raftery (2017) methodology which is an extension of the Chunn et al. (2013) model implemented in \code{\link[bayesLife]{run.e0.mcmc}} of the \pkg{\link[bayesLife]{bayesLife}} package. By default, the \code{\link[bayesLife]{run.e0.mcmc}} function does not include HIV/AIDS countries in the estimation. In fact, in \pkg{\link[bayesLife]{bayesLife}} they have to be treated as \dQuote{extra} countries (using the function \code{\link[bayesLife]{run.e0.mcmc.extra}}) in order to generate projections for the epidemic countries. The \code{run.e0hiv.mcmc} function uses ALL countries in the estimation, so no extra step is needed to generate projections of life expectancy for all countries of the world, and no need to run \pkg{\link[bayesLife]{bayesLife}} separately.

The function \code{run.e0hiv.mcmc} uses a set of global options (for priors, initial values etc.), possibly modified by the \code{mcmc.options} argument. One can also modify these options using \code{\link{e0mcmc.options}}. Call \code{using.bayesLifeHIV(); e0mcmc.options()} for the full set of options. 

See \code{\link[bayesLife]{run.e0.mmc}} for more details. 
}

\value{
	An object of class \code{\link[bayesLife]{bayesLife.mcmc.set}}.
}
\references{
Godwin, J. and Raftery, A.E. (2017): Bayesian projection of life expectancy accounting for the HIV/AIDS epidemic. Demographic Research 37(48):1549-1610.

J. L. Chunn, A. E. Raftery, P. Gerland, H. Sevcikova (2013): Bayesian Probabilis
tic Projections of Life Expectancy for All Countries. Demography 50(3):777-801. 
<doi:10.1007/s13524-012-0193-x>
}
\author{
Hana Sevcikova, Jessica Godwin
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[bayesLife]{run.e0.mmc}},
\code{\link[bayesLife]{get.e0.mcmc}}, \code{\link[bayesLife]{summary.bayesLife.mcmc.set}}, \code{\link[bayesLife]{e0mcmc.options}}, \code{\link{e0hiv.predict}}.
}
\examples{
\dontrun{
m <- run.e0hiv.mcmc(nr.chains = 1, iter = 5, thin = 1, verbose = TRUE)
summary(m)
hiv.countries.est(m$meta)
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{distribution}
\keyword{multivariate}
