\name{e0hiv.predict}
\alias{e0hiv.predict}

\title{
Generating Posterior Trajectories of the Life Expectancy
}
\description{
Using the posterior parameter samples simulated by \code{\link{run.e0hiv.mcmc}} the function  generates posterior trajectories for the life expectancy for all countries of the world, 
while taking into account future HIV/AIDS epidemics.
}
\usage{
e0hiv.predict(mcmc.set = NULL, end.year = 2100, 
    sim.dir = file.path(getwd(), "bayesLife.output"), replace.output = FALSE, 
    predict.jmale = TRUE, nr.traj = NULL, thin = NULL, 
    burnin = 10000, use.diagnostics = FALSE, 
    hiv.countries = NULL, my.art.file = NULL, my.hivtraj.file = NULL,
    scale.hivtraj = FALSE, scale.hivtraj.tofile = NULL, 
    save.as.ascii = 1000, start.year = NULL, output.dir = NULL, 
    low.memory = TRUE, seed = NULL, verbose = TRUE, \dots)
}

\arguments{
  \item{mcmc.set}{Object of class \code{\link{bayesLife.mcmc.set}}. If it is \code{NULL}, the object is loaded from the directory given by \code{sim.dir}.
}
  \item{end.year}{End year of the prediction.}
  \item{sim.dir}{Directory with the MCMC simulation results. It should equal to the \code{output.dir} argument in \code{\link{run.e0hiv.mcmc}}.}
  \item{replace.output}{Logical. If \code{TRUE}, existing predictions in \code{output.dir} will be replaced by results of this run.}
  \item{predict.jmale}{Logical controlling if a joint female-male prediciton should be performed. This is done only if the underlying mcmcs in \code{sim.dir} correspond to a female simulation. In such a case the \code{\link[bayesLife]{e0.jmale.predict}} is invoked. Arguments to this function can be passed in \dots.}
  \item{nr.traj}{Number of trajectories to be generated. If \code{NULL}, the argument \code{thin} is taken to determine the number of trajectories. If both are \code{NULL}, the number of trajectories corresponds to the minimum of the size of the parameter sample and 2000.}
  \item{thin}{Thinning interval used for determining the number of trajectories. Only relevant, if \code{nr.traj} is \code{NULL}.}
  \item{burnin}{Number of iterations to be discarded from the beginning of the parameter traces.}
  \item{use.diagnostics}{Logical determining if an existing convergence diagnostics should be used for choosing the values of \code{thin} and \code{burnin}. In such a case, arguments \code{nr.traj}, \code{thin} and \code{burnin} are ignored. The \sQuote{best} values are chosen from results of running the \code{\link{e0.diagnose}} function. Only diagnostics can be used that suggest a convergence of the underlying MCMCs. If there are more than one such objects, the one is chosen whose recommendation for the number of trajectories is larger and closest to 2000.}
  \item{hiv.countries}{Codes of countries considered as epidemic in the projections. By default, countries are used that have \code{include_code = 3} in the \code{\link[bayesLife]{include}} dataset. All countries included in this argument have to have HIV prevalence and ART coverage projections available, either in the default \code{\link{HIVprevTrajectories}} and \code{\link{ARTcoverage}} datasets, or via the \code{my.art.file} and \code{my.hivtraj.file} arguments.}
  \item{my.art.file}{File name containing user-specified future values of ART coverage for countries with future epidemics. It should have the same format as the default \code{\link{ARTcoverage}} dataset.}
  \item{my.hivtraj.file}{File name containing probabilistic trajectories of HIV prevalence for countries considered as epidemic in the future. It should have the same format as the default \code{\link{HIVprevTrajectories}} dataset.}
  \item{scale.hivtraj}{Logical determining if the HIV trajectories should be scaled so that the median for each country aligns with the values in the \code{\link{HIVprevalence}} dataset (columns corresponding to future time periods). The scaling (which uses adjusted logit) is implemented in the function \code{\link{scale.hiv.trajectories}}.}
  \item{scale.hivtraj.tofile}{File name containing values to which the trajectories should be scaled to. It should have the same format as the \code{\link{HIVprevalence}} dataset.}
  \item{save.as.ascii}{Either a number determining how many trajectories should be converted into an ASCII file, or \dQuote{all} in which case all trajectories are converted. It should be set to 0, if no conversion is desired.}
  \item{start.year}{This argument should be only used if the start year of the prediction is before or at the present year of the MCMC run (see Details below). By default the prediction starts in the next time period after the present year (passed to \code{\link{run.e0hiv.mcmc}}).}
  \item{output.dir}{Directory into which the resulting prediction object and the trajectories are stored. If it is \code{NULL}, it is set to either \code{sim.dir}, or to \code{output.dir} of \code{mcmc.set$meta} if \code{mcmc.set} is given.}
  \item{low.memory}{Logical indicating if the prediction should run in a low-memory mode. If it is \code{FALSE}, the whole traces of all parameters, including the burnin, are loaded into memory. Otherwise, burnins are discarded and parameters are loaded as they are needed and are not kept in the memory.}
  \item{seed}{Seed of the random number generator. If \code{NULL} no seed is set. It can be used to generate reproducible projections.}
  \item{verbose}{Logical switching log messages on and off.}
  \item{\dots}{Additional arguments passed to the \code{\link[bayesLife]{e0.jmale.predict}} function.}
}
\details{
Works similarly as the \code{\link[bayesLife]{e0.predict}} function of the \pkg{\link[bayesLife]{bayesLife} package. For countries with future HIV epidemic, it uses trajectories of projected HIV prevelance as well as projected ART coverage to account for future HIV/AIDS epidemics, as described in Godwin and Raftery (2017).

See \code{\link[bayesLife]{e0.predict}} for more details.

}
\value{
	Object of class \code{\link[bayesLife]{bayesLife.prediction}}.
}
\references{
Godwin, J. and Raftery, A.E. (2017): Bayesian projection of life expectancy accounting for the HIV/AIDS epidemic. Demographic Research 37(48):1549–1610.

J. L. Chunn, A. E. Raftery, P. Gerland, H. Sevcikova (2013): Bayesian Probabilistic Projections of Life Expectancy for All Countries. Demography 50(3):777-801. <doi:10.1007/s13524-012-0193-x>
}
\author{
Hana Sevcikova, Jessica Godwin
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{run.e0hiv.mcmc}}, \code{\link[bayesLife]{e0.jmale.predict}}, \code{\link[bayesLife]{convert.e0.trajectories}}, \code{\link[bayesLife]{e0.predict}},
\code{\link[bayesLife]{get.e0.prediction}}, \code{\link[bayesLife]{summary.bayesLife.prediction}}
}
\examples{
\dontrun{
m <- run.e0hiv.mcmc(nr.chains = 1, iter = 50, thin = 1, verbose = TRUE)
pred <- e0hiv.predict(m, burnin = 25, verbose = TRUE)
summary(pred, country = "Botswana")}
hiv.countries.pred(m$meta)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{distribution}
\keyword{multivariate}
