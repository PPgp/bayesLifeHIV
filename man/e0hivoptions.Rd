\name{using.bayesLifeHIV}
\alias{using.bayesLifeHIV}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Setting global options
}
\description{
Setting and updating global options.
}
\usage{
using.bayesLifeHIV()
}

\details{
Function \code{using.bayesLifeHIV} sets all global options to their default values. It overwrites options in the object \code{\link[bayesLife]{e0options}} set by  \code{\link[bayesLife]{using.bayesLife}}. Among others, it adds the parameter \code{betanonART} and its initial value to the set of MCMC parameters. 

Since \pkg{bayesLifeHIV} shares (and overwrites) some of the global objects of \pkg{bayesLife}, it is recommended not to mix simulations from \pkg{bayesLifeHIV} and \pkg{bayesLife} in the same session. Or, when working on objects created with \pkg{bayesLifeHIV}, one should start with calling \code{using.bayesLifeHIV()} and when switching to a \pkg{bayesLife} simulation, preceed it with \code{using.bayesLife()}.

One can view or modify the global options via any of the functions described in \code{\link[bayesLife]{e0options}}, such as \code{e0options()}, \code{e0mcmc.options()} or \code{e0pred.options()}.
}


\author{
Hana Sevcikova
}


\seealso{
\code{\link[bayesLife]{e0options}}
}
\examples{
# view default options for an annual simulation
e0mcmc.options(annual = TRUE) 

# Set new initial value for betanonART and modify ranges for outliers 
e0mcmc.options(betanonART = list(ini = -0.3), outliers = c(-7, 7),
                    annual = TRUE)
e0mcmc.options("betanonART", "outliers", annual = TRUE)

# revert to defaults
using.bayesLifeHIV()
e0mcmc.options("betanonART", "outliers", annual = TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~manip }% use one of  RShowDoc("KEYWORDS")

