\name{hiv.countries.pred}
\alias{hiv.countries.pred}
\alias{hiv.countries.est}

\title{
Information about HIV countries
}
\description{
Functions for returning information about which countries are used when estimating the HIV/AIDS model and which countries are considered as epidemic when projecting life expectancy.
}
\usage{
hiv.countries.est(meta)

hiv.countries.pred(meta)
}

\arguments{
  \item{meta}{Object of class \code{\link[bayesLife]{bayesLife.mcmc.meta}}.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
Data frame with columns \dQuote{country_code} and \dQuote{country_name} containing 
countries used in the estimation (\code{hiv.countries.est}) or prediction (\code{hiv.countries.pred})
}

\author{
Hana Sevcikova
}

\seealso{
\code{\link{run.e0hiv.mcmc}}, \code{\link{e0hiv.predict}}
}
\examples{
\dontrun{
sim.dir <- tempfile()
m <- run.e0hiv.mcmc(nr.chains = 1, iter = 50, thin = 1, 
    output.dir = sim.dir, verbose = TRUE)
hiv.countries.est(m$meta)
hiv.countries.pred(m$meta)

unlink(sim.dir, recursive=TRUE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ attribute }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
