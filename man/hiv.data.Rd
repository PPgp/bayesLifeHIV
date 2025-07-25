\name{HIVprevalence}
\alias{HIVprevalence}
\alias{ARTcoverage}
\alias{HIVprevTrajectories}
\alias{HIVprevalence1y}
\alias{ARTcoverage1y}
\alias{HIVprevTrajectories1y}

\docType{data}
\title{
Datasets of HIV Prevalence and ART Coverage
}
\description{
Datasets of historical data on HIV prevalence, ART coverage, as well as future trajectories of HIV prevalence. Datasets with the suffix \dQuote{1y} are on the annual scale, otherwise on the 5-year scale.
}
\usage{
data(HIVprevalence)
data(HIVprevalence1y)
data(ARTcoverage)
data(ARTcoverage1y)
data(HIVprevTrajectories)
data(HIVprevTrajectories1y)
}

\format{
  Data frames \code{HIVprevalence*} and \code{ARTcoverage*} have one row per country. Columns include \code{country_code} (numeric), \code{name} (character), \code{include_code} (numeric) and one (numeric) column per time period, starting with \code{1970-1975} and ending with \code{2095-2100} for the 5-year versions, and from \code{1970} to \code{2100} for the annual versions. 
  
  Data frames \code{HIVprevTrajectories*} have columns \code{country_code} (numeric), \code{Trajectory} (numeric) and one (numeric) column per future time period, from 2010 to 2100. 
}
\details{
\describe{
\item{HIVprevalence, HIVprevalence1y}{Contains historical and projected HIV prevalence for countries with past or/and future epidemics. In the estimation, the columns for future years are only used for scaling. In addition, only countries are considered as epidemic in the estimation, if their \dQuote{include_code} is 1.}
\item{ARTcoverage, ARTcoverage1y}{Contains historical and projected values of ART coverage for countries with past or/and future epidemics. Countries considered as epidemic must have their \dQuote{include_code} set to 1 in both datasets, \code{HIVprevalence} and \code{ARTcoverage}. By default, all included countries are considered as epidemic. Note that in the estimation, only columns for historical time periods are used, while future time periods are used in the projections.}
\item{HIVprevTrajectories, HIVprevTrajectories1y}{These are probabilistic trajectories of HIV prevalence for countries considered as epidemic in the future. The default dataset contains 1000 trajectories, scaled to the HIVprevalence dataset.}
}
}
%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}

\examples{
data(HIVprevalence)
head(HIVprevalence)
}
\keyword{datasets}
