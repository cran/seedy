\name{estcoaltime}
\alias{estcoaltime}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate expected time to coalescence for sampled lineages in bottlenecked population
}
\description{
Estimates the expected time to coalescence for two randomly sampled lineages at a particular time, given the past population dynamics.
}
\usage{
estcoaltime(bottlesize, popsize, bottletimes, obstime)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bottlesize}{
  Effective size of population bottlenecks.
}
  \item{popsize}{
  Effective population size between bottlenecks.
}
  \item{bottletimes}{
  Vector of bottleneck times.
}
  \item{obstime}{
  Time of observation.
}
}
\details{
Expected coalescent time is calculated under the assumption that population remains constant at \code{popsize}, but drops to \code{bottlesize} for a single generation at bottleneck times. The probability of coalescence in a particular generation is the reciprocal of the population size at the time.
}
\value{
Returns the expected number of generations since coalescence.
}
\examples{
estcoaltime(bottlesize=10, popsize=10000, bottletimes=c(1000,2000,3000,4000), 
            obstime=5000)
}

