\name{expsnps}
\alias{expsnps}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Distribution of genetic distances
}
\description{
For any pair of epidemiologically linked individuals, returns the distribution of genetic distance separating randomly drawn samples.
}
\usage{
expsnps(x, m.rate, c.rate, tau)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  Vector of (non-negative integer) quantiles.
}
  \item{m.rate}{
  Mutation rate (per genome per generation).
}
  \item{c.rate}{
  Rate of coalescence prior to lineage divergence (constant).
}
  \item{tau}{
  Total time from lineage divergence to observations.
}
}
\details{
  A pair of genomes sampled during an outbreak will have an epidemiological and an evolutionary relationship with each other. 'Lineage divergence' is defined to be the time at which the lineages ceased to exist within the same host (or, the latest possible time of coalescence). Mutations may arise in two distinct time periods: (a) the time between coalescence and lineage divergence, and (b) the time between lineage divergence and observation. The latter follows a Poisson distribution with mean equal to the mutation rate multiplied by the total time between lineage divergence and observations (two branches). The former is not Poisson distributed, since the time to coalescence is typically unknown (but follows an exponential distribution with a constant effective population size). It follows that, with a constant effective population size (and therefore coalescent rate), the number of mutations follows a geometric distribution. Therefore, the total number of SNPs between two samples is a geometric-Poisson mixture distribution, which this function returns.
}
\value{
Returns the probability density for given genetic distances.
}

\seealso{
\link{estcoaltime}
}
\examples{
expsnps(3, m.rate=0.003, c.rate=1/10000, tau=1000)
  
plot(expsnps(0:100, m.rate=0.003, c.rate=1/2000, tau=5000), type="h")
}

