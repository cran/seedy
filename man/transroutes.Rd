\name{transroutes}
\alias{transroutes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Assessment of transmission routes using theoretical SNP distribution
}
\description{
Calculates likelihood and posterior probability of each potential transmission route using a geometric-Poisson approximation of SNP distance.
}
\usage{
transroutes(ID, GD, sample.times, inf.times, rec.times=NULL, mut.rate, eq.size, 
            bottle.size=1, p.level=0.95, geninterval=NULL, summary=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ID}{
  Vector of individual IDs.
}
  \item{GD}{
  Matrix of genetic distances, with the [i,j]-th entry corresponding to the genetic distance between samples from the ith and jth entry of \code{ID}.
}
  \item{sample.times}{
  Vector of genome sampling times (same length as \code{ID}).
}
  \item{inf.times}{
  Vector of infection times (same length as \code{ID}).
}
  \item{rec.times}{
  Vector of recovery times (same length as \code{ID} if specified).
}
  \item{mut.rate}{
  Mutation rate.
}
  \item{eq.size}{
  Equilibrium within-host effective population size.
}
  \item{bottle.size}{
    Size of population bottleneck.
  }
  \item{p.level}{
  Probability level at which to reject potential transmission routes.
}
  \item{geninterval}{
  Generation interval (if \code{bottle.size}>1).
}
  \item{summary}{
  Should a summary for each ID be printed to screen?
}
}
\details{
  Calculates the likelihood and posterior probability (given a flat prior) for each potential transmission route, as well as indicating which transmission routes would be rejected at a given probability level. If recovery times are not specified, each individual is assumed to be infectious for the duration of the outbreak.
}
\value{
  Returns calculated values:
  	\item{maxpostsource}{Vector of maximum posterior probability transmission sources corresponding to \code{ID}.}
  	\item{likelihood}{A matrix of likelihood values from the geometric-Poisson distribution. The [i,j]th entry provides the likelihood that the i-th individual infected the j-th individual.}
  	\item{posterior}{A matrix of posterior transmission probabilities. The [i,j]th entry provides the posterior probability that the i-th individual infected the j-th individual.}
  	\item{closestsource}{A list of the individuals carrying the most genetically similar genotype. Each entry corresponds to \code{ID}.}
  	\item{reject}{Matrix indicating whether a transmission route would be rejected at the specified probability level. The [i,j]th entry is equal to 1 if the route from i to j is rejected at this level.}
}

\seealso{
 \code{\link{expsnps}}
}
\examples{

  data(transmission)
  W <- transmission
  ID <- unique(W$sampledata[,1])
  GD <- gd(W$sampledata[,3], W$libr, W$nuc, W$librstrains)
  
  sample.times <- W$sampledata[,2]
  inf.times <- numeric(length(ID))
  rec.times <- numeric(length(ID))
  truesource <- numeric(length(ID))
  for (i in 1:length(ID)) {
    inf.times[i] <- W$epidata[which(W$epidata[,1]==ID[i]),2]
    rec.times[i] <- W$epidata[which(W$epidata[,1]==ID[i]),3]
    truesource[i] <- W$epidata[which(W$epidata[,1]==ID[i]),4]
  }
  
  K <- transroutes(ID=ID, GD=GD, sample.times=sample.times, inf.times=inf.times, 
                   rec.times=rec.times, mut.rate=0.01, eq.size=5000, 
                   bottle.size=1, p.level=0.95, summary=TRUE)
  truesource
  K$maxpostsource

}

