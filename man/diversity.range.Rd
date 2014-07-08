\name{diversity.range}
\alias{diversity.range}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Range of genetic diversity over time
}
\description{
Generates multiple populations stochastically from an identical source, and measures the resulting diversity over time in each.
}
\usage{
diversity.range(eq.size, m.rate, runtime, iterations = 10, n.points = 100, 
genomelength = 1e+05, bottle.times=0, bottle.size=1, feedback = 1000, 
makeplot = TRUE, area = TRUE, colline = "blue", colarea = rgb(0, 0, 1, 0.4), 
ref.strain = NULL, init.freq = 1, libr=NULL, nuc=NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eq.size}{
  Equilibrium effective population size.
}
  \item{m.rate}{
  Mutation rate (per genome per generation).
}
  \item{runtime}{
  Number of bacterial generations over which to simulate.
}
  \item{iterations}{
  Number of populations to simulate.
}
  \item{n.points}{
  Number of equidistant points to sample diversity during \code{runtime}.
}
  \item{genomelength}{
  Genome length.
}
  \item{bottle.times}{
    Vector of population bottleneck times.
  }
  \item{bottle.size}{
    Size of population bottleneck.
  }
  \item{feedback}{
  Number of generations between each simulation report.
}
  \item{makeplot}{
  Should resulting diversity be plotted?
}
  \item{area}{
  Should 95\% central quantile of genetic diversity be shaded? If \code{FALSE}, then individual diversity trajectories are plotted.
}
  \item{colline}{
  Colour of lines (if \code{makeplot=TRUE}).
}
  \item{colarea}{
  Colour of shaded area (if \code{makeplot=TRUE} and \code{area=TRUE}).
}
  \item{ref.strain}{
  Reference strain, if required.
}
  \item{init.freq}{
  Initial frequency of strains in starting population (if \code{libr} and \code{nuc} specified)
}
  \item{libr}{
  Library of initial sequences.
}
  \item{nuc}{
  Nucleotides at polymorphic sites, corresponding to \code{libr}.
}

  \item{\dots}{
  Additional arguments to be passed to \code{plot}.
}
}
\details{
  Provides an empirical estimate of the expected genetic diversity (pairwise SNP distance) over time, with associated uncertainty. Initial population can be specified using the \code{libr}, \code{nuc} and \code{init.freq} arguments, otherwise population is grown from a single genotype. Resolution can be improved by increasing \code{n.points}, and accuracy by increasing \code{iterations} (at the expense of accuracy).
}
\value{
A \code{iterations} by \code{n.points} matrix with diversity over time for each simulation.
}

\seealso{
 \code{\link{plotdiversity}}
}
\examples{
iterations <- 10
K <- diversity.range(eq.size=1000, m.rate=0.0005, runtime=1000, 
iterations=iterations, n.points=100, genomelength=100000, feedback=100, 
makeplot=TRUE, area=TRUE, colline="blue", colarea=rgb(0,0,1,0.4))
}

