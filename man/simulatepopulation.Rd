\name{simulatepopulation}
\alias{simulatepopulation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate a pathogen population
}
\description{
Simulates a pathogen population undergoing a Wright-Fisher discrete-time evolutionary process.
}
\usage{
simulatepopulation(m.rate, runtime, equi.pop, sample.times, n.samples=1, 
genomelength=100000, shape=flat, bottle.times=0, bottle.size=1, full=FALSE, 
feedback=1000, init.freq=1, libr=NULL, nuc=NULL, ref.strain=NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m.rate}{
	Mutation rate (per sequence per generation).
}
  \item{runtime}{
	Number of generations for simulation to run.
}
  \item{equi.pop}{
  Equilibrium effective population size of pathogens within-host.
}
  \item{sample.times}{
	Vector of times at which to sample population.
}
  \item{n.samples}{
	Number of samples to take at each sampling point (if \code{deepseq=FALSE}).
}
  \item{genomelength}{
	Length of genome.
}
  \item{shape}{
  Function describing the population growth dynamics. See Details.
  }
  \item{bottle.times}{
	Vector of population bottleneck times.
}
  \item{bottle.size}{
	Size of population bottleneck (if \code{bottle.times!=0}).
}
  \item{full}{
	Should complete samples (all genotypes and their frequencies) be returned?
}
  \item{feedback}{
	Intervals between R feedback on simulation progress.
}
  \item{init.freq}{
	Vector of initial frequencies of genotypes, if initial population should be specified. By default, the population grows from a single, randomly generated genotype.
}
  \item{libr}{
	Initial list of genotypes, if initial population should be specified. Must have same length as \code{init.freq}.
}
  \item{nuc}{
	Initial list of mutations, if initial population should be specified. Must have same length as \code{init.freq}.	
}
  \item{ref.strain}{
	Reference strain. By default, this is randomly generated.
}
  \item{...}{
  Additional arguments to be passed to the \code{shape} function.
}

}
\details{
	Simulation of a bacterial population. Population is by default initially clonal, and initiated by a single genotype. Population is prevented from extinction at all times. Population can be passed through repeated bottlenecks using the arguments \code{bottle.times} and \code{bottle.size}.
Population growth dynamics are defined by the function called by 'shape'. This function returns the expected population size at each time step, given the total simulation time. By default, the population is expected to grow exponentially until reaching an equilibrium level, specified by \code{equi.pop} (\code{\link{flat}}). Alternatively, the population can follow a sinusoidal growth curve, peaking at \code{runtime}/2 (\code{\link{hump}}). User-defined functions should be of the form \code{function(time,span,equi.pop,...)}, where \code{span} is equal to \code{runtime} in this setting.
}
\value{
	Returns a list of sampling data;
  	\item{libr}{List of unique genotypes observed. Each entry is a vector of mutant loci relative to the reference strain.}
	\item{nuc}{List of mutation types corresponding to \code{libr}. Each entry is a vector of nucleotides mutated from the reference strain, corresponding to entries in \code{libr}.}
  	\item{librstrains}{Vector of unique genotype IDs corresponding to the \code{libr} object.}
  	\item{obs.strain}{If \code{full=TRUE}, list of observed genotypes, each entry corresponding to the \code{sample.times}. If \code{full=FALSE}, a vector of genotype IDs, corresponding to the returned vector \code{obs.time}.}
  	\item{obs.freq}{List of observed genotype frequencies returned if \code{full=TRUE}.}
  	\item{obs.time}{Vector of observation times returned if \code{full=FALSE}.}
  	\item{ref.strain}{Reference strain used.}  	
}

\examples{
# Generate 5 genome samples at 5 time points
X <- simulatepopulation(m.rate=0.0005, runtime=10000, equi.pop=2000,
                      sample.times=(1:5)*2000, n.samples=5, 
                      genomelength=10000, full=FALSE)
	
# Generate complete observations at 5 time points
X <- simulatepopulation( m.rate=0.0005, runtime=10000, equi.pop=2000,
                      sample.times=(1:5)*2000, genomelength=10000, 
                      bottle.times=5000, bottle.size=1, full=TRUE)
}

