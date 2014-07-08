\name{simulatebacteria}
\alias{simulatebacteria}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate a bacterial population
}
\description{
Simulates a bacterial population undergoing a Wright-Fisher discrete-time evolutionary process.
}
\usage{
simulatebacteria(eq.size, m.rate, runtime, sample.times, n.samples=1, 
genomelength=100000, bottle.times=0, bottle.size=1, fullrecords=FALSE, 
deepseq=FALSE, feedback=1000, init.freq=1, libr=NULL, nuc=NULL, ref.strain=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eq.size}{
	Equilibrium effective population size.
}
  \item{m.rate}{
	Mutation rate (per sequence per generation).
}
  \item{runtime}{
	Number of generations for simulation to run.
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
  \item{bottle.times}{
	Vector of population bottleneck times.
}
  \item{bottle.size}{
	Size of population bottleneck (if \code{bottle.times!=0}).
}
  \item{fullrecords}{
	Should type and frequency of current genomes be recorded at every generation?
}
  \item{deepseq}{
	Should deep sequenced sampled be recorded? All present genotypes and their frequencies are recorded for each sampling time.
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
}
\details{
	Simulation of a bacterial population. Population is by default initially clonal, and initiated by a single genotype. Population is prevented from extinction at all times. Population can be passed through repeated bottlenecks using the arguments \code{bottle.times} and \code{bottle.size}.
}
\value{
	Returns a list of sampling data;
  	\item{libr}{List of unique genotypes observed. Each entry is a vector of mutant loci relative to the reference strain.}
	\item{nuc}{List of mutation types corresponding to \code{libr}. Each entry is a vector of nucleotides mutated from the reference strain, corresponding to entries in \code{libr}.}
  	\item{librstrains}{Vector of unique genotype IDs corresponding to the \code{libr} object.}
  	\item{obs.strain}{If \code{deepseq=TRUE}, list of observed genotypes, each entry corresponding to the \code{sample.times}. If \code{deepseq=FALSE}, a vector of genotype IDs, corresponding to the returned vector \code{obs.time}.}
  	\item{obs.freq}{List of observed genotype frequencies returned if \code{deepseq=TRUE}.}
  	\item{obs.time}{Vector of observation times returned if \code{deepseq=FALSE}.}
  	\item{obs.time}{Vector of observation times returned if \code{deepseq=FALSE}.}
  	\item{ref.strain}{Reference strain used.}  	
}

\examples{
# Generate 5 genome samples at 5 time points
X <- simulatebacteria(eq.size=2000, m.rate=0.0005, runtime=10000, 
                      sample.times=(1:5)*2000, n.samples=5, 
                      genomelength=10000, deepseq=FALSE)
	
# Generate deep-sequenced samples at 5 time points
X <- simulatebacteria(eq.size=2000, m.rate=0.0005, runtime=10000, 
                      sample.times=(1:5)*2000, genomelength=10000, 
                      bottle.times=5000, bottle.size=1, deepseq=TRUE)

}

