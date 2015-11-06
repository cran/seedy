\name{simfixoutbreak}
\alias{simfixoutbreak}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Simulate evolutionary dynamics on a given transmission tree
}
\description{
  Simulate within-host evolutionary dynamics on top of an existing transmission tree and generate genomic samples.
}
\usage{
simfixoutbreak(ID,inf.times, rec.times, inf.source, mut.rate, equi.pop=10000, shape=flat,
           inoc.size=1, imp.var=25, samples.per.time=1, samp.schedule="random", 
           samp.freq=500, full=FALSE, feedback=500, glen=100000, 
           ref.strain=NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ID}{
	Vector of unique IDs.
}
  \item{inf.times}{
	Vector of (integer) infection times.
}
  \item{rec.times}{
	Vector of (integer) removal times.
}
  \item{inf.source}{
	Vector of infection sources. The \code{i}th entry corresponds to the ID of the source of infection. For importations, the source should be 0.
}
  \item{mut.rate}{
	Mutation rate (per genome per generation).
}
  \item{equi.pop}{
	Equilibrium effective population size of pathogens within-host.
}
  \item{shape}{
  Function describing the population growth dynamics. See Details.
  }
  \item{inoc.size}{
	Size of pathogen inoculum at transmission.
}
  \item{imp.var}{
  The expected number of mutations separating unconnected importations.
}
  \item{samples.per.time}{
	Number of samples taken at sampling times.
}
  \item{samp.schedule}{
	How should sampling be conducted? Accepted values are: "calendar" - samples are taken from all current infectives every \code{samp.freq} generations; "individual" - samples are taken from each infective at intervals of \code{samp.freq} after infection; "random" - samples are taken at one time between infection and removal for each infective.
}
  \item{full}{Should `full' genomic sampling be returned? That is, should a vector of genotypes and their respective frequencies be stored from each individual's sampling times?
}
  \item{samp.freq}{
  Number of generations between each sampling time (see \code{samp.schedule}).
}
  \item{feedback}{
	Number of generations between simulation updates returned to R interface.
}
  \item{glen}{
	Length of genome.
}
  \item{ref.strain}{
	Initial sequence. By default, a random sequence of length \code{glen}.
}
  \item{...}{
  Additional arguments to be passed to the \code{shape} function.
}
}
\details{
Population growth dynamics are defined by the function called by 'shape'. This function returns the expected population size at each time step, given the total simulation time. By default, the population is expected to grow exponentially until reaching an equilibrium level, specified by \code{equi.pop} (\code{\link{flat}}). Alternatively, the population can follow a sinusoidal growth curve, peaking at \code{runtime}/2 (\code{\link{hump}}). User-defined functions should be of the form \code{function(time,span,equi.pop,...)}, where \code{span} is equal to the duration of infection in this setting.
}
\value{
	Returns a list of outbreak data:
	\item{epidata}{A matrix of epidemiological data with columns: person ID, infection time, removal time, source of infection.}
	\item{sampledata}{A matrix of genome samples with columns: person ID, sampling time, genome ID.}
	\item{libr}{A list with an entry for each unique genotype observed. Each entry is a vector of mutation positions relative to the reference genome.}
	\item{nuc}{A list with an entry for each unique genotype observed. Each entry is a vector of nucleotide types (integer between 1 and 4).}
	\item{librstrains}{A vector of unique genotype IDs corresponding to the \code{libr} object.}
	\item{endtime}{End time of the outbreak.}	
}

\examples{

# Simulate a transmission chain
inf.times <- (0:20)*100
rec.times <- inf.times + 100 + rpois(21,50)
inf.source <- 0:20
inf.source[c(3,11)] <- 0 # Two importations
mut.rate <- 0.001

# Now simulate evolutionary dynamics and samples on top of this tree
W <- simfixoutbreak(ID=1:21, inf.times, rec.times, inf.source, mut.rate, equi.pop=1000, shape=flat,
                    inoc.size=10, imp.var=25, samples.per.time=5, samp.schedule="random", 
                    samp.freq=500, full=FALSE, feedback=100, glen=100000, 
                    ref.strain=NULL)

sampledata <- W$sampledata
epidata <- W$epidata

# Calculate distance matrix for observed samples
distmat <- gd(sampledata[,3], W$libr, W$nuc, W$librstrains)

# Now pick colors for sampled isolates
colvec <- rainbow(1200)[1:1000] # Color palette
refnode <- 1 # Compare distance to which isolate?
colv <- NULL # Vector of colors for samples
maxD <- max(distmat[,refnode])

for (i in 1:nrow(sampledata)) {
  colv <- c(colv, 
            colvec[floor((length(colvec)-1)*(distmat[refnode,i])/maxD)+1])
}


plotoutbreak(epidata, sampledata, col=colv, stack=TRUE, arr.len=0.1, 
             blockheight=0.5, hspace=60, label.pos="left", block.col="grey", 
             jitter=0.004, xlab="Time", pch=1) 

}

