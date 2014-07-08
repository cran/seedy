\name{simulateoutbreak}
\alias{simulateoutbreak}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Simulate transmission and evolutionary dynamics
}
\description{
	Simulate within-host evolutionary dynamics on top of an SIR transmission process and generate genomic samples.
}
\usage{
simulateoutbreak(init.sus, inf.rate, rem.rate, mut.rate, nmat = NULL, 
                 equi.pop = 10000, init.inf = 1, inoc.size = 1, 
                 samples.per.time = 1, samp.schedule = "random", 
                 samp.freq = 500, mincases = 1, feedback = 500, 
                 glen = 1e+05, ref.strain = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{init.sus}{
	Initial number of susceptible individuals.
}
  \item{inf.rate}{
	SIR rate of infection.
}
  \item{rem.rate}{
	SIR rate of removal.
}
  \item{mut.rate}{
	Mutation rate (per genome per generation).
}
  \item{nmat}{
	Connectivity matrix. Entry [i,j] gives the relative rate at which person i may contact person j.
}
  \item{equi.pop}{
	Equilibrium effective population size of pathogens within-host.
}
  \item{init.inf}{
	Initial number of infected individuals.
}
  \item{inoc.size}{
	Size of pathogen inoculum at transmission.
}
  \item{samples.per.time}{
	Number of samples taken at sampling times.
}
  \item{samp.schedule}{
	How should sampling be conducted? Accepted values are: "calendar" - samples are taken from all current infectives every \code{samp.freq} generations; "individual" - samples are taken from each infective at intervals of \code{samp.freq} after infection; "random" - samples are taken at one time between infection and removal for each infective.
}
  \item{samp.freq}{
	Number of generations between each sampling time (see \code{samp.schedule}).
}
  \item{mincases}{
	Minimum final size of outbreak to output. If final size is less than this value, another outbreak is simulated.
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

W <- simulateoutbreak(init.sus=10, inf.rate=0.002, rem.rate=0.001, mut.rate=0.0001, 
                      equi.pop=2000, inoc.size=1, samples.per.time=10, 
                      samp.schedule="calendar", samp.freq=500, mincases=3) 
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
             blockheight=0.5, hspace=500, label.pos="left", block.col="grey", 
             jitter=0.004, xlab="Time", pch=1) 
}

