\name{plotoutbreak}
\alias{plotoutbreak}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot outbreak
}
\description{
Provides a graphical representation of simulated outbreak and sampled genomes.
}
\usage{
plotoutbreak(epidata, sampledata, col = "red", stack = TRUE, arr.len = 0.1, 
blockheight = 0.5, hspace = max(epidata[,3])/20, labels = NULL, label.pos = "left", 
	block.col = "grey", jitter = 0, pch = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{epidata}{
	Simulated epidemiological data - matrix consisting of person IDs, infection and recovery times, and infection source.
}
  \item{sampledata}{
	Simulated genomic sample data - matrix of person IDs, sample times and genome ID. Genomic sample data simulated with \code{full=TRUE} cannot be used.
}
  \item{col}{
	A vector of colors to represent increasing genetic distance.
}
  \item{stack}{
	Should infectious periods be organized to minimize plot height? Alternatively, each period occupies one row.
}
  \item{arr.len}{
	Arrow length. See \code{\link{arrows}}.
}
  \item{blockheight}{
	The height of each bar representing infectious periods. Takes a value between 0 and 1.
}
  \item{hspace}{
	Minimum horizontal space between two infectious period bars (in units of time).
}
  \item{labels}{
	Labels for each infectious episode. First columns of \code{epidata} by default.
}
  \item{label.pos}{
	Position of labels. Accepted values are "centre", "left" and "right".
}
  \item{block.col}{
	Background color for each infectious period bar.
}
  \item{jitter}{
	Amount of jitter to be applied to genome sample points (as a proportion of plot dimensions).
}
  \item{pch}{
	Point type for genome samples (see \code{\link{par}}).
}
  \item{\dots}{
	Additional arguments to be passed to \code{\link{plot}}.
}
}
\details{
	Graphical representation of transmission dynamics and sampled genomes. For multiple genome samples per time point, set \code{jitter>0}.
}

\examples{
data(outbreak)
sampledata <- outbreak$sampledata
epidata <- outbreak$epidata

distmat <- gd(sampledata[,3], outbreak$libr, outbreak$nuc, outbreak$librstrains)

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
