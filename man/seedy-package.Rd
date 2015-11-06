\name{seedy-package}
\alias{seedy-package}
\alias{seedy}
\docType{package}
\title{
Simulation of Evolutionary and Epidemiological Dynamics
}
\description{
	A package of functions to simulate, visualize and assess epidemiological and pathogen genomic sample data collected during an outbreak.
}
\details{
\tabular{ll}{
Package: \tab seedy\cr
Type: \tab Package\cr
Version: \tab 1.3\cr
Date: \tab 2015-11-06\cr
License: \tab GPL-3\cr
}

}
\author{
Colin Worby (cworby@hsph.harvard.edu)
}
\examples{
# Load within host data
data(withinhost)
	
# Calculate genetic distance matrix
Gmat <- gd(withinhost$obs.strain, withinhost$libr, withinhost$nuc, 
           withinhost$librstrains)

# Set colors
colvec <- rainbow(1200)[1:1000] # Color palette
coltext <- rep("black", length(colvec)) # Corresponding text colors
coltext[680:970] <- "white" # White text for darker background colours

# Plot distance matrix
plotdistmat(Gmat, colvec, coltext, pos="bottomleft", labels=NULL, numbers=TRUE)
	
# Load outbreak data
data(outbreak)
sampledata <- outbreak$sampledata
epidata <- outbreak$epidata

# Calculate distance matrix for observed samples
distmat <- gd(sampledata[,3], outbreak$libr, outbreak$nuc, outbreak$librstrains)

# Now pick colors for sampled isolates
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
