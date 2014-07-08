\name{withinhost}
\alias{withinhost}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Within-host genomic data
}
\description{
	Simulated within-host bacterial genomic samples.
}
\usage{
withinhost
}
%- maybe also 'usage' for other objects documented here.
\format{
  Two bacterial samples were taken every 2000 generations up to time 10000.
}

\examples{
data(withinhost)
Gmat <- gd(withinhost$obs.strain, withinhost$libr, withinhost$nuc, 
           withinhost$librstrains)

colvec <- rainbow(1200)[1:1000] # Color palette
coltext <- rep("black", length(colvec)) # Corresponding text colors
coltext[680:970] <- "white" # White text for darker background colours

plotdistmat(Gmat, colvec, coltext, pos="bottomleft", labels=NULL, numbers=TRUE)
	
}

