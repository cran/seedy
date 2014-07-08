\name{plotdistmat}
\alias{plotdistmat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Plot genetic distance matrix
}
\description{
Provides a graphical representation of the pairwise genetic distance matrix for a collection of genomes.
}
\usage{
plotdistmat(distmat, colvec, coltext, pos = "topleft", labels = NULL, 
            numbers = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{distmat}{
	Symmetrical pairwise distance matrix, the [i,j]th entry corresponding to the genetic distance between genomes i and j.
}
  \item{colvec}{
	Vector of colors to represent increasing genetic distance.
}
  \item{coltext}{
	Vector of colors for numerals representing geentic distance on plot (if \code{numbers=TRUE}).
}
  \item{pos}{
	Position of the matrix in plot. Allowed values are \code{"topleft"}, \code{"topright"}, \code{"bottomleft"} and \code{"bottomright"}.
}
  \item{labels}{
	Axis labels for genomes (by default, marked 1,...,n).
}
  \item{numbers}{
	Should the genetic distance be recorded in each cell?
}
  \item{\dots}{
	Additional arguments to be passed to \code{plot}.
}
}
\details{
	Plots the upper (or lower) diagonal genetic distance matrix, with each entry colored according to geentic distance.
}
\seealso{
	\code{\link{gd}}
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

