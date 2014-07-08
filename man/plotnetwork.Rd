\name{plotnetwork}
\alias{plotnetwork}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Plot weighted transmission network
}
\description{
For a given weighted adjacency matrix, plots a directed network illustrating transmission probabilities.
}
\usage{
plotnetwork(probmat, labels, arrlen = 0.15, scale = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{probmat}{
	Transmission probability matrix. The [i,j]th element represents the probability that person i infected person j.
}
  \item{labels}{
	Labels representing row and columns in the probability matrix.
}
  \item{arrlen}{
	Length of arrow heads (see \code{\link{arrows}}).
}
  \item{scale}{
	Scaling for colors according to probability values.
}
  \item{\dots}{
	Additional arguments to be passed to plot.
}
}

\examples{
data(outbreak)
K <- networkmat(outbreak$epidata[,1], outbreak$epidata[,4])
plotnetwork(K, labels=outbreak$epidata[,1])
		
}
