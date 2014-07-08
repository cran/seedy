\name{networkmat}
\alias{networkmat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Create adjacency matrix
}
\description{
For a given set of infection routes, returns an adjacency matrix.
}
\usage{
networkmat(ID, sources)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ID}{
	Vector of infective IDs
}
  \item{sources}{
	Vector of infection sources corresponding to \code{ID}.
}
}

\value{
	Returns a matrix with the [i,j]th entry equal to 1 if the ith infective in \code{ID} infected the jth infective in \code{ID}.
}

\seealso{
\code{\link{ancestors}}
}
\examples{
data(outbreak)
networkmat(outbreak$epidata[,1], outbreak$epidata[,4])
}

