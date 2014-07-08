\name{ancestors}
\alias{ancestors}

\title{
Vector of infection ancestors
}
\description{
Provides the chain of infection leading to a specified individual.
}
\usage{
ancestors(x, ID, sources)
}

\arguments{
  \item{x}{
  Individual ID.
}
  \item{ID}{
  List of person IDs.
}
  \item{sources}{
  List of infection sources, corresponding to \code{ID}.
}
}
\details{
  First element will be the ID \code{x}, last element will be zero.
}

\examples{
data(outbreak)
ancestors(9, outbreak$epidata[,1], outbreak$epidata[,4])
}

