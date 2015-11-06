\name{deepseqmat}
\alias{deepseqmat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Generate deep sequence observations
}
\description{
  Generates a matrix of polymorphism frequencies
}
\usage{
  deepseqmat(X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
    Simulated complete genomic sampling data from the \code{\link{simulateoutbreak}} or \code{\link{simulatepopulation}} functions.
  }
}
\details{
  The argument \code{full=TRUE} must be passed to \code{\link{simulateoutbreak}} or \code{\link{simulatepopulation}} in order for complete genomic data to be generated for each sample. Data of this format can be passed to this function.
}
\value{
Returns a matrix, each row represents a polymorphic site, each column a sample. Frequencies of the mutant type (relative to the initial strain) are reported.
}
\examples{
data(deepseq)
W <- deepseqmat(deepseq)
}