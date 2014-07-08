\name{plotdiversity}
\alias{plotdiversity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Plot genetic diversity over time
}
\description{
	Calculates and plots the genetic diversity (mean pairwise number of SNPs between samples) in deep-sequenced samples.
}
\usage{
	plotdiversity(X, sample.times, makeplot = TRUE, filter=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
	Simulated deep sequence data over time.
}
  \item{sample.times}{
	Vector of sample collection times.
}
  \item{makeplot}{
	Should diversity samples be plotted?
}
  \item{filter}{
  	Should a subset of the sequence data be plotted? That is, should only the sequences at \code{sample.times} be plotted, out of a larger sample set?
  }
  \item{\dots}{
	Additional arguments to be passed to \code{plot} (if \code{makeplot=TRUE}).
}
}

\value{
	Returns a vector of diversity values for each sample time.
}

\examples{
data(deepseq)
D <- plotdiversity(deepseq, sample.times=(1:50)*1000, xlab="Time", 
                   ylab="Expected pairwise SNPs")

}

