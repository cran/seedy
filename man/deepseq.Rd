\name{deepseq}
\alias{deepseq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Deep-sequenced bacterial samples
}
\description{
	Simulated deep-sequence samples within a host.
}
\usage{
	deepseq
}
%- maybe also 'usage' for other objects documented here.
\format{
  Deep-sequenced samples were taken at 50 time points. 
}

\examples{
data(deepseq)
D <- plotdiversity(deepseq, sample.times=(1:50)*1000, xlab="Time", ylab="Expected pairwise SNPs")

		
}
