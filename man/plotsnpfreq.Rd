\name{plotsnpfreq}
\alias{plotsnpfreq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Plot expected frequency of polymorphic sites in a model deep sequencing project
}
\description{
	Given an average depth of coverage and per base sequencing error rate, estimate the read frequency of most common polymorphic sites.
}
\usage{
plotsnpfreq(data, timepoint=1, coverage=50, error=0.001, iterations=100, maxsnp=50, 
            legend=TRUE, ylim=c(0,1.5*coverage), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
	Full sequence data generated from the \code{link{simulatepopulation}} function.
}
  \item{timepoint}{
	Which sampling time should be used.
}
  \item{coverage}{
	Coverage depth.
}
  \item{error}{
	Sequencing erorr rate per base.
}
  \item{iterations}{
	Number of iterations to generate confidence bounds.
}
  \item{maxsnp}{
	Number of polymorphic sites to show.
}
  \item{legend}{
	Should a legend be plotted?
}
  \item{ylim}{
	Bounds of y axis.
}
  \item{...}{
	Additional arguments to be passed to \code{plot}.
}
}
\details{
  Plots the expected number of reads containing the most frequent polymorphisms, with a 95\% confidence interval, derived from repeated random draws. Additionally shows the frequency of false positive observations arising through sequencing error.
}
\value{
	Returns an \code{iterations} x \code{maxsnp} matrix. Each column represents a polymorphic site, and each row one iteration. Each entry provides the number of reads containing the polymorphism in a given iteration. Columns are ordered by frequency.
}

\author{T. D. Read (tread@emory.edu)}

\examples{

data(deepseq)

# At the 25th time point
plotsnpfreq(data=deepseq, timepoint=25, xaxt="n", xlab="Ranked polymorphic sites", 
      ylab="Reads", yaxs="i", las=1)
# At the 50th time point
plotsnpfreq(data=deepseq, timepoint=50, xaxt="n", xlab="Ranked polymorphic sites", 
      ylab="Reads", yaxs="i", las=1)

}

