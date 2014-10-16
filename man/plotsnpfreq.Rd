\name{plotsnpfreq}
\alias{plotsnpfreq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Plot frequency of polymorphic sites
}
\description{
	Returns the frequency of polymorphic sites above a probability level.
}
\usage{
plotsnpfreq(data, timepoint=1, type="S", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
	Full sequence data generated from the \code{\link{simulatepopulation}} function.
}
  \item{timepoint}{
	Which sampling time should be used.
}
  \item{type}{
  Type of plot desired. See \code{\link{plot}}.
}
  \item{...}{
	Additional arguments to be passed to \code{\link{plot}}.
}
}
\details{
  Plots frequency of polymorphic sites above each probability level.
}
\value{
	Returns a matrix with points used for the plot.
}

\examples{

data(deepseq)

# At the 25th time point
plotsnpfreq(data=deepseq, timepoint=25, xlab="Mutant frequency", ylim=c(0,25),
      ylab="No. sites", yaxs="i", xaxs="i", las=1, bty="l", col="red", lwd=2)
# At the 50th time point
plotsnpfreq(data=deepseq, timepoint=50, xlab="Mutant frequency", ylim=c(0,25),
      ylab="No. sites", yaxs="i", xaxs="i", las=1, bty="l", col="red", lwd=2)

}

