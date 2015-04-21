\name{flat}
\alias{flat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Equilibrium population growth model
}
\description{
Provides the expected pathogen population at any point during infection
}
\usage{
flat(x, span, eq.size)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  Time after infection.
}
  \item{span}{
  Total duration of infection.
}
  \item{eq.size}{
  Expected equilibrium population size.
}

}

\value{
  Returns expected population size at \code{x}.
}

\examples{
flat(1:200, span=200, eq.size=1000)
}

