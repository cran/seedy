\name{hump}
\alias{hump}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sinusoidal population growth model
}
\description{
Provides the expected pathogen population at any point during infection
}
\usage{
hump(x, span, max.size)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  Time after infection.
}
  \item{span}{
  Total duration of infection.
}
  \item{max.size}{
  Expected maximum population size, attained at midpoint of infection.
}

}

\value{
  Returns expected population size at \code{x}.
}

\examples{
hump(1:200, span=200, max.size=1000)
}

