\name{gd}
\alias{gd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Genetic distance matrix
}
\description{
Given a set of genomic samples, returns a pairwise geentic distance matrix
}
\usage{
gd(samp, libr, nuc, key)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{samp}{
  Vector of sample IDs.
}
  \item{libr}{
  Library object from simulation functions. A list in which each entry represents a unique genotype, and is a vector of mutated nucleotide positions relative to the reference seqeunce.
}
  \item{nuc}{
  Nucleotide database from simulation functions. A list (corresponding to \code{libr}) in which each entry represents a unique genotype, and is a vector of mutated nucleotides relative to the reference sequence.
}
  \item{key}{
  Vector of sample IDs in the order they appear in the \code{libr} and \code{nuc} objects.
}
}
\details{
  Each element of \code{samp} represents one row and column in the genetic distance matrix.
}
\value{
  Returns a symmetric genetic distance matrix with rows and columns corresponding to the \code{samp} vector.
}

\examples{
data(withinhost)
gd(withinhost$obs.strain, withinhost$libr, withinhost$nuc, withinhost$librstrains)
}

