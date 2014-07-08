\name{librtoDNA}
\alias{librtoDNA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Convert simulation objects to DNA sequences or Nexus files.
}
\description{
  Creates a character string or matrix of nucleotides (C, A, G, T), output to vector, matrix or Nexus file.
}
\usage{
librtoDNA(sampleID, libr, nuc, ref.strain, key, sampletime=NULL, 
          strings = FALSE, filename = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sampleID}{
  Vector of sample IDs to output.
}
  \item{libr}{
  Library object from simulation functions. A list in which each entry represents a unique genotype, and is a vector of mutated nucleotide positions relative to the reference sequence.
}
  \item{nuc}{
  Nucleotide database from simulation functions. A list (corresponding to \code{libr}) in which each entry represents a unique genotype, and is a vector of mutated nucleotides relative to the reference sequence.
}
  \item{ref.strain}{
  Reference strain to which the \code{libr} and \code{nuc} objects are compared (string of integers in 1,...,4).
}
  \item{key}{
  Vector of sample IDs corresponding to the order of \code{libr}.
}
  \item{sampletime}{
  Vector of sample times. If specified, incorporates sample times into genome name in Nexus file.
}
  \item{strings}{
  If \code{TRUE}, returns a character vector, each element containing one genotype. Otherwise, returns a (number of genotypes)x(length of genome) character matrix.
}
  \item{filename}{
  File to which sequence data should be written. Output format is Nexus. Not written out if \code{NULL}
}
}

\value{
  A character vector or matrix, depending on \code{strings}.
}
\examples{
data(withinhost)
G <- librtoDNA(sampleID=withinhost$obs.strain, libr=withinhost$libr, nuc=withinhost$nuc, 
               ref.strain=withinhost$ref.strain, key=withinhost$librstrains, strings=TRUE)
}

