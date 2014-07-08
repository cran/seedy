\name{meansnps}
\alias{meansnps}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Mean diversity within a single population
}
\description{
  Calculates the expected pairwise genetic distance (number of SNPs) between two randomly sampled genomes in a population.
}
\usage{
meansnps(strain.log, freq.log, libr, nuc, key)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{strain.log}{
  Vector of unique genotype IDs in population
}
  \item{freq.log}{
  Vector of counts of genotypes in population, corresponding to \code{strain.log}.
}
  \item{libr}{
  Library object from simulation functions. A list in which each entry represents a unique genotype, and is a vector of mutated nucleotide positions relative to the reference seqeunce.
}
  \item{nuc}{
  Nucleotide database from simulation functions. A list (corresponding to \code{libr}) in which each entry represents a unique genotype, and is a vector of mutated nucleotides relative to the reference sequence.
}
  \item{key}{
  Vector of genotype IDs corresponding to the \code{libr} and \code{nuc} objects.
}
}

\examples{
data(withinhost)
meansnps(withinhost$obs.strain, rep(1, 10), withinhost$libr, withinhost$nuc, 
         withinhost$librstrains)

}

