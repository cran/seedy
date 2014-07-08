librtoDNA <-
function(sampleID, libr, nuc, ref.strain, key, sampletime=NULL, strings=FALSE, filename=NULL) {
  ngenomes <- length(sampleID)
  DNAoutput <- NULL
  if (!is.null(filename)) {
    strings <- TRUE
  }
  ntides <- c("C", "A", "G", "T")
  if (!is.null(filename)) {
    write(paste("#nexus\nbegin data;\ndimensions ntax=", length(sampleID), 
        " nchar=", length(ref.strain), ";\nformat datatype=dna symbols=\"CAGT\" missing=? gap=-;\nmatrix\n", sep=""),
        file=filename)
    strings <- TRUE
  }
  for (i in 1:ngenomes) {
    K <- ref.strain
    K[libr[[which(key==sampleID[i])]] ] <- nuc[[which(key==sampleID[i])]]
    for (j in 1:length(ntides)) {
      K[which(K==j)] <- ntides[j]
    }
    if (strings) {
      DNAoutput <- c(DNAoutput, paste(K, collapse=""))
    } else {
      DNAoutput <- rbind(DNAoutput, K)
    }
  }
  if (!is.null(filename)) {
    if (is.null(sampletime)) {
      genomenames <- paste("G", 1:length(sampleID), "_", sampleID, sep="")
    } else {
      genomenames <- paste("G", 1:length(sampleID), "_", sampletime, sep="")
    }
    DNAoutput <- cbind(genomenames, DNAoutput)
    write(t(DNAoutput), file=filename, sep="\t", append=TRUE, ncolumns=2)
    write(";\nend;", file=filename, append=TRUE)
  } else {
    return(DNAoutput)
  }
}
