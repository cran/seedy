simulatebacteria <-
function(eq.size, m.rate, runtime, sample.times, n.samples=1, genomelength=100000, 
         bottle.times=0, bottle.size=1, fullrecords=FALSE, deepseq=FALSE, 
         feedback=1000, init.freq=1, libr=NULL, nuc=NULL, ref.strain=NULL) {
  
  time <- 1 # in bacterial generations
  
  if (is.null(ref.strain)) {
    ref.strain <- sample(1:4, genomelength, replace=T) # reference strain
  }
  if (fullrecords) {
    sample.times <- 1:runtime
    deepseq <- TRUE
  }
  
  totcurstrains <- 1:length(init.freq) # strains present or recorded
  sampledstrains <- NULL
  uniquestrains <- length(init.freq) # Number of unique strain types
  
  # Initialize logs
  if (is.null(libr)) {
    libr <- list()
    mut.nuc <- list()
    
    for (i in 1:uniquestrains) {
      if (i==1) {
        libr[[i]] <- sample(genomelength,1) # Pick random location for mutation
        mut.nuc[[i]] <- sample((1:4)[-ref.strain[libr[[i]]]], 1) # mutation type
      } else {
        libr[[i]] <- c(libr[[i-1]], sample(genomelength, 1)) # Each additional strain is one SNP from previous
        if (libr[[i]][i] %in% libr[[i-1]]) { # If choosing previously chosen locus
          src <- which(libr[[i-1]]==libr[[i]][i])
          mut.nuc[[i]] <- c(mut.nuc[[i-1]], sample((1:4)[-mut.nuc[[i-1]][src]], 1)) # Choose different mutation
        }
        mut.nuc[[i]] <- c(mut.nuc[[i-1]], sample((1:4)[-ref.strain[libr[[i]][i]]], 1)) # don't choose same nucleotide as ref.strain
      }
    }
  } else {
    mut.nuc <- nuc
  }
  if (deepseq || fullrecords) {
    obs.freq <- list()
    obs.strain <- list()
  } else {
    obs.strain <- NULL
    obs.time <- NULL
  }
  freq.log <- init.freq
  strain.log <- totcurstrains

  types <- 1 # Cumulative number of strain types
  while (time <= runtime) {
    if (time%%feedback==0) {
      cat("t=", time, ": no. genotypes=", length(strain.log), ", diversity=", 
          round(meansnps(strain.log, freq.log, libr, nuc, totcurstrains),2), "\n", sep="")
    }
    if (time %in% bottle.times) { # bottleneck?
      cat("t=", time, " Bottleneck size ", bottle.size, "\n", sep="")
      if (length(strain.log)==1) { # if source has clonal infection
        inoc.samp <- rep(strain.log, bottle.size)
      } else {
        inoc.samp <- sample(strain.log, bottle.size, 
                            prob=freq.log, replace=T) # take random sample
      }
      strain.log <- unique(inoc.samp) # distinct types in new infection
      if (length(strain.log)>1) {
        f <- numeric(length(strain.log))
        k <- 1
        for (i in strain.log) {
          f[k] <- sum(inoc.samp==i)
          k <- k+1
        }
        freq.log <- f # frequency of types
      } else {
        freq.log <- bottle.size
      }
    }
    # mutate existing strains for each individual
    pop.size <- sum(freq.log)
    death.prob <- 0.5 + 0.5*(pop.size-eq.size)/eq.size
    if (length(freq.log)==0 || sum(is.na(freq.log))>0 || death.prob>1 || death.prob<0) {
      cat("deathprob=", death.prob, "\npop.size=", pop.size, "\nequi.pop=", eq.size, "generation:", time, "\nFreq.log:\n")
      if (length(freq.log)>0) {
        for (k in 1:length(freq.log)) {
          cat(freq.log[[i]][k], "\n")
        }
      }
    }
    freq.log <- 2*rbinom(length(freq.log), freq.log, 1-death.prob)
    if (length(freq.log)==0 || sum(freq.log)==0) {
      freq.log <- 1
      strain.log <- strain.log[1]
    }
    if (0 %in% freq.log) {
      zeros <- which(freq.log==0)
      freq.log <- freq.log[-zeros]
      strain.log <- strain.log[-zeros]
    }
    if (length(freq.log)!=length(strain.log)) {
      stop("prob")
    }
    n.mutations <- rbinom(1, sum(freq.log), m.rate)
    if (n.mutations > 0) {
      for (mt in 1:n.mutations) {
        types <- types+1
        if (length(strain.log)==1) {
          mutate.grp <- strain.log
        } else {
          mutate.grp <- sample(strain.log, 1, prob=freq.log)
        }
        if (mutate.grp %in% totcurstrains) {
          mut.loc <- sample(genomelength, 1)
          mut.nuc[[length(totcurstrains)+1]] <- 
            mut.nuc[[which(totcurstrains==mutate.grp)]]
          if (mut.loc %in% libr[[which(totcurstrains==mutate.grp)]]) { # if mutation at existing location
            kn <- which(libr[[which(totcurstrains==mutate.grp)]]==mut.loc)
            mut.nuc[[length(totcurstrains)+1]][kn] <- sample((1:4)[-mut.nuc[[which(totcurstrains==mutate.grp)]][kn]], 1)
            libr[[length(totcurstrains)+1]] <- libr[[which(totcurstrains==mutate.grp)]]
          } else {
            mut.nuc[[length(totcurstrains)+1]] <- 
              c(mut.nuc[[length(totcurstrains)+1]], 
                sample((1:4)[-ref.strain[mut.loc]], 1))
            libr[[length(totcurstrains)+1]] <- 
              c(libr[[which(totcurstrains==mutate.grp)]], mut.loc)
          }
          strain.log <- c(strain.log, types)
          freq.log <- c(freq.log, 1)
          totcurstrains <- c(totcurstrains, types)
        }
      }
    }
    # take samples, make observations
    if (time%in%sample.times) {
      if (deepseq) {
        n <- length(obs.freq)+1
        obs.freq[[n]] <- freq.log
        obs.strain[[n]] <- strain.log
        sampledstrains <- unique(c(sampledstrains, unique(strain.log)))
      } else {
        if (length(strain.log)==1) {
          obs.strain <- rep(strain.log, n.samples)
        } else {
          obs.strain <- c(obs.strain, sample(strain.log, n.samples, prob=freq.log, replace=TRUE))
          obs.time <- c(obs.time, rep(time, n.samples))
        }
        sampledstrains <- unique(obs.strain)
      }
    } 
    if (time==runtime) {
      strain.log = 0
      freq.log = 0
    }
    if (!fullrecords) {
      deleters <- NULL
#       # clean up libr etc.
#       deleters <- NULL
#       uniquestrains <- 0
#       for (j in 1:length(totcurstrains)) {
#         tottype <- 0
#         if (totcurstrains[j]%in%strain.log) { # if strain is extant
#           tottype <- tottype+1
#         }
#         if (sum(!strain.log%in%totcurstrains)>0 && time!=runtime) {
#           stop("Deleted sequence for observed sample")
#         }
#         recorded <- FALSE
#         if (length(obs.strain)!=0) {
#           if (totcurstrains[j] %in% sampledstrains) {
#             recorded <- TRUE
#           }
# #           if (deepseq) {
# #             for (w in 1:length(obs.strain)) {
# #               if (totcurstrains[j] %in% obs.strain[[w]]) {
# #                 recorded <- TRUE
# #               }
# #             }
# #           } else {
# #             if (totcurstrains[j] %in% obs.strain) {
# #               recorded <- TRUE
# #             }
# #           }
#         }
      if (sum(!totcurstrains%in%strain.log)>0) {
        notpresentpos <- which(!totcurstrains%in%strain.log)
        notpresent <- totcurstrains[notpresentpos]
        if (sum(notpresent%in%sampledstrains)>0) {
          deleters <- notpresentpos[which(!notpresent%in%sampledstrains)]
        } else {
          deleters <- notpresentpos
        }
      }

#         if (tottype>0) { # don't delete if still around
#           uniquestrains <- uniquestrains+1
#         } else if (tottype==0 && !recorded) { # if not around AND not logged
#           deleters <- c(deleters, j) # delete
#         }
#       }
      if (length(deleters)>0) {
        for (i in sort(deleters, decreasing=T)) {
          libr[[i]] <- NULL
          mut.nuc[[i]] <- NULL
        }
        deletegroup <- totcurstrains[deleters]
        totcurstrains <- totcurstrains[-deleters]
      }
    }
    time <- time+1
  }
  if (deepseq) {
    return(invisible(list(libr=libr, nuc=mut.nuc, librstrains=totcurstrains, obs.freq=obs.freq, obs.strain=obs.strain,
                          ref.strain=ref.strain)))
  } else {
    return(invisible(list(libr=libr, nuc=mut.nuc, librstrains=totcurstrains, obs.strain=obs.strain, obs.time=obs.time, 
                          ref.strain=ref.strain)))    
  }
}


