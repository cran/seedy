plotoutbreak <-
function(epidata, sampledata, col="red", stack=TRUE, arr.len=0.1, blockheight=0.5, hspace=max(epidata[,3])/20,
                          labels=NULL, label.pos="left", block.col="grey", jitter=0, pch=1, ...) {
  IDvec <- epidata[,1]
  inf.times <- epidata[,2]
  rec.times <- epidata[,3]
  inf.source <- epidata[,4]
  
  sampleID <- sampledata[,1]
  sample.times <- sampledata[,2]
  if (length(col)==1) {
    col <- rep(col, nrow(sampledata))
  }
  if (length(block.col)==1) {
    block.col <- rep(block.col, max(IDvec))
  }
  if (is.null(labels)) {
    labels <- IDvec
  }
  tot.inf <- length(rec.times)
  if (stack) {
    row <- c(0,rep(NA,tot.inf-1))
    for (i in 2:tot.inf) {
      possrows <- NULL
      for (j in 1:(i-1)) {
        rowrec <- max(rec.times[which(row==row[j])])
        if (inf.times[i]>rowrec+hspace) {
          possrows <- c(possrows, row[j])
        }
      }
      possrows <- c(possrows, max(row, na.rm=TRUE)+1, min(row, na.rm=TRUE)-1)
      rowdist <- abs(possrows-row[which(IDvec==inf.source[i])])
      K <- which(rowdist==min(rowdist))
      if (length(K)>1) {
        K <- K[-which(K>max(row, na.rm=TRUE))]
        K <- K[-which(K<min(row, na.rm=TRUE))]
        if (length(K)==0) {
          K <- which(rowdist==min(rowdist))[1]
        }
      }
      row[i] <- possrows[K]
    }
  } else {
    row <- 1:tot.inf
  }
  plot(NULL, xlim=c(0,max(rec.times)), ylim=c(min(row)-0.5, max(row)+0.5), yaxt="n", ylab="", ...)
  for (i in 1:tot.inf) {
    rect(inf.times[i], row[i]-blockheight/2, rec.times[i], row[i]+blockheight/2, col=block.col[IDvec[i]])
    if (pch>20) {
      points(sample.times[which(sampleID==IDvec[i])]+rnorm(sum(sampleID==IDvec[i]), 0, jitter*max(rec.times)), 
             rep(row[i], sum(sampleID==IDvec[i]))+rnorm(sum(sampleID==IDvec[i]), 0, jitter*(max(row)-min(row))), 
             bg=col[which(sampleID==IDvec[i])], pch=pch)
    } else {
      points(sample.times[which(sampleID==IDvec[i])]+rnorm(sum(sampleID==IDvec[i]), 0, jitter*max(rec.times)), 
             rep(row[i], sum(sampleID==IDvec[i]))+rnorm(sum(sampleID==IDvec[i]), 0, jitter*(max(row)-min(row))), 
             col=col[which(sampleID==IDvec[i])], pch=pch)    
    }
    if (label.pos=="centre") {
      text(rec.times[i]-inf.times[i], row[i],labels[i], cex=0.75)      
    } else if (label.pos=="right") {
      text(rec.times[i], row[i],labels[i], cex=0.75, pos=4, offset=0.2)      
    } else if (label.pos=="left") {
      text(inf.times[i], row[i],labels[i], cex=0.75, pos=2, offset=0.2)      
    }
    if (i>1) {
      if (row[which(IDvec==inf.source[i])]>row[i]) {
        arrows(inf.times[i],row[which(IDvec==inf.source[i])]-blockheight/2, inf.times[i], row[i]+blockheight/2, 
               col="black", angle=20, length=arr.len)
      } else {
        arrows(inf.times[i],row[which(IDvec==inf.source[i])]+blockheight/2, inf.times[i], row[i]-blockheight/2, 
               col="black", angle=20, length=arr.len)
      }
    }
  }
  
}
