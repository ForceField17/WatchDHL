# This function sorts the matrix for better visualization of mutual exclusivity across genes
memoSort <- function(M) {
  geneOrder <- sort(rowSums(M), decreasing=TRUE, index.return=TRUE)$ix;
  scoreCol <- function(x) {
    score <- 0;
    for(i in 1:length(x)) {
      if(x[i]) {
        score <- score + 2^(length(x)-i);
      }
    }
    return(score);
  }
  scores <- apply(M[geneOrder, ], 2, scoreCol);
  sampleOrder <- sort(scores, decreasing=TRUE, index.return=TRUE)$ix;
  return(M[geneOrder, sampleOrder]);
}

# This is the plotting function
oncoPrint <- function(M, sort=TRUE) {
  if(sort) {
    alts <- memoSort(M);		
  } else {
    alts <- M;
  }
  
  ngenes <- nrow(alts);
  nsamples <- ncol(alts);
  coverage <- sum(colSums(alts) > 0);
  
  ### OncoPrint
  numOfOncos <- ngenes*nsamples;
  oncoCords <- matrix( rep(0, numOfOncos * 5), nrow=numOfOncos );
  colnames(oncoCords) <- c("xleft", "ybottom", "xright", "ytop", "altered");
  
  xpadding <- .01;
  ypadding <- .01;
  cnt <- 1;
  for(i in 1:ngenes) {
    for(j in 1:nsamples) {
      xleft <- j-1 + xpadding;
      ybottom <- ((ngenes-i+1) -1) + ypadding;
      xright <- j - xpadding;
      ytop <- (ngenes-i+1) -ypadding;
      altered <- alts[i, j];
      
      oncoCords[cnt, ] <- c(xleft, ybottom, xright, ytop, altered);
      cnt <- cnt+1;
    }
  }
  
  colors <- rep("lightgray", cnt);
  colors[ which(oncoCords[, "altered"] == 1) ] <- "#a1d76a";
  colors[ which(oncoCords[, "altered"] == 2) ] <- "#e9a3c9";
  plot(c(0, nsamples), c(0, ngenes), type="n", main=sprintf("Gene set altered in %.2f%%: %d of %d cases", coverage/nsamples*100, coverage, nsamples), xlab="Samples", ylab="", yaxt="n");
  rect(oncoCords[, "xleft"], oncoCords[, "ybottom"],oncoCords[, "xright"], oncoCords[, "ytop"], col=colors, border="white");
  axis(2, at=(ngenes:1)-.5, labels=rownames(alts), las=2);
}

# M is the alteration matrix with T => Altered, F => Not-altered
# Rows are genes, columns are samples
M <- matrix(rep(FALSE, 100*3), ncol=100, nrow=3);
rownames(M) <- c("Gene1", "Gene2", "Gene2");
colnames(M) <- paste("Sample", 1:100, sep="");
M[1, 1:25] <- TRUE;
M[2, 26:50] <- TRUE;
M[3, 51:75] <- TRUE;

# We don't care about the order, so samples can be ordered in any way before plotting
M <- M[, sample(1:100, 100)];
# End of custom alteration matrix

# Create the plot
oncoPrint(M);

# You can also disable the sorting
oncoPrint(M, sort=FALSE)#