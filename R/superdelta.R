## turn a vector of t-statistics into two sided pvalues
.t2p <- function(x,df) 2-2*pt(abs(x),df)

## this function takes tvec, which is a vector of t-statistics
## associated with the ith gene normalized by other genes.  It
## estimates one best t-stat according to different methods.
.est.t <- function(tvec, freq.thresh=0.8, method=c("freq", "mclust", "prefilter"), ...){
  
}

my.freq <- function(tipi,freq.thresh){
  ind <- m * freq.thresh
  tipi.sorted <- tipi[with(tipi, order(p.value)),]
  c(tipi.sorted[ind,"statistic"],tipi.sorted[ind,"p.value"])
}
my.mclust <- function(tipi, ...){
  ## Use mclust to split the t-stats computed from all pairs into 3
  ## clusters: paired with down regulated, paired with NDEG, paired
  ## with up-regulated genes.
  t.i <- tipi[,"statistic"]
  mm.i <- summary(mclustBIC(t.i, G=3, modelName="V"),t.i)
  NC1 <- sum(mm.i$classification==1); NC2 <- sum(mm.i$classification==2); NC3 <- sum(mm.i$classification==3)
  largest.cluster <- which.max(c(NC1, NC2, NC3))
  ti.est <- as.real(mm.i$parameters$mean[largest.cluster])
  c(ti.est, .t2p(ti.est,nn-2))
}


superdelta <- function(X, classlabel, methods=c("freq", "mclust"),freq.thresh=0.8, ...){
  m <- dim(X)[1]; nn <- dim(X)[2]
  results <- matrix(-1,nrow=m, ncol=2*length(methods))
  rownames(results) <- rownames(X)
  colnames(results) <- as.vector(outer(c("statistic","p.value"),methods,paste,sep="."))
  for (i in 1:m){
    delta.i <- matrix(rep(as.real(X[i,]),m-1),nrow=m-1,byrow=TRUE) -X[-i,]
    tipi <- rowttests(as.matrix(delta.i), factor(classlabel))
    for (mm in methods){
      if (mm=="freq"){
        results.freq <- my.freq(tipi,freq.thresh)
        results[i,paste("statistic",mm,sep=".")] <- results.freq[1]
        results[i,paste("p.value",mm,sep=".")] <- results.freq[2]
      } else if (mm=="mclust"){
        ## Since mclust isn't 100% safe to use, we need to define a
        ## substitude function IN CASE IT FAILS. Here we choose the freq
        ## algorithm.
        my.substitute <- function(x) {print(i); print(x); my.freq(tipi,freq.thresh)}
        results.mclust <- tryCatch(my.mclust(tipi,...), error=my.substitute)
        results[i,paste("statistic",mm,sep=".")] <- results.mclust[1]
        results[i,paste("p.value",mm,sep=".")] <- results.mclust[2]
      } else {stop(paste("Method", mm, "is not implemented!")) }
    }
  }
  results
}
