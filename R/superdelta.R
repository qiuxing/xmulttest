## turn a vector of t-statistics into two sided pvalues
.t2p <- function(x,df) 2-2*pt(abs(x),df)

## cluster tvec by mclust
.tclust <- function(tvec, ...) {
  mm.i <- summary(mclustBIC(tvec, G=3, modelName="V", ...), tvec)
  NC1 <- sum(mm.i$classification==1); NC2 <- sum(mm.i$classification==2); NC3 <- sum(mm.i$classification==3)
  largest.cluster <- which.max(c(NC1, NC2, NC3))
  t.est <- as.real(mm.i$parameters$mean[largest.cluster])
  return(t.est)
}

## this function takes tvec, which is a vector of t-statistics
## associated with the ith gene normalized by other genes.  It
## estimates one best t-stat according to different methods.
.est.t <- function(tvec, method=c("robust", "median", "mclust"), filter=0.2, ...){
  method=match.arg(method)
  if (method=="robust"){
    tvec2 <- abs(tvec)
    t.est <- mean(tvec[tvec2>quantile(tvec2, filter)])
  } else if (method=="median"){
    t.est <- median(tvec)
  } else if (method=="mclust"){
    ## mclust is not that robust, so use the following wrapper as a backup.
    .substitute <- function(x) {
      warning("Mclust does not converge, use the robust method instead.")
      tvec2 <- abs(tvec)
      t.est <- mean(tvec[tvec2>quantile(tvec2, filter)])
      return(t.est)
    }
    t.est <- tryCatch(.tclust(tvec, ...), warning=.substitute, error=.substitute)
  } else {
    stop(paste("Method", method, "is not implemented!"))
  }
  return(t.est)
}

superdelta <- function(X, classlabel, methods=c("robust", "median"), ...){
  m <- dim(X)[1]; nn <- dim(X)[2]
  teststat <- matrix(0, nrow=m, ncol=length(methods))
  colnames(teststat) <- methods; rownames(teststat) <- rownames(X)
  for (i in 1:m){
    delta.i <- matrix(rep(as.real(X[i,]),m-1),nrow=m-1,byrow=TRUE) -X[-i,]
    tvec <- rowttests(as.matrix(delta.i), factor(classlabel), tstatOnly=TRUE)$statistic
    for (mm in methods){
      teststat[i,mm] <- .est.t(tvec, method=mm, ...)
    }
  }
  list(teststat=teststat, rawp=.t2p(teststat, df=nn-2))
}
