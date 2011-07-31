## turn a vector of t-statistics into two sided pvalues
.t2p <- function(x, df, side="abs"){
  if (side=="abs"){
    return(2-2*pt(abs(x),df))
  } else if (side=="lower"){
    return(pt(x, df))
  } else if (side=="upper"){
    return(1-pt(x, df))
  } else {
    stop(paste("side must be one of the following three options: abs (default), lower, or upper."))
  }
}

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
.est.t <- function(tvec, method=c("robust", "mean", "median", "mclust"), trim=0.2, ...){
  method=match.arg(method)
  if (method=="robust"){
    ## 7/30/2011
    ## median; fold; trim; unfold; mean
    med1 <- median(tvec)
    tvec2 <- abs(tvec-med1)
    t.est <- mean(tvec[tvec2<=quantile(tvec2, 1-trim)])
  } else if (method=="median"){
    t.est <- median(tvec)
  } else if (method=="mean"){
    t.est <- mean(tvec)
  } else if (method=="mclust"){
    ## mclust is not that robust, so use the following wrapper as a backup.
    .substitute <- function(x) {
      warning("Mclust does not converge, use the robust method instead.")
      tvec2 <- abs(tvec)
      t.est <- mean(tvec[tvec2<=quantile(tvec2, 1-trim)])
      return(t.est)
    }
    t.est <- tryCatch(.tclust(tvec, ...), warning=.substitute, error=.substitute)
  } else {
    stop(paste("Method", method, "is not implemented!"))
  }
  ## 7/30/2011.  Try t-stat with discounted variance.
  return(sqrt(2)*t.est)
}

superdelta <- function(X, classlabel, test="t", side="abs",
                       methods="robust", trim=0.2, ...){
  m <- dim(X)[1]; nn <- dim(X)[2]
  teststat <- matrix(0, nrow=m, ncol=length(methods))
  colnames(teststat) <- methods; rownames(teststat) <- rownames(X)
  for (i in 1:m){
    delta.i <- matrix(rep(as.real(X[i,]),m-1),nrow=m-1,byrow=TRUE) -X[-i,]
    if (test=="t"){
      tvec <- rowttests(as.matrix(delta.i), factor(classlabel), tstatOnly=TRUE)$statistic
    } else {
      stop(paste("Test statistic", test, "is not implemented yet!"))
    }
    for (mm in methods){
      teststat[i,mm] <- .est.t(tvec, method=mm, trim=trim, ...)
    }
  }
  ## two sided or one sided test; different test statistics, etc
  if (test=="t"){
    rawp=.t2p(teststat, df=nn-2, side=side)
  } else {
    stop(paste("Test statistic", test, "is not implemented yet!"))
  }
  ## two representations
  if (ncol(teststat)==1){
    rr <- data.frame(teststat=teststat, rawp=rawp)
    colnames(rr) <- c("teststat", "rawp")
    return(rr)
  } else {
    return(list(teststat=teststat, rawp=rawp))
  }
}
