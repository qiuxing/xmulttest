\name{superdelta}
\alias{superdelta}
\title{Super-delta sequence method for testing gene significance.}

\description{ The Super-delta sequence method is a combination of a)
  gene normalization procedure based on pairwise normalization; b) and a
  gene selection procedure based on a suitable test statistic.  It must
  be applied to log-transformed, but \emph{not normalized} gene
  expression level data.  Only Student's two sample t-test is
  implemented in the current version.  }

\usage{ superdelta(X, classlabel, test="t", side="abs",
                  methods="robust", trim=0.2, ...) }

\arguments{
  \item{X}{
    A data frame or matrix, with \eqn{m} rows corresponding to
    variables (hypotheses) and \eqn{n} columns to observations. In the
    case of gene expression data, rows correspond to genes and columns to
    mRNA samples. The data can be read using \code{\link{read.table}}.
  }

  \item{classlabel}{
    A vector of integers corresponding to observation (column)
    class labels. For \eqn{k} classes, the labels must be integers
    between 0 and \eqn{k-1}.  
  }	

  \item{test}{A character string specifying the statistic to be
    used to test the null hypothesis of no association between the
    variables and the class labels.\cr
    If \code{test="t"}, the tests are based on Student's two-sample t-statistics
    (equal variances). \cr
    If \code{test="f"}, the tests are based on F-statistics.\cr
    If \code{test="pairt"}, the tests are based on paired t-statistics. 
  }

  \item{side}{A character string specifying the type of rejection region.\cr
    If \code{side="abs"}, two-tailed tests, the null hypothesis is rejected for large absolute values of the test statistic.\cr
    If \code{side="upper"}, one-tailed tests, the null hypothesis is
    rejected for large values of the test statistic.\cr
    If \code{side="lower"}, one-tailed tests,  the null hypothesis is rejected for small values of the test statistic.
  }

  \item{methods}{A list of estimating methods to be used with the
    super-delta sequence method.  The default value "robust" should be the
    best method for most applications.  Other options include "mean",
    "median", and "mclust" (very time consuming), these advanced options
    are implemented for specialized situations only.
  }

  \item{trim}{ A real number between 0 and 1.  It controls the
    proportion of "extreme" test statistics to be exluded from the test
    statistic estimation step.  This parameter is similiar in idea, but
    not identical to, the same parameter in function \code{mean()} for
    computing the robust sample mean.  
  }

  \item{...}{Other options passed to the advanced method
    "mclust". Please see help on \code{Mclust} for more details. }
}

\details{
  The super-delta sequence method is a combination of a)
  gene normalization procedure based on pairwise normalization; b) and a
  gene selection procedure based on a suitable test statistic.  It must
  be applied to log-transformed, but \emph{not normalized} gene
  expression level data.  Only Student's two sample t-test (with equal
  variance assumption) is implemented in the current version.  Future
  plan includes implementing Welch's t-test (with unequal variance
  adjustment), paired t-test, and ANOVA F-test.  For more technical details
  of the super-delta sequence method please see the reference.

  A few words about the \code{trim} parameter which controls how robust
  the super-delta sequence method is.  In theory, a safe \code{trim} value
  should be no less than the proportion of true differentially expressed
  genes.  The default value 0.2 should work for most practical
  applications for two reasons: a) it is not that common to have more
  than 20\% genes selected as DEG; b) even if more than 20\% of genes are
  DEG, chances are some of the "outliers" caused by pairing with
  the up-regulated genes get cancelled with those caused by the
  down-regulated genes; c) trimming is by nature a conservative
  operation, so some residual bias will not cause noticeable inflation
  of type I error.  Please see the reference for more details.    
}

\value{ A list with components
  \item{teststat}{Vector of test statistics.}
  \item{rawp}{Vector of raw (unadjusted) \eqn{p}-values.}
}

\author{Xing Qiu}
\seealso{\code{\link{mean}}, \code{\link{t.test}}, \code{\link{p.adjust}}}
\examples{
## The Golub data.  On my computer, it takes about 2.5sec to run, YMMV.
## The result is a data.frame with two columns: teststat and rawp.

data(golub)
smallgd<-golub[1:500,]
classlabel<-golub.cl
system.time(res1 <- superdelta(smallgd, classlabel))

## use both robust and the mclust methods together with a lower-side
## test to find up-regulated genes.  The results are organized in a
## list of two data.frames instead. This command takes about 40 seconds
## to run because Mclust is an expensive computation.
\dontrun{library(mclust); res2 <- superdelta(smallgd, classlabel, methods=c("robust", "mclust"), side="lower")}
}

\keyword{htest}% __ONLY ONE__ keyword per line
