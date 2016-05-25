\name{prandomForest}
\alias{prandomForest}
\alias{prandomForest.default}
\title{Parallel random forest generation}
\description{
The machine learning function prandomForest() is an ensemble tree classifier 
that constructs a forest of classification trees from bootstrap samples of a 
dataset in parallel. The random forest algorithm can be used to classify both
categorical and continuous variables. This function provides a parallel 
equivalent to the serial randomForest() function from the randomForest 
package. 
Note that the randomForest library must be loaded before calling the prandomForest function. 
library("randomForest")

N.B. Please see the SPRINT User Guide for how to run the code in parallel 
using the mpiexec command.
} 
\usage{
prandomForest(x, ...)
\method{prandomForest}{default}(x, y=NULL, xtest=NULL, ytest=NULL, ntree=500, 
                      mtry = if (!is.null(y) && !is.factor(y))
                                 max(floor(ncol(x)/3), 1) 
                             else floor(sqrt(ncol(x))),
                      replace=TRUE, classwt=NULL, cutoff, strata,
                      sampsize = if (replace) nrow(x) 
                                 else ceiling(.632*nrow(x)),
                      nodesize = if (!is.null(y) && !is.factor(y)) 5 else 1,
                      maxnodes=NULL, importance=FALSE, localImp=FALSE, 
                      nPerm=1, proximity, oob.prox=proximity, norm.votes=TRUE,
                      do.trace=FALSE, 
                      keep.forest = !is.null(y) && is.null(xtest), 
                      corr.bias=FALSE, keep.inbag=FALSE, ...)

}

\arguments{
  \item{x}{array of data}
  \item{...}{optional parameters to be passed to the low level function
               randomForest.default.} 
        \item{y}{vector, if a factor, classification is assumed, otherwise 
                 regression is assumed. If omitted, prandomForest() will run 
                 in unsupervised mode.}
        \item{xtest}{data array of predictors for the test set}
        \item{ytest}{response for the test set}
        \item{ntree}{integer, the number of trees to grow}
        \item{mtry}{integer, the number of variables randomly sampled as 
                    candidates at each split. The default value is sqrt(p) for
                    classification and p/3 for regression, where p is the 
                    number of variables in the data matrix x.}
        \item{replace}{boolean, whether the sampling of cases is done with or
                       without replacement. The default value is TRUE.}
        \item{classwt}{vector if priors of the classes. The default value is 
                       NULL.}
        \item{cutoff}{vector of k elements where k is the number of classes. 
                      The winning class for an observation is the one with the
                      maximum ratio of proportion of votes to cutoff. The 
                      default value is 1/k.}
        \item{strata}{variable used for stratified sampling}
        \item{sampsize}{size of sample to draw. For classification, if 
                        sampsize is a vector of the length of the number of 
                        strata, then sampling is stratified by strata, and the
                        elements of sampsize indicate the numbers to be drawn
                        from the strata.}
        \item{nodesize}{integer, the minimum size of the terminal nodes. The 
                        default value is 1 for classification and 5 for 
                        regression.}
        \item{maxnodes}{integer, maximum number of terminal nodes allowed for
                        the trees. The default value is NULL.}
        \item{importance}{boolean, whether the importance of predictors is 
                          assessed. The default value is FALSE.}
        \item{localImp}{boolean, whether casewise importance measure is to be 
                        computed. The default value is FALSE.}
        \item{nPerm}{integer, the number of times the out-of-bag data are 
                     permuted per tree for assessing variable importance. The 
                     default value is one. Regression only.}
        \item{proximity}{boolean, whether the proximity measure among the rows
                         is to be calculated.}
        \item{oob.prox}{boolean, whether the proximity is to be calculated for 
                        out-of-bag data. The default value is set to be the 
                        same as the value of the proximity parameter.}
        \item{norm.votes}{boolean, whether the final result of votes are 
                          expressed as fractions or whether the raw vote 
                          counts are returned. The default value is TRUE. 
                          Classification only.}
        \item{do.trace}{boolean, whether a verbose output is produced. The 
                        default value is FALSE. If set to an integer i then 
                        the output is printed for every i trees.}
        \item{keep.forest}{boolean, whether the forest is returned in the 
                           output object. The default value is FALSE.}   
        \item{corr.bias}{boolean, whether to perform a bias correction. The 
                        default value is FALSE. Regression only.}
        \item{keep.inbag}{boolean, whether the matrix which keeps track of 
                          which samples are in-bag in which trees should be 
                          returned. The default value is FALSE.}
}

\seealso{
    \code{\link[randomForest]{randomForest}}
    \code{\link{SPRINT}}
}
\author{
University of Edinburgh SPRINT Team
    \email{sprint@ed.ac.uk}
    \url{www.r-sprint.org}
    }
\keyword{utilities}
\keyword{interface}