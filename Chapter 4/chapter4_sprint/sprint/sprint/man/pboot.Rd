\name{pboot}
\alias{pboot}
\title{Parallel Bootstrapping}
\description{
pboot() generates R bootstrap replicates of a statistic applied to data. It 
implements a parallel version of the bootstrapping method boot() from the boot
R package.
}
\arguments{
    \item{data}{array of data, if a 2D array then each row is considered as 
                one multivariate observation}
    \item{statistic}{function, when sim is set to parametric, 
                     the first argument to statistic must be the data. For 
                     each replicate a simulated dataset returned by ran.gen 
                     will be passed. In all other cases, statistic must take 
                     at least two arguments. The first argument passed will 
                     always be the original data. The second will be a vector 
                     of indices, frequencies or weights which define the 
                     bootstrap sample.}
    \item{R}{number of bootstrap replicates}
    \item{sim}{string, indicates the type of simulation. The default value is 
               "ordinary". Other possible values are parametric, balanced,
               permutation, and antithetic. Importance resampling is 
               specified by including importance weights; the type of 
               importance resampling must still be specified but may only be 
               ordinary or balanced in this case.}
    \item{stype}{string, indicates what the second argument of statistic 
                 represents. The default value is i for indices. Other 
                 possible values are f for frequencies and w for weights.
                 It is not used when sim is set to parametric.}
    \item{strata}{vector of integer, specifies the strata for multi-sample 
                  problems. This may be specified for any simulation, but is 
                  ignored when sim is set to parametric. When strata is 
                  supplied for a nonparametric bootstrap, the simulations are 
                  done within the specified strata.}
    \item{L}{vector of influence values evaluated at the observations. This is
             used only when sim is set to antithetic. If not supplied, they
             are calculated through a call to empinf.  This will use the 
             infinitesimal jackknife provided that stype is set to w 
             otherwise the usual jackknife is used.}
    \item{m}{the number of predictions which are to be made at each bootstrap
             replicate. This is most useful for (generalized) linear models. 
             This can only be used when sim is ordinary. m will usually be 
             a single integer but, if there are strata, it may be a vector 
             with length equal to the number of strata, specifying how many 
             of the errors for prediction should come from each strata. The
             actual predictions should be returned as the final part of the 
             output of statistic, which should also take an argument giving 
             the vector of indices of the errors to be used for the 
             predictions.}
    \item{weights}{array of importance weights. If a vector then it should 
                   have as many elements as there are observations in the 
                   input data. When simulation from more than one set of 
                   weights is required, weights should be a matrix where each
                   row of the matrix is one set of importance weights. If 
                   weights is a matrix then the number of bootstrap replicates
                   R must be a vector of length nrow(weights). This parameter
                   is ignored if sim is not set to ordinary or balanced.}
    \item{ran.gen}{function, used only when sim is set to parametric. It 
                   describes how random values are to be generated. It should
                   be a function of two arguments. The first argument should 
                   be the observed data and the second argument consists of 
                   any other information needed (e.g. parameter estimates). 
                   The second argument may be a list, allowing any number of 
                   items to be passed to ran.gen. The returned value should be
                   a simulated data set of the same form as the observed data
                   which will be passed to statistic to get a bootstrap 
                   replicate. It is important that the returned value be of 
                   the same shape and type as the original dataset. If ran.gen
                   is not specified, the default is a function which returns 
                   the original input data in which case all simulation should
                   be included as part of statistic. Setting sim to 
                   parametric and using a suitable ran.gen allows the user 
                   to implement any types of nonparametric resampling which 
                   are not supported directly.}
    \item{mle}{secong argument to ran.gen, typically these will be maximum 
               likelihood estimates of the parameters. For efficiency mle is 
               often a list containing all of the objects needed by ran.gen 
               which can be calculated using the original data set only.}
    \item{simple}{boolean, can only be set to TRUE if sim is set to 
                  ordinary, stype is set to I and n is set to 0. Otherwise
                  it is ignored and generates a warning. By default a n by R 
                  index array is created which can be large. If simple is set
                  to TRUE, this is avoided by sampling separately for each 
                  replication, which is slower but uses less memory.}
    \item{...}{other named arguments for statistic which are passed unchanged
               each time.}                  
}
\details{
This version is an early but fully working prototype. However, it is not 
compatible with other SPRINT functions, i.e. you cannot bootstrap other 
parallel functions from the SPRINT library. It is therefore recommended to use
it only as a standalone function. 

N.B. Please see the SPRINT User Guide for how to run the code in parallel 
using the mpiexec command.
}
\seealso{
    \code{\link{boot}}
    \code{\link{SPRINT}}
}
\author{
University of Edinburgh SPRINT Team
    \email{sprint@ed.ac.uk}
    \url{www.r-sprint.org}
    }
\keyword{utilities}
\keyword{interface}