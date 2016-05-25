\name{ppam}
\alias{ppam}
\title{Parallel Partitioning Around Medoids}
\description{
Parallel implementation of the Partitioning Around Medoids algorithm, based
on the cluster "pam" serial function. 
}
\usage{
ppam(x, k, medoids = NULL, is_dist = inherits(x, "dist"), 
     cluster.only = FALSE, do.swap = TRUE, trace.lev = 0)
}
\arguments{
    \item{x}{input data, either a 2D array or an ff object}
    \item{k}{positive integer, indicating for the number of clusters}
    \item{medoids}{vector, with the initial 'k' medoids or NULL to let the 
                   algorithm select the initial medoids}
    \item{is_dist}{boolean, whether the input data is a distance or 
                   dissimilarity matrix or a symmetric matrix}
    \item{cluster.only}{boolean, whether only the clustering is computed and 
                        returned}
    \item{do.swap}{boolean, whether the swap phase of the algorithm is 
                   required}
    \item{trace.lev}{positive integer for the level of details returned for
                     diagnostics}
}
\details{
The interface and parameters to parallel function ppam() are similar to the 
serial function pam() but not identical. ppam() requires a distance matrix as
input parameters. Although, ppam() does not include the option to calculate 
the distance matrix, this can easily be done using SPRINT pcor() function with the 'distance' parameter set to TRUE. 

N.B. Please see the SPRINT User Guide for how to run the code in parallel 
using the mpiexec command.
}
\seealso{
    \code{\link[cluster]{pam}}
    \code{\link[ff]{ff}}
    \code{\link{pcor}}
    \code{\link{SPRINT}}
}
\author{
University of Edinburgh SPRINT Team
    \email{sprint@ed.ac.uk}
    \url{www.r-sprint.org}
    }
\keyword{utilities}
\keyword{interface}