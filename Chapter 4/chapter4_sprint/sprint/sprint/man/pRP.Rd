\name{pRP}
\alias{pRP}
\title{Parallel Rank Product}
\description{
Parallel rank product function helps identifying differentially 
regulated genes in replicated microarray experiments.
}
\usage{
pRP(data, cl, num.perm, logged = TRUE, na.rm = FALSE, gene.names = NULL, 
    plot = FALSE, rand = NULL, sum = FALSE)
}
\arguments{
    \item{data}{array, input data}
    \item{cl}{vector, class labels of the samples}
    \item{num.perm}{integer, the number of permutations used in the 
                    calculation of the null density. The default value is 
                    100.}
    \item{logged}{boolean, whether the data is logged or not. The default 
                  value is TRUE.}
    \item{na.rm}{boolean, whether missing values are to be replaced by the 
                 gene-wise mean of the non-missing values and used in 
                 computing rank. The default value is FALSE.}
    \item{gene.names}{the gene name to be assigned to the estimated percentage
                      of false positive predictions. The default value is 
                      NULL.}
    \item{plot}{boolean, whether to plot the estimated percentage of false 
                positive predictions against the rank of each gene. The 
                default value is FALSE.}
    \item{rand}{number, the seed for the random number generator if specified.
                The default value is NULL.}
    \item{sum}{boolean, whether to perform a rank sum analysis. The default 
               value is NULL.}            
}
\details{
The SPRINT task parallel implementation of the rank product method is 
approximately twice as fast in serial as the existing RP() function from the RankProd package and it shows excellent scaling.

N.B. Please see the SPRINT User Guide for how to run the code in parallel 
using the mpiexec command.
}
\seealso{
    \code{\link[RankProd]{RP}}
    \code{\link{SPRINT}}
}
\author{
University of Edinburgh SPRINT Team
    \email{sprint@ed.ac.uk}
    \url{www.r-sprint.org}
    }
\keyword{utilities}
\keyword{interface}