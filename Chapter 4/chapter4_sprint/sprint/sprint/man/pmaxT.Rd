\name{pmaxT}
\alias{pmaxT}
\title{Ajusted p-values for step down multiple testing}
\description{
Function which implements a parallel version of the multtest "mt.maxT" 
function. It computes the adjusted p-values for step-down multiple testing 
procedures.

N.B. Please see the SPRINT User Guide for how to run the code in parallel 
using the mpiexec command.
}
\usage{
pmaxT(X, classlabel, test = "t", side = "abs", fixed.seed.sampling = "y", 
      B = 10000, na = .naNUM, nonpara = "n") 
}
\arguments{
    \item{X}{double precision 2D array of data}
    \item{classlabel}{class column labels of the input data array}
    \item{test}{one of the following statistical test: t. t.equalvar, 
                Wilcoxon, F, Pair-T, Block-F}
    \item{side}{Type of rejection region, either abs, upper or lower}
    \item{fixed.seed.sampling}{whether the permutations are calculated on the 
                               fly or save to memory}
    \item{B}{number or permutations}
    \item{na}{missing value tag}
    \item{nonpara}{whether non-parametric test statistics or not}
}
\seealso{
    \code{\link[multtest]{mt.maxT}}
    \code{\link{SPRINT}}
}
\author{
University of Edinburgh SPRINT Team
    \email{sprint@ed.ac.uk}
    \url{www.r-sprint.org}
    }
\keyword{utilities}
\keyword{interface}