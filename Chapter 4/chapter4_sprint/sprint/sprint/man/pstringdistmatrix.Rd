\name{pstringdistmatrix}
\alias{pstringdistmatrix}
\title{Parallel compute hamming distance between strings}
\description{Calculates the hamming distance matrix between strings.
}
\usage{ 
  pstringdistmatrix(a, b, method = "h", filename = NULL, weight = NULL,
                 maxDist = 0, ncores = NULL)
}
\arguments{
  \item{a}{a: R object (target); will be converted by 'as.character'.}
  \item{b}{b: R object (source); will be converted by 'as.character'. Must be the same as argument a in this version of the software.}
  \item{method}{Method for distance calculation - only option 'h' for hamming distance is supported.}
  \item{filename}{Results will be stored here as binary data}
  \item{weight}{Not used in the hamming distance measure.}
  \item{maxDist}{Not used in the hamming distance measure.}
  \item{ncores}{Not used by SPRINT, please see the SPRINT user guide.}
}
\details{
  Calculates the hamming distance between each pair of strings. Returns an ff result matrix.

N.B. Please see the SPRINT User Guide for how to run the code in parallel 
using the mpiexec command.
}
\seealso{
    \code{\link[stringdist]{stringdistmatrix}}
    \code{\link{SPRINT}}
}
\author{
  University of Edinburgh SPRINT Team
  \email{sprint@ed.ac.uk}
  \url{www.r-sprint.org}
}
\keyword{utilities}
\keyword{interface}

