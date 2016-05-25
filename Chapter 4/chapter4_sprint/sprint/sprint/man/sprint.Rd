\name{About SPRINT}
\alias{SPRINT}
\docType{package}
\title{Overview of SPRINT}
\description{
SPRINT (Simple Parallel R INTerface) is a parallel framework for R. It 
provides a High Performance Computing (HPC) harness which allow R scripts to 
run on HPC clusters. SPRINT contains a library of selected R functions that 
have been parallelized. Functions are named after the original R function with
the added prefix 'p', i.e. the parallel version of cor() in SPRINT is called 
pcor(). These parallelized functions are written in C and MPI. Call to these 
functions are included directly in standard R scripts. 

The following functions are implemented in SPRINT 1.5.0:
 - papply
 - pboot
 - pcor
 - pmaxT
 - ppam
 - prandomForest
 - pRP 
 - pstringdistmatrix
 - pterminate
 - ptest

See the User Guide and Release Notes in the sprint folder or the SPRINT web 
page at http://www.r-sprint.org for more information.
}
\details{
To make use of SPRINT it is necessary to include the library first. Then 
include calls to the SPRINT functions you want to use. It is also necessary 
to exit SPRINT using the pterminate function which shutdown MPI as well as
SPRINT.
}
\examples{
library("sprint")
ptest()
pterminate()
quit()
}
\author{
University of Edinburgh SPRINT Team
    \email{sprint@ed.ac.uk}
    \url{www.r-sprint.org}
    }
\keyword{utilities}
\keyword{interface}