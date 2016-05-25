\name{papply}
\alias{papply}
\title{Parallel Apply}
\description{ 
Parallel apply function is used to perform the same operation over all the 
elements of data objects like matrices, or lists. This function 
provides a parallel implementation of both the apply() and lapply() functions
from the core of the R programming language. The papply() function only 
accepts matrices or lists of matrices as input data objects. 
}
\usage{
papply(data, fun, margin = 1, out_filename = NULL)
}
\arguments{
    \item{data}{input data matrix or list of matrices}
    \item{fun}{function to be applied}
    \item{margin}{vector indicating which elements of the matrix the function
                  will be applied to. The default value is 1 and indicates the
                  rows, 2 indicates the columns and the parameter is ignored 
                  if data is a list.}
    \item{out_filename}{string, not used at present.}
}
\details{
The function to be applied can be supplied to papply() either as a function 
name or as a function definition. When only the function name is provided, 
the package implementing the function has to be loaded before the SPRINT 
library is initialised in order to ensure that the name is recognised by all 
the processes involved in the computation.

Note that papply() does not fully implement lapply functionality. It will 
only accept lists of matrices, and not lists made up of other data types. 

N.B. Please see the SPRINT User Guide for how to run the code in parallel 
using the mpiexec command.
}
\seealso{
    \code{\link{apply}}
    \code{\link{lapply}}
    \code{\link[ff]{ff}}
    \code{\link{SPRINT}}
}
\author{
University of Edinburgh SPRINT Team
    \email{sprint@ed.ac.uk}
    \url{www.r-sprint.org}
    }
\keyword{utilities}
\keyword{interface}