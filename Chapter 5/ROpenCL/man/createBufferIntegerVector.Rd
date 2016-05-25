\name{createBufferIntegerVector}
\alias{createBufferIntegerVector}
\title{Creates a Vector of Integers for OpenCL to use as buffer}
\description{
   Creates a Vector of Integers for OpenCL to use as buffer.
   This method is for internal use, users are encouraged to use createBuffer instead.
}
\details{
  Calls clCreateBuffer for an array of integers.
}
\references{
\url{http://www.khronos.org/opencl/}
}
\keyword{OpenCL}