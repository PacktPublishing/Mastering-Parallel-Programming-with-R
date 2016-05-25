\name{createBufferFloatVector}
\alias{createBufferFloatVector}
\title{Creates a Vector of Floats for OpenCL to use as buffer}
\description{
   Creates a Vector of Floats for OpenCL to use as buffer.
   This method is for internal use, users are encouraged to use createBuffer instead.
}
\details{
  Calls clCreateBuffer for an array of floats.
}
\references{
\url{http://www.khronos.org/opencl/}
}
\keyword{OpenCL}