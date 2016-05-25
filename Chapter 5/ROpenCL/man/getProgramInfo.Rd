\name{getProgramInfo}
\alias{getProgramInfo}
\title{Retreives information about the given OpenCL program.}
\description{
   Retreives information about the given OpenCL program.
}
\details{
  Calls clGetProgramInfo and returns the information about the OpenCL program as a character string.
  Currently, only CL_PROGRAM_SOURCE is supported.
}
\references{
\url{http://www.khronos.org/opencl/}
}
\keyword{OpenCL}