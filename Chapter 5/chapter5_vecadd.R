##
## Copyright 2015 Simon Chapple
##
## Packt: "Mastering Parallelism with R"
## Chapter 5 - ROpenCL Vector Addition
##

# You need to set the value of /path below to the location of your chapter6 source code
install.packages("/path/ROpenCL",repos=NULL,type="source")
require(ROpenCL)

# First look-up the GPU and create the OpenCL Context
platformIDs <- getPlatformIDs()
deviceID <- getDeviceIDs(platformIDs[[1]])[[1]]
info <- getDeviceInfo(deviceID) 
context <- createContext(deviceID)

# Initialise the input data in R on the CPU (Host)
# and pre-allocate the output result
aVector <- seq(1.0, 12345678.0, by=1.0)   # Long numeric vector
bVector <- seq(12345678.0, 1.0, by=-1.0)  # Same but in reverse
cVector <- rep(0.0, length(aVector))      # Similar result vector

localWorkSize = 16  # GPU/kernel dependent (explained later)
# globalWorkSize must be integer multiple of localWorkSize
globalWorkSize = ceiling(length(aVector) / localWorkSize) *      
                         localWorkSize

# Allocate the Device’s global memory Buffers: 2x input, 1x output 
aBuffer <- createBuffer(context,"CL_MEM_READ_ONLY",
                        length(aVector),aVector)
bBuffer <- createBuffer(context,"CL_MEM_READ_ONLY",
                        length(aVector),bVector)
cBuffer <- createBufferFloatVector(context,"CL_MEM_WRITE_ONLY",
                        length(aVector))

# Create the OpenCL C Kernel function to add two vectors
kernelSource <- '
__kernel void vectorAdd(__global float *a, __global float *b,        
                        __global float *c, int numDataItems)
{
  int gid = get_global_id(0); // WorkItem index in 1D global range
  if (gid >= numDataItems) return; // Exit fn if beyond data range
  c[gid] = a[gid] + b[gid]; // Perform addition for this WorkItem
}'
vecAddKernel <- buildKernel(context,kernelSource,'vectorAdd',
                            aBuffer,bBuffer,cBuffer,length(aVector))

# Create a device command queue
queue <- createCommandQueue(context,deviceID)
# Prime the two input Buffers
enqueueWriteBuffer(queue,aBuffer,length(aVector),aVector)
enqueueWriteBuffer(queue,bBuffer,length(aVector),bVector)
# Execute the Kernel
enqueueNDRangeKernel(queue,vecAddKernel,
                     globalWorkSize,localWorkSize)
# Retrieve the calculated result
enqueueReadBuffer(queue,cBuffer,length(aVector),cVector)

##
## Packt: "Mastering Parallelism with R"
## Chapter 6 - ROpenCL Vector Addition
##
## Copyright 2015 Simon Chapple
##
