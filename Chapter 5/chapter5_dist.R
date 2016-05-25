##
## Copyright 2015 Simon Chapple
##
## Packt: "Mastering Parallelism with R"
## Chapter 5 - ROpenCL Distance Function
##

# You need to set the value of /path below to the location of your chapter6 source code
install.packages("/path/ROpenCL",repos=NULL,type="source")
filepath <- "/path/chapter6_IMD_data.csv"
require(ROpenCL)

teval <- function(...) {
  gc();
  start <- Sys.time()
  pstart <- proc.time()
  result <- eval(...)
  pfinish <- proc.time()
  finish <- Sys.time()
  wallclock = list(Start=format(start,"%H:%M:%OS3"),Finish=format(finish,"%H:%M:%OS3"))
  return ( list(PID=Sys.getpid(), WallClock=wallclock, Duration=pfinish-pstart, Result=result) )
}

distStart <- function(obsIndex,numObs) {
  return( as.integer(numObs*(obsIndex-1) - obsIndex*(obsIndex-1)/2 ) )
}

# Identify the GPU device
platformIDs <- getPlatformIDs()
deviceID <- getDeviceIDs(platformIDs[[1]])[[2]] 
context <- createContext(deviceID)
queue <- createCommandQueue(context,deviceID)

data <- read.table(file = filepath, header=TRUE, sep=",", row.names=1)
head(data)

# Create a 1D vector of observations X variables from the data table
dvec <- as.vector(t(data))

localWorkSize <- 16
globalWorkSize <- 32768 # closest factor of localWorkSize
maxWorkSize <- 32482 # total number of observations
blockWorkSize <- 2048 # device memory limits max 4096 obs to process in one go
distSizeMax <- distStart(maxWorkSize,maxWorkSize)
distSizeBlock <- distStart(blockWorkSize+1,maxWorkSize)

outIndexes <- integer(maxWorkSize+1)
for (i in 1:maxWorkSize) outIndexes[i] = distStart(i,maxWorkSize)
outIndexes[maxWorkSize+1] = outIndexes[maxWorkSize]

inBuffer <- createBuffer(context,"CL_MEM_READ_ONLY",length(dvec),dvec)
indexBuffer <- createBuffer(context,"CL_MEM_READ_ONLY",length(outIndexes),outIndexes)
outBuffer <- createBufferFloatVector(context,"CL_MEM_WRITE_ONLY",distSizeBlock)

# Euclidean Distance Kernel Function
kernelCode1 <- '
__kernel void kernelDist1(/*1*/__global const float *input,
  /*2*/__global const int *indexes, /*3*/__global float *output,
  /*4*/int numObs, /*5*/int numVars, /*6*/int startObs, /*7*/int stopObs)
{
  int startIndex = get_global_id(0) + startObs;
  if (startIndex >= stopObs) return;
  __global float *sptr = &input[startIndex * numVars];
  __global float *aptr;
  __global float *bptr = sptr + numVars;
  int distIndex = indexes[startIndex] - indexes[startObs]; 
  __global float *optr = &output[distIndex];
  int obsIndex; int i;
  float sum; float diff;
  for (obsIndex = startIndex+1; obsIndex < numObs; obsIndex++, optr++)
  {
    aptr = sptr; sum = 0.0;
    for (i = 0; i < numVars; i++, aptr++, bptr++)
    {
      diff = *aptr - *bptr;
      sum += diff * diff;
    }
    *optr = sqrt(sum);
  }
}
'
kernelCode2 <- '
__kernel void kernelDist2(/*1*/__global const float *input,
  /*2*/__global const int *indexes, /*3*/__global float *output,
  /*4*/int numObs, /*5*/int numVars, /*6*/int startObs, /*7*/int stopObs)
{
  int startIndex = get_global_id(0) + startObs;
  if (startIndex >= stopObs) return;
  __global float *sptr = &input[startIndex * numVars];
  __global float *bptr = sptr + numVars;
  int distIndex = indexes[startIndex] - indexes[startObs]; 
  __global float *optr = &output[distIndex];
  int obsIndex;
  float sum;
  float16 a, b, d, d2;
  a = vload16(0,sptr);
  for (obsIndex = startIndex+1; obsIndex < numObs; obsIndex++, optr++, bptr += numVars)
  {
    b = vload16(0,bptr);
    d = a - b;
    d2 = d * d;
    sum = d2.s0 + d2.s1 + d2.s2 + d2.s3 + d2.s4 + d2.s5 + d2.s6 + d2.s7 + 
          d2.s8 + d2.s9 + d2.sA + d2.sB + d2.sC + d2.sD + d2.sE;
    *optr = sqrt(sum);
  }
}
'
kernel1 <- buildKernel(context,kernelCode1,'kernelDist1', # Name has to match that declared in C
                       inBuffer,indexBuffer,outBuffer,
                       as.integer(maxWorkSize),as.integer(15),
                       as.integer(0),as.integer(blockWorkSize))
kernel2 <- buildKernel(context,kernelCode2,'kernelDist2', # Name has to match that declared in C
                       inBuffer,indexBuffer,outBuffer,
                       as.integer(maxWorkSize),as.integer(15),
                       as.integer(0),as.integer(blockWorkSize))

enqueueWriteBuffer(queue,inBuffer,length(dvec),dvec)
enqueueWriteBuffer(queue,indexBuffer,length(outIndexes),outIndexes)

result <- numeric(distSizeMax) # Large memory allocation can take several seconds

gpuDist <- function(kernel,result)
{
  numBlocks <- globalWorkSize / blockWorkSize
  remainderWorkSize = maxWorkSize
  obsIndex <- 1
  for (b in 1:numBlocks)
  {
    workSize <- blockWorkSize
    if (remainderWorkSize < workSize) workSize <- remainderWorkSize
    
    kernelStartObs <- obsIndex-1 # Index mapped from R:1..n to C:0..n-1
    kernelStopObs <- kernelStartObs + workSize
    assignKernelArg(kernel,6,as.integer(kernelStartObs))
    assignKernelArg(kernel,7,as.integer(kernelStopObs))
    
    # For NDRangeKernel execution, globalWorkSize must be a multiple of localWorkSize
    enqueueNDRangeKernel(queue,kernel,blockWorkSize,localWorkSize)
    
    distOffset <- outIndexes[obsIndex]
    distSize <- outIndexes[obsIndex + workSize] - distOffset
    enqueueReadBuffer(queue,outBuffer,distSize,result,distOffset)
    
    obsIndex <- obsIndex + workSize
    remainderWorkSize <- remainderWorkSize - workSize
  }
  return (TRUE)
}

#teval(gpuDist(kernel1,result))
teval(gpuDist(kernel2,result))

##
## Packt: "Mastering Parallelism with R"
## Chapter 6 - ROpenCL Distance Function
##
## Copyright 2015 Simon Chapple
##