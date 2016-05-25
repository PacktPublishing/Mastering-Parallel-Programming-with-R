#ifndef _ROpenCL_CREATECONTEXT_H
#define _ROpenCL_CREATECONTEXT_H

#include <CL/opencl.h>
#include <Rcpp.h>

/*
 * note : RcppExport is an alias to `extern "C"` defined by Rcpp.
 *
 * It gives C calling convention to the rcpp_hello_world function so that 
 * it can be called from .Call in R. Otherwise, the C++ compiler mangles the 
 * name of the function and .Call can't find it.
 *
 * It is only useful to use RcppExport when the function is intended to be called
 * by .Call. See the thread http://thread.gmane.org/gmane.comp.lang.r.rcpp/649/focus=672
 * on Rcpp-devel for a misuse of RcppExport
 */
RcppExport SEXP getPlatformIDs();
RcppExport SEXP getPlatformInfo(SEXP platformID);
RcppExport SEXP getPlatformName(SEXP sPlatformID);
RcppExport SEXP getDeviceIDs(SEXP platformID);
RcppExport SEXP getGPUDeviceIDs(SEXP platformID);
RcppExport SEXP getCPUDeviceIDs(SEXP platformID);
RcppExport SEXP getDeviceInfo(SEXP deviceID);
RcppExport SEXP getDeviceDoubleFP(SEXP deviceID);
RcppExport SEXP getDeviceSingleFP(SEXP deviceID);
RcppExport SEXP getDeviceHalfFP(SEXP deviceID);
RcppExport SEXP getDeviceType(SEXP sDeviceID);
RcppExport SEXP getKernelWorkGroupInfo(SEXP sKernel, SEXP sDeviceID);
RcppExport SEXP createContextFromType(SEXP deviceType);
RcppExport SEXP createCommandQueue(SEXP sContext, SEXP sDeviceID);
RcppExport SEXP createBufferFloatVector(SEXP sContext, SEXP sMemFlag, SEXP sGlobalWorkSize);
RcppExport SEXP createBufferIntegerVector(SEXP sContext, SEXP sMemFlag, SEXP sGlobalWorkSize);
RcppExport SEXP createProgramWithSource(SEXP sContext, SEXP sKernelSrc);
RcppExport SEXP buildProgram(SEXP sProgram);
RcppExport SEXP createKernel(SEXP sProgram, SEXP sKernelName);
RcppExport SEXP setKernelArgMem(SEXP sKernel, SEXP sIndex, SEXP sBuffer);
RcppExport SEXP setKernelArgInt(SEXP sKernel, SEXP sIndex, SEXP sIntegerValue);
RcppExport SEXP setKernelArgFloat(SEXP sKernel, SEXP sIndex, SEXP sFloatValue);
RcppExport void enqueueWriteBufferFloatVector(SEXP sQueue, SEXP sMemBuffer, SEXP sGlobalWorkSize, SEXP sObject);
RcppExport void enqueueWriteBufferIntegerVector(SEXP sQueue, SEXP sMemBuffer, SEXP sGlobalWorkSize, SEXP sObject);
RcppExport void enqueue1DRangeKernel(SEXP sQueue, SEXP sKernel, SEXP sGlobalWorkSize, SEXP sLocalWorkSize);
RcppExport void enqueue2DRangeKernel(SEXP sQueue, SEXP sKernel, SEXP sGlobalWorkSize, SEXP sLocalWorkSize);
RcppExport SEXP enqueueReadBufferFloatVector(SEXP sQueue, SEXP sMemBuffer, SEXP offset, SEXP sGlobalWorkSize, SEXP sObject);
RcppExport SEXP enqueueReadBufferIntegerVector(SEXP sQueue, SEXP sMemBuffer, SEXP offset, SEXP sGlobalWorkSize, SEXP sObject);
RcppExport SEXP getProgramInfo(SEXP sProgram, SEXP sProgramInfo);
RcppExport void getContextInfo(SEXP sContext, SEXP sContextInfo);
RcppExport SEXP createContext(SEXP sDeviceID);


/*
static void intObjFinalizer(SEXP ref);
SEXP int2EXP(int *o);
int *SEXP2int(SEXP o);
SEXP getIntPointer();
SEXP doubleIntPointer(SEXP test);*/

#endif
