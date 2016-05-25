#include <CL/opencl.h>
#include <Rcpp.h>
#include "createContext.h"

RCPP_MODULE(ropencl){
	using namespace Rcpp ;
	function( "getPlatformIDs"  , &getPlatformIDs   , "getPlatformIDs help" ) ;
    function( "getPlatformName"  , &getPlatformName   , "getPlatformName help" ) ;
    function( "getPlatformInfo"  , &getPlatformInfo   , "getPlatformInfo help" ) ;
    function( "getDeviceIDs"  , &getDeviceIDs   , "getDeviceIDs help" ) ;
    function( "getGPUDeviceIDs"  , &getGPUDeviceIDs   , "getGPUDeviceIDs help" ) ;
    function( "getCPUDeviceIDs"  , &getCPUDeviceIDs   , "getCPUDeviceIDs help" ) ;
    function( "getDeviceType"  , &getDeviceType   , "getDeviceType help" ) ;
    function( "getDeviceInfo"  , &getDeviceInfo   , "getDeviceInfo help" ) ;
    function( "getDeviceDoubleFP"  , &getDeviceDoubleFP   , "getDeviceDoubleFP help" ) ;
    function( "getDeviceSingleFP"  , &getDeviceSingleFP   , "getDeviceSingleFP help" ) ;
    function( "getDeviceHalfFP"  , &getDeviceHalfFP   , "getDeviceHalfFP help" ) ;
    function( "getKernelWorkGroupInfo"  , &getKernelWorkGroupInfo   , "getKernelWorkGroupInfo help" ) ;
    function( "createContextFromType"  , &createContextFromType   , "createContextFromType help" ) ;
    function( "createCommandQueue"  , &createCommandQueue   , "createCommandQueue help" ) ;
    function( "createBufferFloatVector"  , &createBufferFloatVector   , "createBufferFloatVector help" ) ;
    function( "createBufferIntegerVector"  , &createBufferIntegerVector , "createBufferIntegerVector help" ) ;
    function( "createProgramWithSource"  , &createProgramWithSource   , "createProgramWithSource help" ) ;
    function( "buildProgram"  , &buildProgram   , "buildProgram help" ) ;
    function( "createKernel"  , &createKernel   , "createKernel help" ) ;
    function( "setKernelArgMem"  , &setKernelArgMem , "setKernelArgMem help" ) ;
    function( "setKernelArgInt"  , &setKernelArgInt , "setKernelArgInt help" ) ;
    function( "setKernelArgFloat"  , &setKernelArgFloat , "setKernelArgFloat help" ) ;
    function( "enqueueWriteBufferFloatVector"  , &enqueueWriteBufferFloatVector , "enqueueWriteBufferFloatVector help" ) ;
    function( "enqueueWriteBufferIntegerVector"  , &enqueueWriteBufferIntegerVector , "enqueueWriteBufferIntegerVector help" ) ;
    function( "enqueue1DRangeKernel"  , &enqueue1DRangeKernel , "enqueue1DRangeKernel help" ) ;
    function( "enqueue2DRangeKernel"  , &enqueue2DRangeKernel , "enqueue2DRangeKernel help" ) ;
    function( "enqueueReadBufferFloatVector"  , &enqueueReadBufferFloatVector , "enqueueReadBufferFloatVector help" ) ;
    function( "enqueueReadBufferIntegerVector"  , &enqueueReadBufferIntegerVector , "enqueueReadBufferIntegerVector help" ) ;
    function( "getProgramInfo"  , &getProgramInfo , "getProgramInfo help" ) ;
    function( "getContextInfo"  , &getContextInfo , "getContextInfo help" ) ;
    function( "createContext"  , &createContext, "createContext help" ) ;
}




