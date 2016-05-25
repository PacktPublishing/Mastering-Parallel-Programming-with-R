#include "createContext.h"
#include <iostream>
#include <CL/opencl.h>
#include <Rcpp.h>
#include <vector>
#include <string>
#include <sstream>


static cl_int last_error_code = CL_SUCCESS;
static const char *last_error_function = NULL;

static int logError(cl_int code)
{
    if (code == CL_SUCCESS) return CL_FALSE;
    
    last_error_code = code;
    
    return CL_TRUE;
}

static void resetError(const char *fn)
{
    last_error_code = CL_SUCCESS;
    last_error_function = fn;
}

static cl_int getLastError(void)
{
    return last_error_code;
}

static const char *getLastFunction(void)
{
    return last_error_function;
}

/*
 * char *getErrorMessage(cl_int code)
 *
 * Maps the return value from a cl_API call to its error message equivalent.
 * code=0 is CL_SUCCESS so NULL is returned in this case.
 * If code is negative but not a recognised error code in terms of the compilation constants defined here,
 * then the string "ERROR%d" is returned.
 */
static const char *getErrorMessage(cl_int code)
{
    static char cbuf[64]; // small constant buffer memory for dynamic error message.
    const char *msg = NULL;
    cl_int unknown = CL_FALSE;
    
    switch (code)
    {
#ifdef CL_SUCCESS
        case CL_SUCCESS: return NULL;
#endif
#ifdef CL_DEVICE_NOT_FOUND
        case CL_DEVICE_NOT_FOUND: msg = "DEVICE_NOT_FOUND"; break;
#endif
#ifdef CL_DEVICE_NOT_AVAILABLE
        case CL_DEVICE_NOT_AVAILABLE: msg = "DEVICE_NOT_AVAILABLE"; break;
#endif
#ifdef CL_COMPILER_NOT_AVAILABLE
        case CL_COMPILER_NOT_AVAILABLE: msg = "COMPILER_NOT_AVAILABLE"; break;
#endif
#ifdef CL_MEM_OBJECT_ALLOCATION_FAILURE
        case CL_MEM_OBJECT_ALLOCATION_FAILURE: msg = "MEM_OBJECT_ALLOCATION_FAILURE"; break;
#endif
#ifdef CL_OUT_OF_RESOURCES
        case CL_OUT_OF_RESOURCES: msg = "OUT_OF_RESOURCES"; break;
#endif
#ifdef CL_OUT_OF_HOST_MEMORY
        case CL_OUT_OF_HOST_MEMORY: msg = "OUT_OF_HOST_MEMORY"; break;
#endif
#ifdef CL_PROFILING_INFO_NOT_AVAILABLE
        case CL_PROFILING_INFO_NOT_AVAILABLE: msg = "PROFILING_INFO_NOT_AVAILABLE"; break;
#endif
#ifdef CL_MEM_COPY_OVERLAP
        case CL_MEM_COPY_OVERLAP: msg = "MEM_COPY_OVERLAP"; break;
#endif
#ifdef CL_IMAGE_FORMAT_MISMATCH
        case CL_IMAGE_FORMAT_MISMATCH: msg = "IMAGE_FORMAT_MISMATCH"; break;
#endif
#ifdef CL_IMAGE_FORMAT_NOT_SUPPORTED
        case CL_IMAGE_FORMAT_NOT_SUPPORTED: msg = "IMAGE_FORMAT_NOT_SUPPORTED"; break;
#endif
#ifdef CL_BUILD_PROGRAM_FAILURE
        case CL_BUILD_PROGRAM_FAILURE: msg = "BUILD_PROGRAM_FAILURE"; break;
#endif
#ifdef CL_MAP_FAILURE
        case CL_MAP_FAILURE: msg = "MAP_FAILURE"; break;
#endif
#ifdef CL_MISALIGNED_SUB_BUFFER_OFFSET
        case CL_MISALIGNED_SUB_BUFFER_OFFSET: msg = "MISALIGNED_SUB_BUFFER_OFFSET"; break;
#endif
#ifdef CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST
        case CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST: msg = "EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST"; break;
#endif
#ifdef CL_COMPILE_PROGRAM_FAILURE
        case CL_COMPILE_PROGRAM_FAILURE: msg = "COMPILE_PROGRAM_FAILURE"; break;
#endif
#ifdef CL_LINKER_NOT_AVAILABLE
        case CL_LINKER_NOT_AVAILABLE: msg = "LINKER_NOT_AVAILABLE"; break;
#endif
#ifdef CL_LINK_PROGRAM_FAILURE
        case CL_LINK_PROGRAM_FAILURE: msg = "LINK_PROGRAM_FAILURE"; break;
#endif
#ifdef CL_DEVICE_PARTITION_FAILED
        case CL_DEVICE_PARTITION_FAILED: msg = "DEVICE_PARTITION_FAILED"; break;
#endif
#ifdef CL_KERNEL_ARG_INFO_NOT_AVAILABLE
        case CL_KERNEL_ARG_INFO_NOT_AVAILABLE: msg = "KERNEL_ARG_INFO_NOT_AVAILABLE"; break;
#endif
#ifdef CL_INVALID_VALUE
        case CL_INVALID_VALUE: msg = "INVALID_VALUE"; break;
#endif
#ifdef CL_INVALID_DEVICE_TYPE
        case CL_INVALID_DEVICE_TYPE: msg = "INVALID_DEVICE_TYPE"; break;
#endif
#ifdef CL_INVALID_PLATFORM
        case CL_INVALID_PLATFORM: msg = "INVALID_PLATFORM"; break;
#endif
#ifdef CL_INVALID_DEVICE
        case CL_INVALID_DEVICE: msg = "INVALID_DEVICE"; break;
#endif
#ifdef CL_INVALID_CONTEXT
        case CL_INVALID_CONTEXT: msg = "INVALID_CONTEXT"; break;
#endif
#ifdef CL_INVALID_QUEUE_PROPERTIES
        case CL_INVALID_QUEUE_PROPERTIES: msg = "INVALID_QUEUE_PROPERTIES"; break;
#endif
#ifdef CL_INVALID_COMMAND_QUEUE
        case CL_INVALID_COMMAND_QUEUE: msg = "INVALID_COMMAND_QUEUE"; break;
#endif
#ifdef CL_INVALID_HOST_PTR
        case CL_INVALID_HOST_PTR: msg = "INVALID_HOST_PTR"; break;
#endif
#ifdef CL_INVALID_MEM_OBJECT
        case CL_INVALID_MEM_OBJECT: msg = "INVALID_MEM_OBJECT"; break;
#endif
#ifdef CL_INVALID_IMAGE_FORMAT_DESCRIPTOR
        case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR: msg = "INVALID_IMAGE_FORMAT_DESCRIPTOR"; break;
#endif
#ifdef CL_INVALID_IMAGE_SIZE
        case CL_INVALID_IMAGE_SIZE: msg = "INVALID_IMAGE_SIZE"; break;
#endif
#ifdef CL_INVALID_SAMPLER
        case CL_INVALID_SAMPLER: msg = "INVALID_SAMPLER"; break;
#endif
#ifdef CL_INVALID_BINARY
        case CL_INVALID_BINARY: msg = "INVALID_BINARY"; break;
#endif
#ifdef CL_INVALID_BUILD_OPTIONS
        case CL_INVALID_BUILD_OPTIONS: msg = "INVALID_BUILD_OPTIONS"; break;
#endif
#ifdef CL_INVALID_PROGRAM
        case CL_INVALID_PROGRAM: msg = "INVALID_PROGRAM"; break;
#endif
#ifdef CL_INVALID_PROGRAM_EXECUTABLE
        case CL_INVALID_PROGRAM_EXECUTABLE: msg = "INVALID_PROGRAM_EXECUTABLE"; break;
#endif
#ifdef CL_INVALID_KERNEL_NAME
        case CL_INVALID_KERNEL_NAME: msg = "INVALID_KERNEL_NAME"; break;
#endif
#ifdef CL_INVALID_KERNEL_DEFINITION
        case CL_INVALID_KERNEL_DEFINITION: msg = "INVALID_KERNEL_DEFINITION"; break;
#endif
#ifdef CL_INVALID_KERNEL
        case CL_INVALID_KERNEL: msg = "INVALID_KERNEL"; break;
#endif
#ifdef CL_INVALID_ARG_INDEX
        case CL_INVALID_ARG_INDEX: msg = "INVALID_ARG_INDEX"; break;
#endif
#ifdef CL_INVALID_ARG_VALUE
        case CL_INVALID_ARG_VALUE: msg = "INVALID_ARG_VALUE"; break;
#endif
#ifdef CL_INVALID_ARG_SIZE
        case CL_INVALID_ARG_SIZE: msg = "INVALID_ARG_SIZE"; break;
#endif
#ifdef CL_INVALID_KERNEL_ARGS
        case CL_INVALID_KERNEL_ARGS: msg = "INVALID_KERNEL_ARGS"; break;
#endif
#ifdef CL_INVALID_WORK_DIMENSION
        case CL_INVALID_WORK_DIMENSION: msg = "INVALID_WORK_DIMENSION"; break;
#endif
#ifdef CL_INVALID_WORK_GROUP_SIZE
        case CL_INVALID_WORK_GROUP_SIZE: msg = "INVALID_WORK_GROUP_SIZE"; break;
#endif
#ifdef CL_INVALID_WORK_ITEM_SIZE
        case CL_INVALID_WORK_ITEM_SIZE: msg = "INVALID_WORK_ITEM_SIZE"; break;
#endif
#ifdef CL_INVALID_GLOBAL_OFFSET
        case CL_INVALID_GLOBAL_OFFSET: msg = "INVALID_GLOBAL_OFFSET"; break;
#endif
#ifdef CL_INVALID_EVENT_WAIT_LIST
        case CL_INVALID_EVENT_WAIT_LIST: msg = "INVALID_EVENT_WAIT_LIST"; break;
#endif
#ifdef CL_INVALID_EVENT
        case CL_INVALID_EVENT: msg = "INVALID_EVENT"; break;
#endif
#ifdef CL_INVALID_OPERATION
        case CL_INVALID_OPERATION: msg = "INVALID_OPERATION"; break;
#endif
#ifdef CL_INVALID_GL_OBJECT
        case CL_INVALID_GL_OBJECT: msg = "INVALID_GL_OBJECT"; break;
#endif
#ifdef CL_INVALID_BUFFER_SIZE
        case CL_INVALID_BUFFER_SIZE: msg = "INVALID_BUFFER_SIZE"; break;
#endif
#ifdef CL_INVALID_MIP_LEVEL
        case CL_INVALID_MIP_LEVEL: msg = "INVALID_MIP_LEVEL"; break;
#endif
#ifdef CL_INVALID_GLOBAL_WORK_SIZE
        case CL_INVALID_GLOBAL_WORK_SIZE: msg = "INVALID_GLOBAL_WORK_SIZE"; break;
#endif
#ifdef CL_INVALID_PROPERTY
        case CL_INVALID_PROPERTY: msg = "INVALID_PROPERTY"; break;
#endif
#ifdef CL_INVALID_IMAGE_DESCRIPTOR
        case CL_INVALID_IMAGE_DESCRIPTOR: msg = "INVALID_IMAGE_DESCRIPTOR"; break;
#endif
#ifdef CL_INVALID_COMPILER_OPTIONS
        case CL_INVALID_COMPILER_OPTIONS: msg = "INVALID_COMPILER_OPTIONS"; break;
#endif
#ifdef CL_INVALID_LINKER_OPTIONS
        case CL_INVALID_LINKER_OPTIONS: msg = "INVALID_LINKER_OPTIONS"; break;
#endif
#ifdef CL_INVALID_DEVICE_PARTITION_COUNT
        case CL_INVALID_DEVICE_PARTITION_COUNT: msg = "INVALID_DEVICE_PARTITION_COUNT"; break;
#endif
        default:
            unknown = CL_TRUE;
            sprintf(cbuf,"%d",code);
            msg = cbuf;
    }
    
    //if (unknown) return msg;
    //else return "CL_" + msg;
    return msg;
}

/*
 * int printOnError(code,preamble)
 *
 * Currently this function prints an error message to stdout if code != CL_SUCCESS.
 * TODO: errors ought to be handled better than this...
 *
 * Function returns 1 if it was an error code, 0 otherwise.
 */
static int printOnError(cl_int code, const char *preamble)
{
    if (code == CL_SUCCESS) return CL_FALSE;
    
    if (preamble != NULL) { std::cout << preamble << "\n"; };
    std::cout << "ERROR: " << getErrorMessage(code) << "\n";
    
    return CL_TRUE;
}

#define PRINTONERROR() printOnError(getLastError(),getLastFunction())


SEXP getPlatformIDs(){
    //returns a list of platform ids
    resetError("ROpenCL::getPlatformIDs");
    cl_uint num_platforms = 0;
    logError( clGetPlatformIDs(0, 0, &num_platforms) );
    std::vector<cl_platform_id> platforms(num_platforms);
    logError( clGetPlatformIDs(num_platforms, platforms.empty() ? NULL : &platforms.front(), &num_platforms) );
    //for each platform in platforms add its pointer to the return list
    Rcpp::List result(platforms.size());
    for (int i=0; i<platforms.size(); i++){
        cl_platform_id *tempPlatformID = new cl_platform_id;
        *tempPlatformID = platforms[i];
         Rcpp::XPtr<cl_platform_id> tempXPtr(tempPlatformID);
        result[i] = tempXPtr;
    }
    PRINTONERROR();
    return result;
}


SEXP getPlatformName(SEXP sPlatformID){
    resetError("ROpenCL::getPlatformName");
    char cBuffer[1024];
    Rcpp::XPtr<cl_platform_id> platformID(sPlatformID);
    cBuffer[0] = '\0';
    logError( clGetPlatformInfo(*platformID, CL_PLATFORM_NAME, sizeof(cBuffer), cBuffer, NULL) );
    std::string retVal = cBuffer;
    PRINTONERROR();
    return Rcpp::wrap(retVal);
}


SEXP getDeviceIDs(SEXP sPlatformID){
    resetError("ROpenCL::getDeviceIDs");
    Rcpp::XPtr<cl_platform_id> platformID(sPlatformID);

    cl_uint num_GPUs = 0;
    logError( clGetDeviceIDs(*platformID, CL_DEVICE_TYPE_GPU, 0, 0, &num_GPUs) );
    std::vector<cl_device_id> gpus(num_GPUs);
    logError( clGetDeviceIDs(*platformID, CL_DEVICE_TYPE_GPU, num_GPUs, gpus.empty() ? NULL : &gpus.front(), &num_GPUs) );

    cl_uint num_nonGPUs = 0;
    logError( clGetDeviceIDs(*platformID, ~CL_DEVICE_TYPE_GPU, 0, 0, &num_nonGPUs) );
    std::vector<cl_device_id> devices(num_nonGPUs);
    logError( clGetDeviceIDs(*platformID, ~CL_DEVICE_TYPE_GPU, num_nonGPUs, devices.empty() ? NULL : &devices.front(), &num_nonGPUs) );

    //for each device in devices add its pointer to the return list
    Rcpp::List result(gpus.size() + devices.size());
    //To be more compatible with the previous version of this function, put GPUs first in the result list.
    for (int i=0; i<gpus.size(); i++){
        cl_device_id *tempDeviceID = new cl_device_id;
        *tempDeviceID = gpus[i];
        Rcpp::XPtr<cl_device_id> tempXPtr(tempDeviceID);
        result[i] = tempXPtr;
    }
    for (int i=0; i<devices.size(); i++){
        cl_device_id *tempDeviceID = new cl_device_id;
        *tempDeviceID = devices[i];
        Rcpp::XPtr<cl_device_id> tempXPtr(tempDeviceID);
        result[i+gpus.size()] = tempXPtr;
    }
    PRINTONERROR();
    return result;
}


SEXP getGPUDeviceIDs(SEXP sPlatformID){
    resetError("ROpenCL::getGPUDeviceIDs");
    Rcpp::XPtr<cl_platform_id> platformID(sPlatformID);
    cl_uint num_devices = 0;

    logError( clGetDeviceIDs(*platformID, CL_DEVICE_TYPE_GPU, 0, 0, &num_devices) );
    std::vector<cl_device_id> devices(num_devices);
    logError( clGetDeviceIDs(*platformID, CL_DEVICE_TYPE_GPU, num_devices, devices.empty() ? NULL : &devices.front(), &num_devices) );
    //for each platform in platforms add its pointer to the return list
    Rcpp::List result(devices.size());
    for (int i=0; i<devices.size(); i++){
        cl_device_id *tempDeviceID = new cl_device_id;
        *tempDeviceID = devices[i];
        Rcpp::XPtr<cl_device_id> tempXPtr(tempDeviceID);
        result[i] = tempXPtr;
    }
    PRINTONERROR();
    return result;
}


SEXP getCPUDeviceIDs(SEXP sPlatformID){
    resetError("ROpenCL::getCPUDeviceIDs");
    Rcpp::XPtr<cl_platform_id> platformID(sPlatformID);
    cl_uint num_devices = 0;

    logError( clGetDeviceIDs(*platformID, CL_DEVICE_TYPE_CPU, 0, 0, &num_devices) );
    std::vector<cl_device_id> devices(num_devices);
    logError( clGetDeviceIDs(*platformID, CL_DEVICE_TYPE_CPU, num_devices, devices.empty() ? NULL : &devices.front(), &num_devices) );
    //for each platform in platforms add its pointer to the return list
    Rcpp::List result(devices.size());
    for (int i=0; i<devices.size(); i++){
        cl_device_id *tempDeviceID = new cl_device_id;
        *tempDeviceID = devices[i];
        Rcpp::XPtr<cl_device_id> tempXPtr(tempDeviceID);
        result[i] = tempXPtr;
    }
    PRINTONERROR();
    return result;
}


SEXP getPlatformInfo(SEXP sPlatformID){
    resetError("ROpenCL::getPlatformInfo");
    static char cBuffer[1024];
    Rcpp::XPtr<cl_platform_id> platformID(sPlatformID);
    std::string str;
    Rcpp::List result;
    
#ifdef CL_PLATFORM_PROFILE /* TODO CHECK THIS IS STRING RETURN VALUE */
    cBuffer[0] = '\0';
    logError( clGetPlatformInfo(*platformID, CL_PLATFORM_PROFILE, sizeof(cBuffer), cBuffer, NULL) );
    str = cBuffer;
    result["CL_PLATFORM_PROFILE"] = Rcpp::wrap(str);
#endif
#ifdef CL_PLATFORM_VERSION
    cBuffer[0] = '\0';
    logError( clGetPlatformInfo(*platformID, CL_PLATFORM_VERSION, sizeof(cBuffer), cBuffer, NULL) );
    str = cBuffer;
    result["CL_PLATFORM_VERSION"] = Rcpp::wrap(str);
#endif
#ifdef CL_PLATFORM_NAME
    cBuffer[0] = '\0';
    logError( clGetPlatformInfo (*platformID, CL_PLATFORM_NAME, sizeof(cBuffer), cBuffer, NULL) );
    str = cBuffer;
    result["CL_PLATFORM_NAME"] = Rcpp::wrap(str);
#endif
#ifdef CL_PLATFORM_VENDOR
    cBuffer[0] = '\0';
    logError( clGetPlatformInfo(*platformID, CL_PLATFORM_VENDOR, sizeof(cBuffer), cBuffer, NULL) );
    str = cBuffer;
    result["CL_PLATFORM_VENDOR"] = Rcpp::wrap(str);
#endif
#ifdef CL_PLATFORM_EXTENSIONS
    cBuffer[0] = '\0';
    logError( clGetPlatformInfo(*platformID, CL_PLATFORM_EXTENSIONS, sizeof(cBuffer), cBuffer, NULL) );
    str = cBuffer;
    std::string word;
    std::stringstream ss(str); // Insert the string into a stream
    Rcpp::List extensions;
    while (ss >> word) { // Fetch next word from stream
        extensions.push_back(Rcpp::wrap(word));
    }
    result["CL_PLATFORM_EXTENSIONS"] = Rcpp::wrap(extensions);
#endif

    PRINTONERROR();
    return result;
}


SEXP getDeviceType(SEXP sDeviceID){
    resetError("ROpenCL::getDeviceType");
    Rcpp::XPtr<cl_device_id> deviceID(sDeviceID);
    std::string str = "UNKNOWN";
    
#ifdef CL_DEVICE_TYPE
    cl_device_type dt = 0;
    logError( clGetDeviceInfo(*deviceID,CL_DEVICE_TYPE, sizeof(cl_device_type), &dt, NULL) );
# ifdef CL_DEVICE_TYPE_DEFAULT
    if (dt & CL_DEVICE_TYPE_DEFAULT) str = "CL_DEVICE_TYPE_DEFAULT";
# endif
# ifdef CL_DEVICE_TYPE_CPU
    if (dt & CL_DEVICE_TYPE_CPU) str = "CL_DEVICE_TYPE_CPU";
# endif
# ifdef CL_DEVICE_TYPE_GPU
    if (dt & CL_DEVICE_TYPE_GPU) str = "CL_DEVICE_TYPE_GPU";
# endif
# ifdef CL_DEVICE_TYPE_ACCELERATOR
    if (dt & CL_DEVICE_TYPE_ACCELERATOR) str = "CL_DEVICE_TYPE_ACCELERATOR";
# endif
# ifdef CL_DEVICE_TYPE_CUSTOM
    if (dt & CL_DEVICE_TYPE_CUSTOM) str = "CL_DEVICE_TYPE_CUSTOM";
# endif
#endif

    PRINTONERROR();
    return Rcpp::wrap(str);
}


static Rcpp::List getDeviceFPConfigAsList(cl_device_fp_config flags)
{
    Rcpp::List list;
#ifdef CL_FP_DENORM
    if (flags & CL_FP_DENORM) list.push_back(Rcpp::wrap("CL_FP_DENORM"));
#endif
#ifdef CL_FP_INF_NAN
    if (flags & CL_FP_INF_NAN) list.push_back(Rcpp::wrap("CL_FP_INF_NAN"));
#endif
#ifdef CL_FP_ROUND_TO_NEAREST
    if (flags & CL_FP_ROUND_TO_NEAREST) list.push_back(Rcpp::wrap("CL_FP_ROUND_TO_NEAREST"));
#endif
#ifdef CL_FP_ROUND_TO_ZERO
    if (flags & CL_FP_ROUND_TO_ZERO) list.push_back(Rcpp::wrap("CL_FP_ROUND_TO_ZERO"));
#endif
#ifdef CL_FP_ROUND_TO_INF
    if (flags & CL_FP_ROUND_TO_INF) list.push_back(Rcpp::wrap("CL_FP_ROUND_TO_INF"));
#endif
#ifdef CL_FP_FMA
    if (flags & CL_FP_FMA) list.push_back(Rcpp::wrap("CL_FP_FMA"));
#endif
#ifdef CL_FP_SOFT_FLOAT
    if (flags & CL_FP_SOFT_FLOAT) list.push_back(Rcpp::wrap("CL_FP_SOFT_FLOAT"));
#endif
#ifdef CL_FP_CORRECTLY_ROUNDED_DIVIDE_SQRT
    if (flags & CL_FP_CORRECTLY_ROUNDED_DIVIDE_SQRT) list.push_back(Rcpp::wrap("CL_FP_CORRECTLY_ROUNDED_DIVIDE_SQRT"));
#endif
    return list;
}

static Rcpp::List getDeviceMemCacheTypeAsList(cl_device_mem_cache_type flags)
{
    Rcpp::List list;
#ifdef CL_READ_ONLY_CACHE
    if (flags & CL_READ_ONLY_CACHE) list.push_back(Rcpp::wrap("CL_READ_ONLY_CACHE"));
#endif
#ifdef CL_READ_WRITE_CACHE
    if (flags & CL_READ_WRITE_CACHE) list.push_back(Rcpp::wrap("CL_READ_WRITE_CACHE"));
#endif
    return list;
}

SEXP getDeviceDoubleFP(SEXP sDeviceID){ /* 64-bit floating-point */
    Rcpp::XPtr<cl_device_id> deviceID(sDeviceID);
    Rcpp::List result;
#ifdef CL_DEVICE_DOUBLE_FP_CONFIG
    cl_device_fp_config fpc = 0;
    if ( CL_SUCCESS == clGetDeviceInfo(*deviceID,CL_DEVICE_DOUBLE_FP_CONFIG, sizeof(cl_device_fp_config), &fpc, NULL) ) {
        result = getDeviceFPConfigAsList(fpc);
    }
#endif
    return Rcpp::wrap(result);
}

SEXP getDeviceSingleFP(SEXP sDeviceID){ /* 32-bit floating-point */
    Rcpp::XPtr<cl_device_id> deviceID(sDeviceID);
    Rcpp::List result;
#ifdef CL_DEVICE_SINGLE_FP_CONFIG
    cl_device_fp_config fpc = 0;
    if ( CL_SUCCESS == clGetDeviceInfo(*deviceID,CL_DEVICE_SINGLE_FP_CONFIG, sizeof(cl_device_fp_config), &fpc, NULL) ) {
        result = getDeviceFPConfigAsList(fpc);
    }
#endif
    return Rcpp::wrap(result);
}

SEXP getDeviceHalfFP(SEXP sDeviceID){ /* 16-bit floating-point */
    Rcpp::XPtr<cl_device_id> deviceID(sDeviceID);
    Rcpp::List result;
#ifdef CL_DEVICE_HALF_FP_CONFIG
    cl_device_fp_config fpc = 0;
    if ( CL_SUCCESS == clGetDeviceInfo(*deviceID,CL_DEVICE_HALF_FP_CONFIG, sizeof(cl_device_fp_config), &fpc, NULL) ) {
        result = getDeviceFPConfigAsList(fpc);
    }
#endif
    return Rcpp::wrap(result);
}

SEXP getDeviceInfo(SEXP sDeviceID){
    resetError("ROpenCL::getDeviceInfo");
    static char cBuffer[1024];
    Rcpp::XPtr<cl_device_id> deviceID(sDeviceID);
    std::string str;
    cl_uint uival;
    cl_ulong ulval;
    size_t sval;
    Rcpp::List result;
    
#ifdef CL_DEVICE_TYPE
    cl_device_type dt = 0;
    logError( clGetDeviceInfo(*deviceID,CL_DEVICE_TYPE, sizeof(cl_device_type), &dt, NULL) );
    str = "UNKNOWN";
# ifdef CL_DEVICE_TYPE_DEFAULT
    if (dt & CL_DEVICE_TYPE_DEFAULT) str = "CL_DEVICE_TYPE_DEFAULT";
# endif
# ifdef CL_DEVICE_TYPE_CPU
    if (dt & CL_DEVICE_TYPE_CPU) str = "CL_DEVICE_TYPE_CPU";
# endif
# ifdef CL_DEVICE_TYPE_GPU
    if (dt & CL_DEVICE_TYPE_GPU) str = "CL_DEVICE_TYPE_GPU";
# endif
# ifdef CL_DEVICE_TYPE_ACCELERATOR
    if (dt & CL_DEVICE_TYPE_ACCELERATOR) str = "CL_DEVICE_TYPE_ACCELERATOR";
# endif
# ifdef CL_DEVICE_TYPE_CUSTOM
    if (dt & CL_DEVICE_TYPE_CUSTOM) str = "CL_DEVICE_TYPE_CUSTOM";
# endif
    result["CL_DEVICE_TYPE"] = Rcpp::wrap(str);
#endif
#ifdef CL_DEVICE_VENDOR_ID
    uival = 0;
    logError( clGetDeviceInfo(*deviceID, CL_DEVICE_VENDOR_ID, sizeof(uival), &uival, NULL) );
    result["CL_DEVICE_VENDOR_ID"] = Rcpp::wrap(uival);
#endif
#ifdef CL_DEVICE_MAX_COMPUTE_UNITS
    uival = 0;
    logError( clGetDeviceInfo(*deviceID, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(uival), &uival, NULL) );
    result["CL_DEVICE_MAX_COMPUTE_UNITS"] = Rcpp::wrap(uival);
#endif
#ifdef CL_DEVICE_MAX_WORK_GROUP_SIZE
    sval = 0;
    logError( clGetDeviceInfo(*deviceID, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(sval), &sval, NULL) );
    result["CL_DEVICE_MAX_WORK_GROUP_SIZE"] = Rcpp::wrap(sval);
#endif
#ifdef CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS
    uival = 0;
    logError( clGetDeviceInfo(*deviceID, CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, sizeof(uival), &uival, NULL) );
    result["CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS"] = Rcpp::wrap(uival);
    size_t *sizes = (size_t *) malloc(sizeof(size_t) * uival);
# ifdef CL_DEVICE_MAX_WORK_ITEM_SIZES
    logError( clGetDeviceInfo(*deviceID, CL_DEVICE_MAX_WORK_ITEM_SIZES, sizeof(size_t) * uival, sizes, NULL) );
    Rcpp::List wis;
    for (int i = 0; i < uival; i++) {
        int size = sizes[i];
        wis.push_back(Rcpp::wrap(size));
    }
    result["CL_DEVICE_MAX_WORK_ITEM_SIZES"] = Rcpp::wrap(wis);
# endif
    free(sizes);
#endif
#ifdef CL_DEVICE_PREFERRED_VECTOR_WIDTH_CHAR
    uival = 0;
    logError( clGetDeviceInfo(*deviceID, CL_DEVICE_PREFERRED_VECTOR_WIDTH_CHAR, sizeof(uival), &uival, NULL) );
    result["CL_DEVICE_PREFERRED_VECTOR_WIDTH_CHAR"] = Rcpp::wrap(uival);
#endif
#ifdef CL_DEVICE_PREFERRED_VECTOR_WIDTH_SHORT
    uival = 0;
    logError( clGetDeviceInfo(*deviceID, CL_DEVICE_PREFERRED_VECTOR_WIDTH_SHORT, sizeof(uival), &uival, NULL) );
    result["CL_DEVICE_PREFERRED_VECTOR_WIDTH_SHORT"] = Rcpp::wrap(uival);
#endif
#ifdef CL_DEVICE_PREFERRED_VECTOR_WIDTH_INT
    uival = 0;
    logError( clGetDeviceInfo(*deviceID, CL_DEVICE_PREFERRED_VECTOR_WIDTH_INT, sizeof(uival), &uival, NULL) );
    result["CL_DEVICE_PREFERRED_VECTOR_WIDTH_INT"] = Rcpp::wrap(uival);
#endif
#ifdef CL_DEVICE_PREFERRED_VECTOR_WIDTH_LONG
    uival = 0;
    logError( clGetDeviceInfo(*deviceID, CL_DEVICE_PREFERRED_VECTOR_WIDTH_LONG, sizeof(uival), &uival, NULL) );
    result["CL_DEVICE_PREFERRED_VECTOR_WIDTH_LONG"] = Rcpp::wrap(uival);
#endif
#ifdef CL_DEVICE_PREFERRED_VECTOR_WIDTH_FLOAT
    uival = 0;
    logError( clGetDeviceInfo(*deviceID, CL_DEVICE_PREFERRED_VECTOR_WIDTH_FLOAT, sizeof(uival), &uival, NULL) );
    result["CL_DEVICE_PREFERRED_VECTOR_WIDTH_FLOAT"] = Rcpp::wrap(uival);
#endif
#ifdef CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE
    uival = 0;
    logError( clGetDeviceInfo(*deviceID, CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE, sizeof(uival), &uival, NULL) );
    result["CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE"] = Rcpp::wrap(uival);
#endif
#ifdef CL_DEVICE_MAX_CLOCK_FREQUENCY
    uival = 0;
    logError( clGetDeviceInfo(*deviceID, CL_DEVICE_MAX_CLOCK_FREQUENCY, sizeof(uival), &uival, NULL) );
    result["CL_DEVICE_MAX_CLOCK_FREQUENCY"] = Rcpp::wrap(uival);
#endif
#ifdef CL_DEVICE_ADDRESS_BITS
#endif
#ifdef CL_DEVICE_MAX_READ_IMAGE_ARGS
#endif
#ifdef CL_DEVICE_MAX_WRITE_IMAGE_ARGS
#endif
#ifdef CL_DEVICE_MAX_MEM_ALLOC_SIZE
    ulval = 0;
    logError( clGetDeviceInfo(*deviceID, CL_DEVICE_MAX_MEM_ALLOC_SIZE, sizeof(ulval), &ulval, NULL) );
    result["CL_DEVICE_MAX_MEM_ALLOC_SIZE"] = Rcpp::wrap(ulval);
#endif
#ifdef CL_DEVICE_IMAGE2D_MAX_WIDTH
#endif
#ifdef CL_DEVICE_IMAGE2D_MAX_HEIGHT
#endif
#ifdef CL_DEVICE_IMAGE3D_MAX_WIDTH
#endif
#ifdef CL_DEVICE_IMAGE3D_MAX_HEIGHT
#endif
#ifdef CL_DEVICE_IMAGE3D_MAX_DEPTH
#endif
#ifdef CL_DEVICE_IMAGE_SUPPORT
#endif
#ifdef CL_DEVICE_MAX_PARAMETER_SIZE
#endif
#ifdef CL_DEVICE_MEM_BASE_ADDR_ALIGN
#endif
#ifdef CL_DEVICE_MIN_DATA_TYPE_ALIGN_SIZE
#endif
#ifdef CL_DEVICE_SINGLE_FP_CONFIG
    result["CL_DEVICE_SINGLE_FP_CONFIG"] = getDeviceSingleFP(sDeviceID);
#endif
#ifdef CL_DEVICE_GLOBAL_MEM_CACHE_TYPE
#endif
#ifdef CL_DEVICE_GLOBAL_MEM_CACHELINE_SIZE
#endif
#ifdef CL_DEVICE_GLOBAL_MEM_CACHE_SIZE
    ulval = 0;
    logError( clGetDeviceInfo(*deviceID, CL_DEVICE_GLOBAL_MEM_CACHE_SIZE, sizeof(ulval), &ulval, NULL) );
    result["CL_DEVICE_GLOBAL_MEM_CACHE_SIZE"] = Rcpp::wrap(ulval);
#endif
#ifdef CL_DEVICE_GLOBAL_MEM_SIZE
    ulval = 0;
    logError( clGetDeviceInfo(*deviceID, CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(ulval), &ulval, NULL) );
    result["CL_DEVICE_GLOBAL_MEM_SIZE"] = Rcpp::wrap(ulval);
#endif
#ifdef CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE
    ulval = 0;
    logError( clGetDeviceInfo(*deviceID, CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE, sizeof(ulval), &ulval, NULL) );
    result["CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE"] = Rcpp::wrap(ulval);
#endif
#ifdef CL_DEVICE_MAX_CONSTANT_ARGS
#endif
#ifdef CL_DEVICE_LOCAL_MEM_TYPE
#endif
#ifdef CL_DEVICE_LOCAL_MEM_SIZE
    ulval = 0;
    logError( clGetDeviceInfo(*deviceID, CL_DEVICE_LOCAL_MEM_SIZE, sizeof(ulval), &ulval, NULL) );
    result["CL_DEVICE_LOCAL_MEM_SIZE"] = Rcpp::wrap(ulval);
#endif
#ifdef CL_DEVICE_ERROR_CORRECTION_SUPPORT
#endif
#ifdef CL_DEVICE_PROFILING_TIMER_RESOLUTION
#endif
#ifdef CL_DEVICE_ENDIAN_LITTLE
#endif
#ifdef CL_DEVICE_AVAILABLE
#endif
#ifdef CL_DEVICE_COMPILER_AVAILABLE
#endif
#ifdef CL_DEVICE_EXECUTION_CAPABILITIES
#endif
#ifdef CL_DEVICE_QUEUE_PROPERTIES
#endif
#ifdef CL_DEVICE_ERROR_CORRECTION_SUPPORT
#endif
#ifdef CL_DEVICE_NAME
    cBuffer[0] = '\0';
    logError( clGetDeviceInfo (*deviceID, CL_DEVICE_NAME, sizeof(cBuffer), cBuffer, NULL) );
    str = cBuffer;
    result["CL_DEVICE_NAME"] = Rcpp::wrap(str);
#endif
#ifdef CL_DEVICE_VENDOR
    cBuffer[0] = '\0';
    logError( clGetDeviceInfo(*deviceID, CL_DEVICE_VENDOR, sizeof(cBuffer), cBuffer, NULL) );
    str = cBuffer;
    result["CL_DEVICE_VENDOR"] = Rcpp::wrap(str);
#endif
#ifdef CL_DRIVER_VERSION
    cBuffer[0] = '\0';
    logError( clGetDeviceInfo(*deviceID, CL_DRIVER_VERSION, sizeof(cBuffer), cBuffer, NULL) );
    str = cBuffer;
    result["CL_DRIVER_VERSION"] = Rcpp::wrap(str);
#endif
#ifdef CL_DEVICE_PROFILE
#endif
#ifdef CL_DEVICE_VERSION
    cBuffer[0] = '\0';
    logError( clGetDeviceInfo(*deviceID, CL_DEVICE_VERSION, sizeof(cBuffer), cBuffer, NULL) );
    str = cBuffer;
    result["CL_DEVICE_VERSION"] = Rcpp::wrap(str);
#endif
#ifdef CL_DEVICE_EXTENSIONS
    cBuffer[0] = '\0';
    logError( clGetDeviceInfo(*deviceID, CL_DEVICE_EXTENSIONS, sizeof(cBuffer), cBuffer, NULL) );
    str = cBuffer;
    std::string word;
    std::stringstream ss(str); // Insert the string into a stream
    Rcpp::List extensions;
    while (ss >> word) { // Fetch next word from stream
        extensions.push_back(Rcpp::wrap(word));
    }
    result["CL_DEVICE_EXTENSIONS"] = Rcpp::wrap(extensions);
#endif
#ifdef CL_DEVICE_PLATFORM
#endif
#ifdef CL_DEVICE_DOUBLE_FP_CONFIG
    result["CL_DEVICE_DOUBLE_FP_CONFIG"] = getDeviceDoubleFP(sDeviceID);
#endif
#ifdef CL_DEVICE_HALF_FP_CONFIG
    result["CL_DEVICE_HALF_FP_CONFIG"] = getDeviceHalfFP(sDeviceID);
#endif
#ifdef CL_DEVICE_PREFERRED_VECTOR_WIDTH_HALF
    uival = 0;
    logError( clGetDeviceInfo(*deviceID, CL_DEVICE_PREFERRED_VECTOR_WIDTH_HALF, sizeof(uival), &uival, NULL) );
    result["CL_DEVICE_PREFERRED_VECTOR_WIDTH_HALF"] = Rcpp::wrap(uival);
#endif
#ifdef CL_DEVICE_HOST_UNIFIED_MEMORY
#endif
#ifdef CL_DEVICE_NATIVE_VECTOR_WIDTH_CHAR
    uival = 0;
    logError( clGetDeviceInfo(*deviceID, CL_DEVICE_NATIVE_VECTOR_WIDTH_CHAR, sizeof(uival), &uival, NULL) );
    result["CL_DEVICE_NATIVE_VECTOR_WIDTH_CHAR"] = Rcpp::wrap(uival);
#endif
#ifdef CL_DEVICE_NATIVE_VECTOR_WIDTH_SHORT
    uival = 0;
    logError( clGetDeviceInfo(*deviceID, CL_DEVICE_NATIVE_VECTOR_WIDTH_SHORT, sizeof(uival), &uival, NULL) );
    result["CL_DEVICE_NATIVE_VECTOR_WIDTH_SHORT"] = Rcpp::wrap(uival);
#endif
#ifdef CL_DEVICE_NATIVE_VECTOR_WIDTH_INT
    uival = 0;
    logError( clGetDeviceInfo(*deviceID, CL_DEVICE_NATIVE_VECTOR_WIDTH_INT, sizeof(uival), &uival, NULL) );
    result["CL_DEVICE_NATIVE_VECTOR_WIDTH_INT"] = Rcpp::wrap(uival);
#endif
#ifdef CL_DEVICE_NATIVE_VECTOR_WIDTH_LONG
    uival = 0;
    logError( clGetDeviceInfo(*deviceID, CL_DEVICE_NATIVE_VECTOR_WIDTH_LONG, sizeof(uival), &uival, NULL) );
    result["CL_DEVICE_NATIVE_VECTOR_WIDTH_LONG"] = Rcpp::wrap(uival);
#endif
#ifdef CL_DEVICE_NATIVE_VECTOR_WIDTH_FLOAT
    uival = 0;
    logError( clGetDeviceInfo(*deviceID, CL_DEVICE_NATIVE_VECTOR_WIDTH_FLOAT, sizeof(uival), &uival, NULL) );
    result["CL_DEVICE_NATIVE_VECTOR_WIDTH_FLOAT"] = Rcpp::wrap(uival);
#endif
#ifdef CL_DEVICE_NATIVE_VECTOR_WIDTH_DOUBLE
    uival = 0;
    logError( clGetDeviceInfo(*deviceID, CL_DEVICE_NATIVE_VECTOR_WIDTH_DOUBLE, sizeof(uival), &uival, NULL) );
    result["CL_DEVICE_NATIVE_VECTOR_WIDTH_DOUBLE"] = Rcpp::wrap(uival);
#endif
#ifdef CL_DEVICE_NATIVE_VECTOR_WIDTH_HALF
    uival = 0;
    logError( clGetDeviceInfo(*deviceID, CL_DEVICE_NATIVE_VECTOR_WIDTH_HALF, sizeof(uival), &uival, NULL) );
    result["CL_DEVICE_NATIVE_VECTOR_WIDTH_HALF"] = Rcpp::wrap(uival);
#endif
#ifdef CL_DEVICE_OPENCL_C_VERSION
#endif
#ifdef CL_DEVICE_LINKER_AVAILABLE
#endif
#ifdef CL_DEVICE_BUILT_IN_KERNELS
#endif
#ifdef CL_DEVICE_IMAGE_MAX_BUFFER_SIZE
#endif
#ifdef CL_DEVICE_IMAGE_MAX_ARRAY_SIZE
#endif
#ifdef CL_DEVICE_PARENT_DEVICE
#endif
#ifdef CL_DEVICE_PARTITION_MAX_SUB_DEVICES
#endif
#ifdef CL_DEVICE_PARTITION_PROPERTIES
#endif
#ifdef CL_DEVICE_PARTITION_AFFINITY_DOMAIN
#endif
#ifdef CL_DEVICE_PARTITION_TYPE
#endif
#ifdef CL_DEVICE_REFERENCE_COUNT
#endif
#ifdef CL_DEVICE_PREFERRED_INTEROP_USER_SYNC
#endif
#ifdef CL_DEVICE_PRINTF_BUFFER_SIZE
#endif
#ifdef CL_DEVICE_IMAGE_PITCH_ALIGNMENT
#endif
#ifdef CL_DEVICE_IMAGE_BASE_ADDRESS_ALIGNMENT
#endif
             
    PRINTONERROR();
    return result;
}


SEXP getKernelWorkGroupInfo(SEXP sKernel, SEXP sDeviceID){
    resetError("ROpenCL::getKernelWorkGroupInfo");
    Rcpp::XPtr<cl_kernel> kernel(sKernel);
    Rcpp::XPtr<cl_device_id> deviceID(sDeviceID);
    size_t sval;
    size_t sizes[3];
    cl_ulong ulval;
    Rcpp::List result;
    
#ifdef CL_KERNEL_WORK_GROUP_SIZE /* size_t */
    sval = 0;
    if (!logError( clGetKernelWorkGroupInfo(*kernel, *deviceID, CL_KERNEL_WORK_GROUP_SIZE, sizeof(sval), &sval, NULL) )) {
        result["CL_KERNEL_WORK_GROUP_SIZE"] = Rcpp::wrap(sval);
    }
#endif
#ifdef CL_KERNEL_COMPILE_WORK_GROUP_SIZE /* 3x size_t */
    sizes[0] = sizes[1] = sizes[2] = 0;
    if (!logError( clGetKernelWorkGroupInfo(*kernel, *deviceID, CL_KERNEL_COMPILE_WORK_GROUP_SIZE, 3*sizeof(size_t), &sizes, NULL) )) {
        Rcpp::List wgs;
        wgs.push_back(Rcpp::wrap(sizes[0]));
        wgs.push_back(Rcpp::wrap(sizes[1]));
        wgs.push_back(Rcpp::wrap(sizes[2]));
        result["CL_KERNEL_COMPILE_WORK_GROUP_SIZE"] = Rcpp::wrap(wgs);
    }
#endif
#ifdef CL_KERNEL_LOCAL_MEM_SIZE /* ulong */
    ulval = 0;
    if (!logError( clGetKernelWorkGroupInfo(*kernel, *deviceID, CL_KERNEL_LOCAL_MEM_SIZE, sizeof(ulval), &ulval, NULL) )) {
        result["CL_KERNEL_LOCAL_MEM_SIZE"] = Rcpp::wrap(ulval);
    }
#endif
#ifdef CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE /* size_t */
    sval = 0;
    if (!logError( clGetKernelWorkGroupInfo(*kernel, *deviceID, CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE, sizeof(sval), &sval, NULL) )) {
        result["CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE"] = Rcpp::wrap(sval);
    }
#endif
#ifdef CL_KERNEL_PRIVATE_MEM_SIZE /* ulong */
    ulval = 0;
    if (!logError( clGetKernelWorkGroupInfo(*kernel, *deviceID, CL_KERNEL_PRIVATE_MEM_SIZE, sizeof(ulval), &ulval, NULL) )) {
        result["CL_KERNEL_PRIVATE_MEM_SIZE"] = Rcpp::wrap(ulval);
    }
#endif
//#ifdef CL_KERNEL_GLOBAL_WORK_SIZE /* 3x size_t */ TODO FIX THIS
/*
    sizes[0] = sizes[1] = sizes[2] = 0;
    if (!logError( clGetKernelWorkGroupInfo(*kernel, *deviceID, CL_KERNEL_GLOBAL_WORK_SIZE, 3*sizeof(size_t), &sizes, NULL) )) {
        Rcpp::List gws;
        gws.push_back(Rcpp::wrap(sizes[0]));
        gws.push_back(Rcpp::wrap(sizes[1]));
        gws.push_back(Rcpp::wrap(sizes[2]));
        result["CL_KERNEL_GLOBAL_WORK_SIZE"] = Rcpp::wrap(gws);
    }
*/
//#endif
    
    PRINTONERROR();
    return result;
}


SEXP createContext(SEXP sDeviceID){
    Rcpp::XPtr<cl_device_id> deviceID(sDeviceID);
    
    cl_context *context = new cl_context;
    *context = clCreateContext(0, 1, deviceID, NULL, NULL, NULL);
    Rcpp::XPtr<cl_context> tempXPtr(context);
    return tempXPtr;
}

SEXP createContextFromType(SEXP deviceType){
    //We need to look at this, because this does not seem to work...
    std::string deviceString = Rcpp::as<std::string>(deviceType);
    cl_context *context = new cl_context;
    if(deviceString == "CL_DEVICE_TYPE_GPU"){
        *context = clCreateContextFromType(NULL, CL_DEVICE_TYPE_GPU, NULL, NULL, NULL);
    }
    Rcpp::XPtr<cl_context> tempXPtr(context);
    return tempXPtr;
}

void getContextInfo(SEXP sContext, SEXP sContextInfo){
    //We need to look at this, because this does not seem to work...
    Rcpp::XPtr<cl_context> context(sContext);
    std::string contextInfo = Rcpp::as<std::string>(sContextInfo);
    
    if(contextInfo == "CL_CONTEXT_NUM_DEVICES"){
        cl_uint *nrDevices = new cl_uint;
        std::cout << "get Context number of devices\n";
        cl_int ciErr1 = clGetContextInfo(*context, CL_CONTEXT_NUM_DEVICES, sizeof(nrDevices), nrDevices, NULL);
        if (ciErr1 == CL_INVALID_CONTEXT){std::cout << "invalid context";};
        if (ciErr1 == CL_INVALID_VALUE){std::cout << "invalid value";};
        if (ciErr1 == CL_OUT_OF_RESOURCES){std::cout << "OOR";};
        if (ciErr1 == CL_OUT_OF_HOST_MEMORY){std::cout << "OOHM";};
        std::cout << *nrDevices;
        delete nrDevices;
    }
}

SEXP createCommandQueue(SEXP sContext, SEXP sDeviceID){
    Rcpp::XPtr<cl_context> context(sContext);
    Rcpp::XPtr<cl_device_id> deviceID(sDeviceID);
    
    cl_command_queue *queue = new cl_command_queue;
    *queue = clCreateCommandQueue(*context, *deviceID, 0, NULL);
    Rcpp::XPtr<cl_command_queue> tempXPtr(queue);
    return tempXPtr;
}

static cl_mem_flags getMemFlags(SEXP sMemFlags){
    cl_mem_flags flags = 0L;
    std::string str = Rcpp::as<std::string>(sMemFlags);
    std::string flag;
    // TODO replace any embedded '|' in string with space character
    std::stringstream ss(str); // Insert the string into a stream
    while (ss >> flag) { // Fetch next word from stream (space separated)
#ifdef CL_MEM_READ_WRITE
        if (flag == "CL_MEM_READ_WRITE") flags |= CL_MEM_READ_WRITE;
#endif
#ifdef CL_MEM_WRITE_ONLY
        if (flag == "CL_MEM_WRITE_ONLY") flags |= CL_MEM_WRITE_ONLY;
#endif
#ifdef CL_MEM_READ_ONLY
        if (flag == "CL_MEM_READ_ONLY") flags |= CL_MEM_READ_ONLY;
#endif
#ifdef CL_MEM_USE_HOST_PTR
        if (flag == "CL_MEM_USE_HOST_PTR") flags |= CL_MEM_USE_HOST_PTR;
#endif
#ifdef CL_MEM_ALLOC_HOST_PTR
        if (flag == "CL_MEM_ALLOC_HOST_PTR") flags |= CL_MEM_ALLOC_HOST_PTR;
#endif
#ifdef CL_MEM_COPY_HOST_PTR
        if (flag == "CL_MEM_COPY_HOST_PTR") flags |= CL_MEM_COPY_HOST_PTR;
#endif
#ifdef CL_MEM_HOST_WRITE_ONLY
        if (flag == "CL_MEM_HOST_WRITE_ONLY") flags |= CL_MEM_HOST_WRITE_ONLY;
#endif
#ifdef CL_MEM_HOST_READ_ONLY
        if (flag == "CL_MEM_HOST_READ_ONLY") flags |= CL_MEM_HOST_READ_ONLY;
#endif
#ifdef CL_MEM_HOST_NO_ACCESS
        if (flag == "CL_MEM_HOST_NO_ACCESS") flags |= CL_MEM_HOST_NO_ACCESS;
#endif
    }
    return flags;
}
             
SEXP createBufferFloatVector(SEXP sContext, SEXP sMemFlag, SEXP sGlobalWorkSize){
    Rcpp::XPtr<cl_context> context(sContext);
    int globalWorkSize = Rcpp::as<int>(sGlobalWorkSize);
    
    cl_mem_flags flags = getMemFlags(sMemFlag);
    cl_mem *memBuffer = new cl_mem;
    cl_int ciErr1 = CL_SUCCESS;
    *memBuffer = clCreateBuffer(*context, flags, sizeof(float) * globalWorkSize, NULL, &ciErr1);
    printOnError(ciErr1,"ROpenCL::createBufferFloatVector()");

    Rcpp::XPtr<cl_mem> tempXPtr(memBuffer);
    return tempXPtr;
}

SEXP createBufferIntegerVector(SEXP sContext, SEXP sMemFlag, SEXP sGlobalWorkSize){
    Rcpp::XPtr<cl_context> context(sContext);
    int globalWorkSize = Rcpp::as<int>(sGlobalWorkSize);
    
    cl_mem_flags flags = getMemFlags(sMemFlag);
    cl_mem *memBuffer = new cl_mem;
    cl_int ciErr1 = CL_SUCCESS;
    *memBuffer = clCreateBuffer(*context, flags, sizeof(int) * globalWorkSize, NULL, &ciErr1);
    printOnError(ciErr1,"ROpenCL::createBufferIntegerVector()");
    
    Rcpp::XPtr<cl_mem> tempXPtr(memBuffer);
    return tempXPtr;
}

SEXP createProgramWithSource(SEXP sContext, SEXP sKernelSrc){
    Rcpp::XPtr<cl_context> context(sContext);
    std::string kernelSrc = Rcpp::as<std::string>(sKernelSrc);
    const char* tmpKernelSrc = kernelSrc.data();
    cl_program *program = new cl_program;
    *program = clCreateProgramWithSource(*context, 1, &tmpKernelSrc, NULL, NULL);
    Rcpp::XPtr<cl_program> tempXPtr(program);
    return tempXPtr;
}

SEXP getProgramInfo(SEXP sProgram, SEXP sProgramInfo){
    Rcpp::XPtr<cl_program> program(sProgram);
    std::string programInfo = Rcpp::as<std::string>(sProgramInfo);
    char cBuffer[1024];
    
    if(programInfo == "CL_PROGRAM_SOURCE"){
        std::cout << "get Program source\n";
        clGetProgramInfo(*program, CL_PROGRAM_SOURCE,	sizeof(cBuffer), cBuffer, NULL);
    }
    std::string retVal = cBuffer;
    return Rcpp::wrap(retVal);
}

SEXP buildProgram(SEXP sProgram){
    Rcpp::XPtr<cl_program> program(sProgram);
    cl_int ciErr1 = clBuildProgram(*program, 0, NULL, NULL, NULL, NULL);
    printOnError(ciErr1,"ROpenCL::buildProgram()");
    Rcpp::XPtr<cl_program> tempXPtr(program);
    return tempXPtr;
}

SEXP createKernel(SEXP sProgram, SEXP sKernelName){
    Rcpp::XPtr<cl_program> program(sProgram);
    std::string kernelName = Rcpp::as<std::string>(sKernelName);
    
    cl_kernel *kernel = new cl_kernel;
    *kernel = clCreateKernel(*program, kernelName.data(), NULL);
    Rcpp::XPtr<cl_kernel> tempXPtr(kernel);
    return tempXPtr;
}

SEXP setKernelArgMem(SEXP sKernel, SEXP sIndex, SEXP sBuffer){
    Rcpp::XPtr<cl_kernel> kernel(sKernel);
    int argNr = Rcpp::as<int>(sIndex);
    Rcpp::XPtr<cl_mem> memObject(sBuffer);
    
    cl_int ciErr1 = clSetKernelArg(*kernel, argNr, sizeof(cl_mem), memObject);
    printOnError(ciErr1,"ROpenCL::setKernelArgMem()");
    Rcpp::XPtr<cl_kernel> tempXPtr(kernel);
    return tempXPtr;
}

SEXP setKernelArgInt(SEXP sKernel, SEXP sIndex, SEXP sIntegerValue){
    Rcpp::XPtr<cl_kernel> kernel(sKernel);
    int argNr = Rcpp::as<int>(sIndex);
    int integerValue = Rcpp::as<int>(sIntegerValue);
    
    cl_int ciErr1 = clSetKernelArg(*kernel, argNr, sizeof(cl_int), (void*)&integerValue);
    printOnError(ciErr1,"ROpenCL::setKernelArgInt()");
    Rcpp::XPtr<cl_kernel> tempXPtr(kernel);
    return tempXPtr;
}

SEXP setKernelArgFloat(SEXP sKernel, SEXP sIndex, SEXP sFloatValue){
    Rcpp::XPtr<cl_kernel> kernel(sKernel);
    int argNr = Rcpp::as<int>(sIndex);
    float floatValue = Rcpp::as<float>(sFloatValue);

    cl_int ciErr1 = clSetKernelArg(*kernel, argNr, sizeof(float), (void*)&floatValue);
    printOnError(ciErr1,"ROpenCL::setKernelArgFloat()");
    Rcpp::XPtr<cl_kernel> tempXPtr(kernel);
    return tempXPtr;
}

void enqueueWriteBufferFloatVector(SEXP sQueue, SEXP sMemBuffer, SEXP sGlobalWorkSize, SEXP sObject){
    Rcpp::XPtr<cl_command_queue> queue(sQueue);
    Rcpp::XPtr<cl_mem> clMemBuffer(sMemBuffer);
    int globalWorkSize = Rcpp::as<int>(sGlobalWorkSize);
    Rcpp::NumericVector vec(sObject);
    
    float *object = new float[globalWorkSize];
    for (int i=0; i<vec.size(); i++) {
        object[i] = vec[i];
    }
    cl_int ciErr1 = clEnqueueWriteBuffer(*queue, *clMemBuffer, CL_FALSE, 0, sizeof(float) * globalWorkSize, object, 0, NULL, NULL);
    printOnError(ciErr1,"ROpenCL::enqueueWriteBufferFloatVector()");
}

void enqueueWriteBufferIntegerVector(SEXP sQueue, SEXP sMemBuffer, SEXP sGlobalWorkSize, SEXP sObject){
    Rcpp::XPtr<cl_command_queue> queue(sQueue);
    Rcpp::XPtr<cl_mem> clMemBuffer(sMemBuffer);
    int globalWorkSize = Rcpp::as<int>(sGlobalWorkSize);
    Rcpp::IntegerVector vec(sObject);
    
    int *object = new int[globalWorkSize];
    for (int i=0; i<vec.size(); i++) {
        object[i] = vec[i];
    }
    cl_int ciErr1 = clEnqueueWriteBuffer(*queue, *clMemBuffer, CL_FALSE, 0, sizeof(int) * globalWorkSize, object, 0, NULL, NULL);
    printOnError(ciErr1,"ROpenCL::enqueueWriteBufferIntegerVector()");
}

void enqueue1DRangeKernel(SEXP sQueue, SEXP sKernel, SEXP sGlobalWorkSize, SEXP sLocalWorkSize){
    Rcpp::XPtr<cl_command_queue> queue(sQueue);
    Rcpp::XPtr<cl_kernel> kernel(sKernel);
    size_t globalWorkSize = Rcpp::as<size_t>(sGlobalWorkSize);
    size_t localWorkSize = Rcpp::as<size_t>(sLocalWorkSize);
    size_t *localWorkSizePtr = &localWorkSize;
    if (localWorkSize == (size_t)0) localWorkSizePtr = NULL;
     
    cl_int ciErr1 = clEnqueueNDRangeKernel(*queue, *kernel, 1, NULL, &globalWorkSize, localWorkSizePtr, 0, NULL, NULL);
    printOnError(ciErr1,"ROpenCL::enqueueNDRangeKernel()");
}

void enqueue2DRangeKernel(SEXP sQueue, SEXP sKernel, SEXP sGlobalWorkSize, SEXP sLocalWorkSize){
    Rcpp::XPtr<cl_command_queue> queue(sQueue);
    Rcpp::XPtr<cl_kernel> kernel(sKernel);
     
    Rcpp::NumericVector gws(sGlobalWorkSize);
    Rcpp::NumericVector lws(sLocalWorkSize);
    size_t *globalWorkSizePtr = (size_t *) ::operator new(sizeof(size_t)*2);
    size_t *localWorkSizePtr  = (size_t *) ::operator new(sizeof(size_t)*2);
     
    globalWorkSizePtr[0] = gws[0]; // TODO: check bounds on gws
    globalWorkSizePtr[1] = gws[1];
     
    if (lws[0] == 0 && lws[1] == 0) { // TODO: check bounds on lws
        localWorkSizePtr = NULL;
    } else {
        localWorkSizePtr[0] = lws[0];
        localWorkSizePtr[1] = lws[1];
    }
     
    cl_int ciErr1 = clEnqueueNDRangeKernel(*queue, *kernel, 2, NULL, globalWorkSizePtr, localWorkSizePtr, 0, NULL, NULL);
    printOnError(ciErr1,"ROpenCL::enqueueNDRangeKernel()");
}
    
SEXP enqueueReadBufferFloatVector(SEXP sQueue, SEXP sMemBuffer, SEXP sOffset, SEXP sSize, SEXP sObject){
    Rcpp::XPtr<cl_command_queue> queue(sQueue);
    Rcpp::XPtr<cl_mem> clMemBuffer(sMemBuffer);
    int offset = Rcpp::as<int>(sOffset);
    int size = Rcpp::as<int>(sSize);
    Rcpp::NumericVector vec(sObject);
    
    float *object = new float[size];
    cl_int ciErr1 = clEnqueueReadBuffer(*queue, *clMemBuffer, CL_TRUE, 0, sizeof(float) * size, object, 0, NULL, NULL);
    printOnError(ciErr1,"ROpenCL::enqueueReadBufferFloatVector()");
    for (int i=0; i<size; i++)
    {
        vec[i+offset] = *object++;
    }
    return Rcpp::wrap(vec);
}

SEXP enqueueReadBufferIntegerVector(SEXP sQueue, SEXP sMemBuffer, SEXP sOffset, SEXP sSize, SEXP sObject){
    Rcpp::XPtr<cl_command_queue> queue(sQueue);
    Rcpp::XPtr<cl_mem> clMemBuffer(sMemBuffer);
    int offset = Rcpp::as<int>(sOffset);
    int size = Rcpp::as<int>(sSize);
    Rcpp::IntegerVector vec(sObject);
    
    int *object = new int[size];
    cl_int ciErr1 = clEnqueueReadBuffer(*queue, *clMemBuffer, CL_TRUE, 0, sizeof(int) * size, object, 0, NULL, NULL);
    printOnError(ciErr1,"ROpenCL::enqueueReadBufferIntegerVector()");
    for (int i=0; i<size; i++)
    {
        vec[i+offset] = *object++;
    }
    return Rcpp::wrap(vec);
}


/* // Stuff to expose the int to R
static void intObjFinalizer(SEXP ref){
       if(TYPEOF(ref) == EXTPTRSXP){
               int *o = static_cast<int*> (R_ExternalPtrAddr(ref));
               if (o) delete o;
       }
}

SEXP int2EXP(int *o){
       SEXP xp = R_MakeExternalPtr(o, R_NilValue, R_NilValue);
       R_RegisterCFinalizerEx(xp, intObjFinalizer, TRUE);
       return xp;
}

int *SEXP2int(SEXP o){
       if(TYPEOF(o) != EXTPTRSXP)
               Rf_error("invalid object");
       return (int*) R_ExternalPtrAddr(o);
}

SEXP getIntPointer(){
    int *test = new int;
    *test = 6;
    SEXP retVal = int2EXP(test);
    int test2 = *SEXP2int(retVal);
    return retVal;
}

SEXP doubleIntPointer(SEXP test){
    int test2 = *SEXP2int(test);
    return Rcpp::wrap(test2*2);
}*/
