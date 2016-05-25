##
## Copyright 2015 Simon Chapple
##
## Packt: "Mastering Parallelism with R"
## Chapter 5 - OpenCL System Configuration
##

library("inline")
cbody <- 'cl_platform_id pfm[1]; cl_uint np; clGetPlatformIDs(1,pfm,&np);
for (int p = 0; p < np; p++) {/* Outer: Loop over platforms */
    char cb1[128]; char cb2[128]; cl_device_id dev[2];
    cl_uint nd; size_t siz;
    clGetPlatformInfo(pfm[p],CL_PLATFORM_VENDOR,128,cb1,NULL);
    clGetPlatformInfo(pfm[p],CL_PLATFORM_NAME,128,cb2,NULL);
    printf("### Platforms[%d]: %s-%s\\n",p+1,cb1,cb2);
    clGetPlatformInfo(pfm[p],CL_PLATFORM_VERSION,128,cb1,NULL);
    printf("CL_PLATFORM_VERSION: %s\\n",cb1);
    clGetDeviceIDs(pfm[p],CL_DEVICE_TYPE_GPU|CL_DEVICE_TYPE_CPU,2,dev,&nd);
    for (int d = 0; d < nd; d++) {/* Inner: Loop over devices */
        cl_uint uival; cl_ulong ulval; cl_device_type dt;
        size_t szs[10]; cl_device_fp_config fp;
        clGetDeviceInfo(dev[d],CL_DEVICE_VENDOR,128,cb1,NULL);
        clGetDeviceInfo(dev[d],CL_DEVICE_NAME,128,cb2,NULL);
        printf("*** Devices[%d]: %s-%s\\n",d+1,cb1,cb2);
        clGetDeviceInfo(dev[d],CL_DEVICE_TYPE,
                        sizeof(cl_device_type),&dt,NULL);
        printf("CL_DEVICE_TYPE: %s\\n",
               dt & CL_DEVICE_TYPE_GPU ? "GPU" : "CPU");
        clGetDeviceInfo(dev[d],CL_DEVICE_VERSION,128,cb1,NULL);
        printf("CL_DEVICE_VERSION: %s\\n",cb1);
        clGetDeviceInfo(dev[d],CL_DEVICE_MAX_COMPUTE_UNITS,
                        sizeof(cl_uint),&uival,NULL);
        printf("CL_DEVICE_MAX_COMPUTE_UNITS: %u\\n",uival);
        clGetDeviceInfo(dev[d],CL_DEVICE_MAX_CLOCK_FREQUENCY,
                        sizeof(cl_uint),&uival,NULL);
        printf("CL_DEVICE_MAX_CLOCK_FREQUENCY: %u MHz\\n",uival);
        clGetDeviceInfo(dev[d],CL_DEVICE_GLOBAL_MEM_SIZE,
                        sizeof(cl_ulong),&ulval,NULL);
        printf("CL_DEVICE_GLOBAL_MEM_SIZE: %llu Mb\\n",
               ulval/(1024L*1024L));
        clGetDeviceInfo(dev[d],CL_DEVICE_LOCAL_MEM_SIZE,
                        sizeof(cl_ulong),&ulval,NULL);
        printf("CL_DEVICE_LOCAL_MEM_SIZE: %llu Kb\\n",ulval/1024L);
        clGetDeviceInfo(dev[d],CL_DEVICE_DOUBLE_FP_CONFIG,
                        sizeof(cl_device_fp_config),&fp,NULL);
        printf("Supports double precision floating-point? %s\\n",
               fp != 0 ? "yes" : "no");
        }
}'
clfn <- cfunction(signature(), cbody, convention=".C",
        includes=list("#include <stdio.h>","#include <OpenCL/opencl.h>"))
clfn

##
## Packt: "Mastering Parallelism with R"
## Chapter 6 - OpenCL System Configuration
##
## Copyright 2015 Simon Chapple
##
