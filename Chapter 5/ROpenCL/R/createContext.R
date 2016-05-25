#' Function to create a context
#' @return pointer to the context
#' @author Willem Ligtenberg
#' @export
createContext <- function(){
	.Call( "createContext", PACKAGE = "ROpenCL" )
}

#' Create a buffer
#' @param sContext context
#' @param sMemFlag memflag
#' @param sGlobalWorkSize globalworksize
#' @param sObject object
#' @return pointer to the buffer
#' @author Willem Ligtenberg
#' @export
createBuffer <- function(sContext, sMemFlag, sGlobalWorkSize, sObject){
    type <- class(sObject)
    if (type == "integer") {
        buffer <- createBufferIntegerVector(sContext, sMemFlag, sGlobalWorkSize)
    }
    if (type == "numeric") {
        buffer <- createBufferFloatVector(sContext, sMemFlag, sGlobalWorkSize)
    }
    if (!(type %in% c("integer", "numeric"))){
      stop(paste("Objects of class ", type, " are not supported yet", sep = ""))
    }
    return(buffer)
}

#' Build a kernel
#' @param context context in which to execute the kernel
#' @param kernelSrc source code of the kernel function as a string
#' @param kernelName name of the kernel
#' @param ... further arguments
#' @return pointer to the kernel
#' @author Willem Ligtenberg
#' @export
buildKernel <- function(context, kernelSrc, kernelName, ...){
    program <- createProgramWithSource(context, kernelSrc)
    program <- buildProgram(program)
    kernel <- createKernel(program, kernelName)
    
    dotList <- list(...)
    index <- 0
    for (item in dotList){
        if(class(item) == "externalptr"){
            kernel <- setKernelArgMem(kernel, index, item)
        }
        if(class(item) == "integer"){
            kernel <- setKernelArgInt(kernel, index, item)
        }
        if(class(item) == "numeric"){
          logdebug('Adding a Float argument to the kernel')
          kernel <- setKernelArgFloat(kernel, index, item)
        }
        index <- index + 1
    }
    return(kernel)
}

#' Change the value of a kernel function parameter
#' @param kernel previously created kernel object
#' @param index 1 for first parameter, 2 for second parameter (R based 1-indexing)
#' @param item new value for the parameter (could be a cl mem buffer, integer or numeric)
#' @return kernel updated kernel object for re-enqueueing
#' @author Simon Chapple
#' @export
assignKernelArg <- function(kernel,index,item){
    index <- index - 1
    if(class(item) == "externalptr"){
        kernel <- setKernelArgMem(kernel, index, item)
    }
    else if(class(item) == "integer"){
        kernel <- setKernelArgInt(kernel, index, item)
    }
    else if(class(item) == "numeric"){
        kernel <- setKernelArgFloat(kernel, index, item)
    }
    return(kernel)
}

#' Enqueue a write buffer
#' @param sQueue queue
#' @param sMemBuffer membufffer
#' @param size (typically) GlobalWorkSize
#' @param sObject object
#' @return pointer
#' @author Willem Ligtenberg
#' @export
enqueueWriteBuffer <- function(sQueue, sMemBuffer, size, sObject){
    type <- class(sObject)
    if(type == "integer") {
        enqueueWriteBufferIntegerVector(sQueue, sMemBuffer, size, sObject)
    }
    if(type == "numeric") {
        enqueueWriteBufferFloatVector(sQueue, sMemBuffer, size, sObject)
    }
    if(!(type %in% c("integer", "numeric"))){
      stop(paste("Objects of class ", type, " are not supported yet", sep = ""))
    }
}

#' Enqueue the read buffer
#' @param sQueue queue
#' @param sMemBuffer membuffer
#' @param size (typically) GlobalWorkSize
#' @param sObject object
#' @param offset (optional) index into vector for copy
#' @return pointer
#' @author Willem Ligtenberg & Simon Chapple
#' @export
enqueueReadBuffer <- function(sQueue, sMemBuffer, size, sObject, offset=0){
    type <- class(sObject)
    if(type == "integer") {
        buffer <- enqueueReadBufferIntegerVector(sQueue, sMemBuffer, as.integer(offset), as.integer(size), sObject)
    }
    if(type == "numeric") {
        buffer <- enqueueReadBufferFloatVector(sQueue, sMemBuffer, as.integer(offset), as.integer(size), sObject)
    }
    if(!(type %in% c("integer", "numeric"))){
        stop(paste("Objects of class ", type, " are not supported yet", sep = ""))
    }
    return(buffer)
}


#' Enqueue the kernel
#' @param queue queue
#' @param kernel kernel
#' @param globalWorkSize global worksize
#' @param localWorkSize local worksize (optional)
#' @return ?
#' @author Simon Chapple
#' @export
enqueueNDRangeKernel <- function(queue, kernel, globalWorkSize, localWorkSize = NULL){
    len <- length(globalWorkSize)
    if (len == 1) {
        if (is.null(localWorkSize)) { localWorkSize <- 0 }
        enqueue1DRangeKernel(queue, kernel, globalWorkSize, localWorkSize)
    }
    if (len == 2) {
        if (is.null(localWorkSize)) { localWorkSize <- c(0,0) }
        enqueue2DRangeKernel(queue, kernel, globalWorkSize, localWorkSize)
    }
}

#' Determine if the device supports double precision
#' @param deviceID device ID
#' @param listName list to which will be appended the cl floating-point configuration - if its defined (optional)
#' @return boolean
#' @author Simon Chapple
#' @export
deviceSupportsDoublePrecision <- function(deviceID, listName = NULL){
    list <- getDeviceDoubleFP(deviceID)
    if (length(list) > 0) {
        if (!is.null(listName)) {
            cfg <- get(listName, envir=parent.frame())
            for (item in list) {
                cfg[length(cfg)+1] <- item
            }
            assign(listName, cfg, envir=parent.frame())
        }
        return(TRUE)
    }
    return(FALSE)
}

#' Determine if the device supports single precision
#' @param deviceID device ID
#' @param listName list to which will be appended the cl floating-point configuration - if its defined (optional)
#' @return boolean
#' @author Simon Chapple
#' @export
deviceSupportsSinglePrecision <- function(deviceID, listName = NULL){
    list <- getDeviceSingleFP(deviceID)
    if (length(list) > 0) {
        if (!is.null(listName)) {
            cfg <- get(listName, envir=parent.frame())
            for (item in list) {
                cfg[length(cfg)+1] <- item
            }
            assign(listName, cfg, envir=parent.frame())
        }
        return(TRUE)
    }
    return(FALSE)
}

#' Determine if the device supports half precision
#' @param deviceID device ID
#' @param listName list to which will be appended the cl floating-point configuration - if its defined (optional)
#' @return boolean
#' @author Simon Chapple
#' @export
deviceSupportsHalfPrecision <- function(deviceID, listName = NULL){
    list <- getDeviceHalfFP(deviceID)
    if (length(list) > 0) {
        if (!is.null(listName)) {
            cfg <- get(listName, envir=parent.frame())
            for (item in list) {
                cfg[length(cfg)+1] <- item
            }
            assign(listName, cfg, envir=parent.frame())
        }
        return(TRUE)
    }
    return(FALSE)
}


