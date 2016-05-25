##' Parallel lapply() function using Amazon's EMR service.
##'
##' Parallel lapply() function for applying a function to every item in a list
##' using Amazon's EMR service.
##' 
##' 
##' @param X list to which the function will be applied
##' @param FUN function to apply
##' @param clusterObject cluster on which to run the process
##' @param taskTimeout maximum time a single unit of work can run (in minutes)
##' @param \dots other params to pass to FUN
##' @return Output as a list
##' 
##' @export
emrlapply <- function(clusterObject, X, FUN, taskTimeout=10, ... ) {
  #set up a local temp directory
  myTempDir <- clusterObject$localTempDir

  #the function to apply gets put into myFunction
  myFun <- FUN
  funArgs <-  as.list(substitute(list(...)))[-1L]

  ## create a string vector of the customPackages short names for use in loading in the mapper
  customPackages <- NULL
  for (name in clusterObject$customPackages) {
    name <- basename(name)
    customPackages <- c(customPackages, name)
  }
  
  cranPackages <- clusterObject$cranPackages 
  rObjectsOnNodes <- clusterObject$rObjectsOnNodes
  
  #save the objects
  objectsFileName <-paste(myTempDir ,"/emrData.RData", sep="") 
  save(cranPackages,
       myFun,
       funArgs,
       rObjectsOnNodes,
       customPackages,
       file = objectsFileName,
       compress=TRUE)

  #make sure the bucket exists, and is empty
  #try(makeS3Bucket(clusterObject$s3TempDir)
  emptyS3Bucket(clusterObject$s3TempDir)

  #the out director must NOT exist
  try(deleteS3Bucket(clusterObject$s3TempDirOut), silent=TRUE)
  
  #upload the datafile to S3
  uploadS3File(clusterObject$s3TempDir, paste(objectsFileName, sep=""))
    
  #upload the mapper to S3
  uploadS3File(clusterObject$s3TempDir, system.file("mapper.R", package="segueSC"))

  #serialize the X list to a temp file
  streamFile <- paste(myTempDir, "/stream.txt", sep="")
  listToCsv(X, streamFile)
  
  #now upload stream.txt to EMR
  uploadS3File(clusterObject$s3TempDir, streamFile)
  
  finalStatus <- submitJob(clusterObject=clusterObject, taskTimeout=taskTimeout) 
  myTempDirOut <- clusterObject$localTempDirOut

  
  #if (finalStatus %in% c("COMPLETED", "WAITING")) {
  #  system(paste("mkdir ", myTempDirOut, sep="" ))
  #  system(paste("rm ", myTempDirOut, "/*", sep=""))

  downloadS3File(clusterObject$s3TempDirOut, ".all", myTempDirOut)
  ##  the results are going to be in the subdirectory "results"
  myTempDirOut <- paste(myTempDirOut, "results", sep="/")

    #open files
  returnedFiles <- list.files(path=myTempDirOut, pattern="part")
    #yes, I read all the results into R then write them out to a text file
    #There was a reason for doing this, but I don't remember it
    #this could all be done in one step
  combinedOutputFile <- file(paste(myTempDirOut, "/combinedOutput.csv", sep=""), "w")
  unparsedOutput <- NULL
    for (file in returnedFiles){
        lines <- readLines(paste(myTempDirOut, "/", file, sep="")) 
        for (line in lines) {
          if (substr(line, 1, 9) == "<result>,") {
            write(substr(line, 10, nchar(line)), file=combinedOutputFile)
          }
        }
    }
    close(combinedOutputFile)
    
    #require(caTools)
    lines <- strsplit(readLines(paste(myTempDirOut, "/combinedOutput.csv", sep="")),
                      split=",")
    output <- list()
    
    for (i in 1:length(lines)){
      output[[as.numeric(lines[[i]][[1]])]] <- (unserialize(
                                                 base64decode(
                                                   substr(
                                                     lines[[i]][[2]],
                                                          1,
                                                          nchar(lines[[i]][[2]])-1), "raw")))
    }
    return(as.list(output))
}

##' Internal function used for assembling the output produced by the EMR process
##'
##' Internal function used for assembling the output produced by the EMR process
##' 
##' @param myTempDirOut the temp directory where the output files from the EMR
##'                     process are kept. Must be local, not S3.
##' @return Returns a list of output that matches the input list in length
##' @author James "JD" Long
assembleOutput <- function(myTempDirOut) {

    #open files
  returnedFiles <- list.files(path=myTempDirOut, pattern="part")
    #yes, I read all the results into R then write them out to a text file
    #There was a reason for doing this, but I don't remember it
    #this could all be done in one step
  combinedOutputFile <- file(paste(myTempDirOut, "/combinedOutput.csv", sep=""), "w")
  unparsedOutput <- NULL
    for (file in returnedFiles){
        lines <- readLines(paste(myTempDirOut, "/", file, sep="")) 
        for (line in lines) {
          if (substr(line, 1, 9) == "<result>,") {
            write(substr(line, 10, nchar(line)), file=combinedOutputFile)
          }
        }
    }
    close(combinedOutputFile)
    
    #require(caTools)
    lines <- strsplit(readLines(paste(myTempDirOut, "/combinedOutput.csv", sep="")),
                      split=",")
    output <- list()
    
    for (i in 1:length(lines)){
      output[[as.numeric(lines[[i]][[1]])]] <- (unserialize(
                                                 base64decode(
                                                   substr(
                                                     lines[[i]][[2]],
                                                          1,
                                                          nchar(lines[[i]][[2]])-1), "raw")))
    }
    return(as.list(output))
}



