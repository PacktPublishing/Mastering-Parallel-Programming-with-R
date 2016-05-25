
## A little simplification would be the first step toward rational living, I think.
## Eleanor Roosevelt 


## Lower logging level: LogManager.getLogManager().getLogger("com.amazonaws.request").setLevel(Level.OFF);
## ref: https://forums.aws.amazon.com/thread.jspa?messageID=186655&#186655

##' ##' AWS Support Function: set up credentials
##'
##' sets up the credentials needed to access AWS and optionally sets environment
##' variables for auto loading of credentials in the future
##' @param awsAccessKeyText your AWS Access Key as a string
##' @param awsSecretKeyText your AWS Secret Key as a string
##' @param setEnvironmentVariables T/F would you like environment variables to be set so
##' Segue will read the credentials on load
##' @author James "JD" Long
##' @export
setCredentials <- function(awsAccessKeyText, awsSecretKeyText, setEnvironmentVariables = TRUE){
    awsCreds <- new(com.amazonaws.auth.BasicAWSCredentials, awsAccessKeyText, awsSecretKeyText)
    assign("awsCreds", awsCreds, envir = .GlobalEnv)

    if (setEnvironmentVariables == TRUE) {
      Sys.setenv(AWSACCESSKEY = awsAccessKeyText, AWSSECRETKEY = awsSecretKeyText)
    }
}
##' AWS Support Function: Delete an S3 Key (a.k.a file)
##'
##' Deteles a key in a given bucket on S3
##' @param bucketName name of the bucket
##' @param keyName the key in the bucket
##' @author James "JD" Long
##' @export
deleteS3Key <- function(bucketName, keyName){
  tx <- new(com.amazonaws.services.s3.transfer.TransferManager, awsCreds)
  s3 <- tx$getAmazonS3Client()
  if (s3$doesBucketExist(bucketName)) { 
    s3$deleteObject(bucketName, keyName)
  }
}

##' AWS Support Function: Empty an S3 bucket
##'
##' Deletes all keys in the designated bucket
##' @param bucketName Name of the bucket to be emptied
##' @author James "JD" Long
##' @export
emptyS3Bucket <- function(bucketName){
  tx <- new(com.amazonaws.services.s3.transfer.TransferManager, awsCreds)
  s3 <- tx$getAmazonS3Client()

  # TODO: need a check to make sure the current user owns the bucket
  #       before trying to delete everything in it
  #       there's some risk this might loop forever if they don't own the bucket
  if (s3$doesBucketExist(bucketName)) {  
    lst <- s3$listObjects(bucketName)
    objSums <- lst$getObjectSummaries()
    listJavaObjs <- .jevalArray(objSums$toArray())
    if (length(listJavaObjs)>0){
      for (i in 1:length(listJavaObjs)) {
        deleteS3Key(bucketName, listJavaObjs[[i]]$getKey()[[1]])
      }
    }
    if (lst$isTruncated()){
      #recursion FTW!
      emptyS3Bucket(bucketName)
    }
  }
}

##' AWS Support Function: Delete an S3 Bucket
##'
##' Does nothing if the bucketName does not exist. If bucket contains Keys,
##' all keys are deleted.
##' @param bucketName the bucket to be deleted
##' @author James "JD" Long
##' @export
deleteS3Bucket <- function(bucketName){
  tx <- new(com.amazonaws.services.s3.transfer.TransferManager, awsCreds)
  s3 <- tx$getAmazonS3Client()
  if (s3$doesBucketExist(bucketName) == TRUE) {
    emptyS3Bucket(bucketName)
    tx <- new(com.amazonaws.services.s3.transfer.TransferManager, awsCreds)
    s3 <- tx$getAmazonS3Client()
    s3$deleteBucket(bucketName)
  }
}

##' AWS Support Function: Creates an S3 Bucket
##'
##' Creates an S3 bucket. If the bucket already exists, no warning is returned.
##' @param bucketName string of the name of the bucket to be created
##' @author James "JD" Long
##' @export
makeS3Bucket <- function(bucketName){
    #awsCreds <- get("awsCreds", envir = segue.env)
    tx <- new(com.amazonaws.services.s3.transfer.TransferManager, awsCreds)
    s3 <- tx$getAmazonS3Client()
    #test if the bucket exists; if not,  make bucket
    if (s3$doesBucketExist(bucketName) == FALSE) {
      s3$createBucket(bucketName)
    } else {
      warning("Unable to Create Bucket. Bucket with same name already exists.", call. = FALSE)
    }
}

##' AWS Support Function: Uploads a local file to an S3 Bucket
##'
##' If buckName does not exist, it is created and a warning is issued. 
##' @param bucketName destination bucket
##' @param localFile local file to be uploaded
##' @author James "JD" Long
##' @export
uploadS3File <- function(bucketName, localFile){
    tx <- new(com.amazonaws.services.s3.transfer.TransferManager, awsCreds)
    s3 <- tx$getAmazonS3Client()
    fileToUpload <-  new(File, localFile)
    request <- new(com.amazonaws.services.s3.model.PutObjectRequest, bucketName, fileToUpload$getName(), fileToUpload)
    s3$putObject(request)
}

##' AWS Support Function: Downloads a key from an S3 Bucket into a local file.
##'
##' Pulls a key (file) from a bucket into a localFile. If the keyName = ".all" then
##' all files from the bucket are pulled and localFile should be a
##' directory name. Ignores "sub directories" in buckets. 
##' @param bucketName destination bucket
##' @param keyName key to download. ".all" to pull all keys
##' @param localFile local file name or path if ".all" is called for keyName
##' @author James "JD" Long
##' @export
downloadS3File <- function(bucketName, keyName, localFile){
    tx <- new(com.amazonaws.services.s3.transfer.TransferManager, awsCreds)
    s3 <- tx$getAmazonS3Client()
    if (keyName != ".all") {
      request <- new(com.amazonaws.services.s3.model.GetObjectRequest, bucketName, keyName)
      theObject <- s3$getObject(request, new(java.io.File, localFile))
    } else {
     # this will only pull the first page of listings
     # so if there are a lot of files it won't grab them all
     # 
     # TODO: make it pull multiple pages of files
     # TODO: pull subdirectories too
      system(paste("mkdir", localFile), ignore.stderr = TRUE)
      lst <- s3$listObjects(bucketName)
      objSums <- lst$getObjectSummaries()
      listJavaObjs <- .jevalArray(objSums$toArray())
      if (length(listJavaObjs)>0){
        for (i in 1:length(listJavaObjs)) {
          # if statement here just to filter out subdirs
          key <- listJavaObjs[[i]]$getKey()[[1]]
          #if ( length( unlist(strsplit(key, split="/")) ) == 1) {
            if (substring( key, nchar( key ) - 7, nchar( key ) )  != "$folder$") {
              localFullFile <- paste(localFile, "/", listJavaObjs[[i]]$getKey()[[1]], sep="")
              downloadS3File(bucketName, listJavaObjs[[i]]$getKey()[[1]], localFullFile)
            }
          #}
        }
      }
    }
  }
  
##' Creates the configuration object, uploads needed files, and starts
##' a Segue Hadoop cluster on Elastic Map Reduce. 
##'
##' The the needed files are uploaded to S3 and the EMR nodes are started.
##' @param numInstances number of nodes (EC2 instances)
##' @param cranPackages vector of string names of CRAN packages to load on each cluster node
##' @param customPackages vector of string file names of custom packages to load on each cluster node
##' @param filesOnNodes vector of string names of full path of files to be loaded on each node.
##' Files will be loaded into the local
##' path (i.e. ./file) on each node. 
##' @param rObjectsOnNodes a named list of R objects which will be passed to the R
##' session on the worker nodes. Be sure the list has names. The list will be attached
##' on the remote nodes using attach(rObjectsOnNodes). If you list does not have names,
##' this will fail.
##' @param enableDebugging T/F whether EMR debugging should be enabled
##' @param instancesPerNode Number of R instances per node. Default of NULL uses AWS defaults.
##' @param masterInstanceType EC2 instance type for the master node
##' @param slaveInstanceType EC2 instance type for the slave nodes
##' @param location EC2 location name for the cluster
##' @param ec2KeyName EC2 Key used for logging into the main node. Use the user name 'hadoop'
##' @param copy.image T/F whether to copy the entire local environment to the nodes. If this feels
##' fast and loose... you're right! It's nuts. Use it with caution. Very handy when you really need it.
##' @param otherBootstrapActions a list-of-lists of other bootstrap actions to run; chlid list members
##    are: "name" == unique identifier of this bootstrap action ; "localFile" == path to local script
##    to be uploaded to the temp area in S3; "s3file" == path to an existing script in S3 (won't be
##    uploaded to the temp area); "args" == vector of character arguments.   "localFile" and "s3file"
##    are mutually exclusive but one is required; "args" is optional.
##' @param sourcePackagesToInstall vector of full paths to source packages to be installed on each node
##' @param masterBidPrice Bid price for master server
##' @param slaveBidPrice Bid price for slave (task) server
##' @return an emrlapply() cluster object with appropriate fields
##'   populated. Keep in mind that this creates the cluster and starts the cluster running.
##' @author James "JD" Long
##' @examples
##' \dontrun{
##' myCluster   <- createCluster(numInstances=2,
##'  cranPackages=c("Hmisc", "plyr"))
##' }
##' @export
createCluster <- function(numInstances=2,
                          cranPackages=NULL,
                          customPackages=NULL,
                          filesOnNodes=NULL,
                          rObjectsOnNodes=NULL, 
                          enableDebugging=FALSE,
                          instancesPerNode=NULL,
                          masterInstanceType="m1.large",
                          slaveInstanceType="m1.large",
                          location = "us-east-1c",
                          ec2KeyName=NULL,
                          copy.image=FALSE ,
                          otherBootstrapActions=NULL,
                          sourcePackagesToInstall=NULL,
                          masterBidPrice=NULL,
                          slaveBidPrice=NULL
                          ){
  ## this used to be an argument but not bootstrapping
  ## caused too many problems
  bootStrapLatestR=TRUE

  clusterObject <- list(numInstances = numInstances,
                        cranPackages = cranPackages,
                        customPackages = customPackages, 
                        enableDebugging = enableDebugging,
                        bootStrapLatestR = bootStrapLatestR,
                        filesOnNodes = filesOnNodes,
                        rObjectsOnNodes = rObjectsOnNodes, 
                        enableDebugging = enableDebugging,
                        instancesPerNode = instancesPerNode,
                        masterInstanceType = masterInstanceType,
                        slaveInstanceType = slaveInstanceType,
                        location = location,
                        ec2KeyName = ec2KeyName ,
                        copy.image = copy.image ,
                        otherBootstrapActions = otherBootstrapActions,
                        sourcePackagesToInstall = sourcePackagesToInstall,
                        masterBidPrice = masterBidPrice,
                        slaveBidPrice = slaveBidPrice
                        )

  if ( tolower(masterInstanceType) == "m1.small") {
    clusterObject$masterInstanceType <- "m1.large"
    print("WARNING: masterInstanceType set to m1.small. Segue requires 64 bit OS so the masterInstanceType is being changed to m1.large. You will be billed by Amazon accordingly.")
  }
  if ( tolower(slaveInstanceType) == "m1.small") {
    clusterObject$slaveInstanceType <- "m1.large"
    print("WARNING: slaveInstanceType set to m1.small. Segue requires 64 bit OS so the slaveInstanceType is being changed to m1.large. You will be billed by Amazon accordingly.")
  }


  
  localTempDir <- paste(tempdir(),
                        paste(sample(c(0:9, letters), 10, rep=T), collapse=""),
                        "-segue",
                        sep="")
  
  clusterObject$localTempDir <- localTempDir
  clusterObject$localTempDirOut <- paste(localTempDir, "/out", sep="")

  dir.create(localTempDir, showWarnings = TRUE, recursive = TRUE, mode = "0777")
  dir.create(clusterObject$localTempDirOut, showWarnings = TRUE, recursive = TRUE, mode = "0777")
  
  s3TempDir <- tolower(unlist(strsplit(localTempDir, "/"))[length(unlist(strsplit(localTempDir, "/")))])
  deleteS3Bucket(s3TempDir)
  clusterObject$s3TempDir <- s3TempDir
  
  s3TempDirOut <- tolower(paste(s3TempDir , "out", sep=""))
  deleteS3Bucket(s3TempDirOut)
  clusterObject$s3TempDirOut <- s3TempDirOut

  #create the s3 bucket
  ## TODO: error check this
  makeS3Bucket(s3TempDir)
  print(paste0("segueSC: s3 temp = s3://",s3TempDir))
  
  s3TempDirLogs <- paste(s3TempDir, "-logs", sep="")
  makeS3Bucket(s3TempDirLogs)
  print(paste0("segueSC: s3 temp logs = s3://",s3TempDirLogs))
  
  #upload the bootstrapper to S3 
  if (bootStrapLatestR==TRUE) {
    ##TODO: error checking in the uploadS3File function
    uploadS3File(s3TempDir, system.file("bootstrapR.sh", package="segueSC") )
    #uploadS3File(s3TempDir, system.file("update.R", package="segueSC") )
    
  }
  clusterObject$bootStrapLatestR <- bootStrapLatestR

  ## if copy.image is TRUE then save an image and  use the fileOnNodes
  ## feature to add the saved image to the nodes
  if (copy.image == TRUE) {
    imageFile <- paste( localTempDir, "/local-workspace-image.RData", sep="" )
    save.image( file=imageFile, compress=TRUE )
    clusterObject$filesOnNodes = c(clusterObject$filesOnNodes, imageFile)
  }
  
  ## if customPackages are present, add them to the filesOnNodes
  if (is.null(customPackages) == FALSE) {
    clusterObject$filesOnNodes = c(clusterObject$filesOnNodes, customPackages)
  }
  
  # start cluster
  jobFlowId <- startCluster(clusterObject)
  clusterObject$jobFlowId <- jobFlowId
  return(clusterObject)
}

##' AWS Support Function: Checks the status of a given job on EMR
##'
##' Checks the status of a previously issued job.
##' @param jobFlowId the Job Flow Id of the job to check
##' @return Job Status 
##' @author Simon R Chapple
##' @export
checkStatus <- function(jobFlowId){
  service <- new( com.amazonaws.services.elasticmapreduce.AmazonElasticMapReduceClient, awsCreds )
  request <- new( com.amazonaws.services.elasticmapreduce.model.DescribeClusterRequest )
  request$setClusterId(jobFlowId)
  result <- service$describeCluster(request)
  result$getCluster()$getStatus()$getState()
}

##' AWS Support Function: Checks the status of a given job on EMR
##'
##' Checks the status of a previously issued step.
##' @param jobFlowId the Job Flow Id of the job to check
##' @return Status of the last step 
##' @author James "JD" Long
##' @export
checkLastStepStatus <- function(jobFlowId){
  service <- new( com.amazonaws.services.elasticmapreduce.AmazonElasticMapReduceClient, awsCreds )
  request <- new( com.amazonaws.services.elasticmapreduce.model.DescribeJobFlowsRequest )
  detailsList <- new( java.util.ArrayList )
  detailsList$add(jobFlowId)
  request$setJobFlowIds(detailsList)
  descriptions <- as.list(service$describeJobFlows(request)$getJobFlows())
  #descriptions[[1]]$getExecutionStatusDetail()$getState()

  steps <- as.list(descriptions[[1]]$getSteps())
  step <- steps[[length(steps)]] #grab the last step only
  status <- step$getExecutionStatusDetail()
  status$getState()
}


##' Starts a cluster on Amazon's EMR service
##' 
##' After a cluster has been defined with createCluster() this function actually
##' starts the machines running. Currently exported, but soon will be internal only.
##' 
##' @param clusterObject cluster object to start
##' @return a Job Flow ID
##' 
##' @export
startCluster <- function(clusterObject){
  numInstances     <- clusterObject$numInstances
  s3TempDir        <- clusterObject$s3TempDir 
  s3TempDirOut     <- clusterObject$s3TempDirOut
  bootStrapLatestR <- clusterObject$bootStrapLatestR
  verbose          <- TRUE
 
  service <- new( com.amazonaws.services.elasticmapreduce.AmazonElasticMapReduceClient, awsCreds )
  request <- new( com.amazonaws.services.elasticmapreduce.model.RunJobFlowRequest )
  conf    <- new( com.amazonaws.services.elasticmapreduce.model.JobFlowInstancesConfig )

  #creates the bootstrap list
  bootStrapList <- new( java.util.ArrayList )
  
  if (bootStrapLatestR == TRUE) {
   scriptBootActionConfig <- new(com.amazonaws.services.elasticmapreduce.model.ScriptBootstrapActionConfig)
   scriptBootActionConfig$setPath(paste("s3://", s3TempDir, "/bootstrapR.sh", sep=""))
  
   bootStrapConfig <- new( com.amazonaws.services.elasticmapreduce.model.BootstrapActionConfig)
     with( bootStrapConfig, setScriptBootstrapAction(scriptBootActionConfig))
     with( bootStrapConfig, setName("R-InstallLatest"))
 
   bootStrapList$add(bootStrapConfig)

   ## update packages
   #   scriptBootActionConfig <- new(com.amazonaws.services.elasticmapreduce.model.ScriptBootstrapActionConfig)
   #scriptBootActionConfig$setPath(paste("s3://", s3TempDir, "/update.R", sep=""))
  
  #bootStrapConfig <- new( com.amazonaws.services.elasticmapreduce.model.BootstrapActionConfig)
  #with( bootStrapConfig, setScriptBootstrapAction(scriptBootActionConfig))
  #with( bootStrapConfig, setName("R-UpdatePackages"))
 
 #bootStrapList$add(bootStrapConfig)
   
  }

  ## handle additional bootstrap actions, if requested.
  if ( ! is.null(clusterObject$otherBootstrapActions) ){
    ## TODO: more graceful exit here? or would stopifnot() be appropriate, in this case?
    stopifnot( "list" == class(clusterObject$otherBootstrapActions) )

    invisible( sapply( clusterObject$otherBootstrapActions , function( action ){
      scriptBootActionConfig <- new(com.amazonaws.services.elasticmapreduce.model.ScriptBootstrapActionConfig)

      ## are we uploading a local file to run? or will we use a script that already exists in
      ## (a non-temporary) S3 bucket?
      if( ! is.null( action$localFile ) ){
        uploadS3File(clusterObject$s3TempDir , action$localFile)
        scriptBootActionConfig$setPath(paste("s3://", clusterObject$s3TempDir, "/" , basename( action$localFile ) , sep=""))
      }else if( ! is.null( action$s3file ) ){
        scriptBootActionConfig$setPath(action$s3file)
      }

      if( ! is.null( action$args ) ){
        ## TODO: proper quoting around args? or leave that for caller?
        argsAsList <- new( java.util.ArrayList )
        sapply( action$args , function(item){ argsAsList$add(item) } )
        scriptBootActionConfig$withArgs(argsAsList)
      }

  
      bootStrapConfig <- new( com.amazonaws.services.elasticmapreduce.model.BootstrapActionConfig)
        with( bootStrapConfig, setScriptBootstrapAction(scriptBootActionConfig))
        with( bootStrapConfig, setName(action$name))

      bootStrapList$add(bootStrapConfig)
    
    } ) )
  }

  if (is.null(clusterObject$filesOnNodes) == FALSE) { # putting files on each node
    print("INFO: You have selected files to be put on each node. These files are being uploaded to S3.")
    ## build a batch file that includes each element of filesOnNodes
    ## then add the batch file as a boot strap action

    ## open the output file (bootStrapFiles.sh) in clusterObject$tempDir
    ## open an output file connection
    outfile <- file( paste( clusterObject$localTempDir, "/bootStrapFiles.sh", sep="" ), "w" )  
    cat("#!/bin/bash", "", file = outfile, sep = "\n")
    cat("mkdir /tmp/segue-upload/", "", file = outfile, sep = "\n")
     ## for each element in filesOnNodes add a hadoop -fs line
    for ( file in clusterObject$filesOnNodes ){
      remotePath <- paste( "/tmp/segue-upload/", tail(strsplit(file,"/")[[1]], 1), sep="" )
      fileName <- tail(strsplit(file,"/")[[1]], 1)
      s3Path <- paste( "s3://", clusterObject$s3TempDir, "/", fileName, sep="" )
      cat( paste( "hadoop fs -get ", s3Path, remotePath)
          , file = outfile, sep = "\n" )
      cat( "\n", file = outfile )
      
      # copy each file to S3
      uploadS3File( clusterObject$s3TempDir, file )
    }
    close( outfile )
     # copy bootStrapFiles.sh to clusterObject$s3TempDir
    uploadS3File( clusterObject$s3TempDir, paste( clusterObject$localTempDir, "/bootStrapFiles.sh", sep="" ) )

     # add a bootstrap action that runs bootStrapFiles.sh
   scriptBootActionConfig <- new(com.amazonaws.services.elasticmapreduce.model.ScriptBootstrapActionConfig)
   scriptBootActionConfig$setPath(paste("s3://", s3TempDir, "/bootStrapFiles.sh", sep=""))
  
   bootStrapConfig <- new( com.amazonaws.services.elasticmapreduce.model.BootstrapActionConfig)
     with( bootStrapConfig, setScriptBootstrapAction(scriptBootActionConfig))
     with( bootStrapConfig, setName("RBootStrapFiles"))
 
   bootStrapList$add(bootStrapConfig)
   print("INFO: Upload of files to S3 is complete.")
  }

  
  if (is.null(clusterObject$sourcePackagesToInstall) == FALSE) {
    print("INFO: Now building sources packages to install and uploading them based on the sourcePackagesToInstall list.")
    ## build a batch file that includes each file in sourcePackagesToInstall
    ## then add the batch file as a boot strap action

    ## open the output file (installSourcePackages.sh) in clusterObject$tempDir
    ## open an output file connection
    outfile <- file( paste( clusterObject$localTempDir, "/installSourcePackages.sh", sep="" ), "w" )  
    cat("#!/bin/bash", "", file = outfile, sep = "\n")
    cat("mkdir /tmp/segue-source-packages/", "", file = outfile, sep = "\n")
     ## for each element in sourcePackagesToInstall add a hadoop -fs line
    for ( file in clusterObject$sourcePackagesToInstall ){
      remotePath <- paste( "/tmp/segue-source-packages/", tail(strsplit(file,"/")[[1]], 1), sep="" )
      fileName <- tail(strsplit(file,"/")[[1]], 1)
      s3Path <- paste( "s3://", clusterObject$s3TempDir, "/", fileName, sep="" )
      cat( paste( "hadoop fs -get ", s3Path, remotePath)
          , file = outfile, sep = "\n" )
      cat( "\n", file = outfile )
      cat( "sudo R CMD INSTALL ", remotePath, "\n", file = outfile, sep = "" )
      
      # copy each file to S3
      uploadS3File( clusterObject$s3TempDir, file )
    }
    close( outfile )
     # copy installSourcePackages.sh) to clusterObject$s3TempDir
    uploadS3File( clusterObject$s3TempDir, paste( clusterObject$localTempDir, "/installSourcePackages.sh", sep="" ) )

     # add a bootstrap action that runs bootStrapFiles.sh
   scriptBootActionConfig <- new(com.amazonaws.services.elasticmapreduce.model.ScriptBootstrapActionConfig)
   scriptBootActionConfig$setPath(paste("s3://", s3TempDir, "/installSourcePackages.sh", sep=""))
  
   bootStrapConfig <- new( com.amazonaws.services.elasticmapreduce.model.BootstrapActionConfig)
     with( bootStrapConfig, setScriptBootstrapAction(scriptBootActionConfig))
     with( bootStrapConfig, setName("RinstallSourcePackages"))
 
   bootStrapList$add(bootStrapConfig)
   print("INFO: Source packages uploaded.")
  }

  

  if (is.null(clusterObject$instancesPerNode) == FALSE) { #sersiously... test this
   scriptBootActionConfig <- new(com.amazonaws.services.elasticmapreduce.model.ScriptBootstrapActionConfig)
   scriptBootActionConfig$setPath("s3://elasticmapreduce/bootstrap-actions/configure-hadoop")

   argList <- new( java.util.ArrayList )
   argList$add( "-s" )
   argList$add( paste( "mapred.tasktracker.map.tasks.maximum=", clusterObject$instancesPerNode, sep="") )
   argList$add( "-s" )
   argList$add( paste( "mapred.tasktracker.reduce.tasks.maximum=", clusterObject$instancesPerNode, sep="") )
 
   scriptBootActionConfig$setArgs( argList )
                                  
   bootStrapConfig <- new( com.amazonaws.services.elasticmapreduce.model.BootstrapActionConfig)
     with( bootStrapConfig, setScriptBootstrapAction(scriptBootActionConfig))
     with( bootStrapConfig, setName("SetInstancePerNode"))
 
   bootStrapList$add(bootStrapConfig)    
  }
          
   ## this adds the bootstrap to the request
  
   request$setBootstrapActions(bootStrapList)
  
  if ( is.null( clusterObject$ec2KeyName ) != TRUE ) {
      conf$setEc2KeyName(clusterObject$ec2KeyName)
  }
  
  
  instanceGroups <- .jnew("java/util/Vector")
  if (!is.null(clusterObject$masterBidPrice)) {
      masterGroupConf <- new( com.amazonaws.services.elasticmapreduce.model.InstanceGroupConfig )
      masterGroupConf$setInstanceCount(new(Integer, as.character(1)))
      masterGroupConf$setInstanceRole(new(String, as.character("MASTER"))) 
      masterGroupConf$setInstanceType(new(String, as.character(clusterObject$masterInstanceType))) 
      masterGroupConf$setMarket(new(String, as.character("SPOT"))) 
      masterGroupConf$setBidPrice(new(String, as.character(clusterObject$masterBidPrice))) 
      instanceGroups$add(masterGroupConf)
  }
  if (!is.null(clusterObject$slaveBidPrice)) {
      slaveGroupConf <- new( com.amazonaws.services.elasticmapreduce.model.InstanceGroupConfig )
      slaveGroupConf$setInstanceCount(new(Integer, as.character(numInstances - 1)))
      slaveGroupConf$setInstanceRole(new(String, as.character("CORE"))) 
      slaveGroupConf$setInstanceType(new(String, as.character(clusterObject$slaveInstanceType))) 
      slaveGroupConf$setMarket(new(String, as.character("SPOT"))) 
      slaveGroupConf$setBidPrice(new(String, as.character(clusterObject$slaveBidPrice))) 
      instanceGroups$add(slaveGroupConf)
  }
  if (!is.null(clusterObject$masterBidPrice) || !is.null(clusterObject$slaveBidPrice)) {
      conf$setInstanceGroups(instanceGroups)
  } else {
      # Must configure instances either using instance count, 
      # master and slave instance type or instance groups but not both
      conf$setInstanceCount(new(Integer, as.character(numInstances)))
      conf$setMasterInstanceType( clusterObject$masterInstanceType )      
      conf$setSlaveInstanceType( clusterObject$slaveInstanceType )
  }
  conf$setKeepJobFlowAliveWhenNoSteps(new(Boolean, TRUE))

  conf$setPlacement(new(com.amazonaws.services.elasticmapreduce.model.PlacementType, clusterObject$location))
  conf$setHadoopVersion("2.4.0") ## conf$setHadoopVersion("0.20.205")
  request$setInstances(conf)
  request$setLogUri(paste("s3://", s3TempDir, "-logs", sep=""))
  jobFlowName <- paste("RJob-", date(), sep="")
  request$setName(jobFlowName)
  request$setAmiVersion("3.11") ## request$setAmiVersion("2.0.4")
  request$setJobFlowRole("EMR_EC2_DefaultRole")
  request$setServiceRole("EMR_DefaultRole")

  result <- service$runJobFlow(request)

  ## seems like this sleep should not be needed... but otherwise
  ## getJobFlowId() does not get the correct jobflowid

  Sys.sleep(15)
  jobFlowId <- result$getJobFlowId()
  print(paste0("segueSC: Cluster JobFlowID=",jobFlowId))
  
  currentStatus <- checkStatus(jobFlowId)
  #currentStatus <- "WAITING"
  
  while (currentStatus  %in% c("COMPLETED", "FAILED", "TERMINATED", "WAITING", "CANCELLED")  == FALSE) {
    Sys.sleep(30)
    currentStatus <- checkStatus(jobFlowId)
    message(paste(currentStatus, " - ", Sys.time(), sep="" ))
  }
 
  if (currentStatus == "WAITING") {
    message("Your Amazon EMR Hadoop Cluster is ready for action. \nRemember to terminate your cluster with stopCluster().\nAmazon is billing you!")
  }
  return(jobFlowId)
  
  ## TODO: need to catch situations where the cluster failed
}

##' Stops a running cluster
##'
##' Stops a running cluster and deletes temporary directories from EC2
##' 
##' @return nothing
##' @author James "JD" Long
##' @param clusterObject a cluster object to stop
##' @export
stopCluster <- function(clusterObject){
  jobFlowId <- clusterObject$jobFlowId

  service <- new( com.amazonaws.services.elasticmapreduce.AmazonElasticMapReduceClient, awsCreds )
  request <- new( com.amazonaws.services.elasticmapreduce.model.TerminateJobFlowsRequest )
  detailsList <- new( java.util.ArrayList )
  detailsList$add(jobFlowId)
  request$withJobFlowIds(detailsList)
  service$terminateJobFlows(request)

  ## I have no idea why AWS needs sleep before
  ## I can delete the temp dirs, but these fail
  ## if I don't have the sleep
  Sys.sleep(10)
  try( deleteS3Bucket( clusterObject$s3TempDir ), silent=TRUE )
  try( deleteS3Bucket( paste( clusterObject$s3TempDir, "-logs", sep="" ) ), silent=TRUE )
  try( deleteS3Bucket( clusterObject$s3TempDirOut ), silent=TRUE )

  ## something weird is going on... I have to do this twice or it
  ## does not fully delete the s3TempDir's subdirectory
  ## will need to give this some attention later
  Sys.sleep(15)
  try( deleteS3Bucket( clusterObject$s3TempDir ), silent=TRUE )
  try( deleteS3Bucket( paste( clusterObject$s3TempDir, "-logs", sep="" ) ), silent=TRUE )
  try( deleteS3Bucket( clusterObject$s3TempDirOut ), silent=TRUE )
 
}

##' Submits a job to a running cluster
##' 
##' After a cluster has been started this function submits jobs to that cluster.
##' If a job is submitted with enableDebugging=TRUE, all jobs submitted to that
##' cluster will also have debugging enabled. To turn debugging off, the cluster
##' must be stopped and restarted.
##' 
##' 
##' @param clusterObject a cluster object to submit to
##' @param stopClusterOnComplete set to true if you want the cluster to be shut down
##' after job completes
##' @param taskTimeout maximum time a single unit of work can run (in minutes)
##' @return Execution status of this job
##' 
##' @export
submitJob <- function(clusterObject, stopClusterOnComplete=FALSE, taskTimeout=10){

  jobFlowId       <- clusterObject$jobFlowId
  s3TempDir       <- clusterObject$s3TempDir
  s3TempDirOut    <- clusterObject$s3TempDirOut
  enableDebugging <- clusterObject$enableDebugging
 
  try(deleteS3Bucket(s3TempDirOut), silent=TRUE)
  unlink(clusterObject$localTempDirOut, recursive = TRUE)
  dir.create(clusterObject$localTempDirOut)
  
  jobFlowId <- clusterObject$jobFlowId

  if (enableDebugging==TRUE){
    service <- new( com.amazonaws.services.elasticmapreduce.AmazonElasticMapReduceClient, awsCreds )
    hadoopJarStep <- new(com.amazonaws.services.elasticmapreduce.model.HadoopJarStepConfig)
    hadoopJarStep$setJar("s3://us-east-1.elasticmapreduce/libs/script-runner/script-runner.jar")
    argList <- new( java.util.ArrayList )
    argList$add( "s3://us-east-1.elasticmapreduce/libs/state-pusher/0.1/fetch" )
    hadoopJarStep$setArgs(argList)
    stepName <- format(Sys.time(), "%Y-%m-%d_%H:%M:%OS5") 
    stepConfig <- new(com.amazonaws.services.elasticmapreduce.model.StepConfig, stepName, hadoopJarStep)
    stepConfig$setActionOnFailure("CANCEL_AND_WAIT")
    stepList <- new( java.util.ArrayList )
    stepList$add( stepConfig )
    request <- new( com.amazonaws.services.elasticmapreduce.model.AddJobFlowStepsRequest, jobFlowId, stepList)
    requestResult <- service$addJobFlowSteps(request)
  }
  
  service <- new( com.amazonaws.services.elasticmapreduce.AmazonElasticMapReduceClient, awsCreds )

  hadoopJarStep <- new(com.amazonaws.services.elasticmapreduce.model.HadoopJarStepConfig)
  hadoopJarStep$setJar("/home/hadoop/contrib/streaming/hadoop-streaming.jar")
  argList <- new( java.util.ArrayList )

  # the task timeout is passed to us in minutes, but AWS/EMR expects it in milliseconds
  taskTimeoutMilliseconds <- taskTimeout * 60 * 1000
  argList$add( "-D" )
  argList$add( paste( "mapreduce.task.timeout=" , as.integer(taskTimeoutMilliseconds) , sep="" ) )

  # s3n requires explicit access, let's try s3...
  argList$add( "-files" )
  argList$add( paste("s3://", s3TempDir, "/emrData.RData",  #"/emrData.RData#emrData.RData",
                     ",s3://", s3TempDir, "/mapper.R", sep=""))  #"/mapper.R#mapper.R"
  argList$add( "-input" )
  argList$add( paste("s3n://", s3TempDir, "/stream.txt", sep="") )
  argList$add( "-output" )
  argList$add( paste("s3n://", s3TempDirOut, "/results", sep="") )
  argList$add( "-mapper" )
  argList$add( "cat" )
  argList$add( "-reducer" )
  argList$add( "mapper.R" )


  hadoopJarStep$setArgs(argList)

  stepName <- format(Sys.time(), "%Y-%m-%d_%H:%M:%OS5") 
  
  stepConfig <- new(com.amazonaws.services.elasticmapreduce.model.StepConfig, stepName, hadoopJarStep)
  stepConfig$setActionOnFailure("CANCEL_AND_WAIT")
  
  stepList <- new( java.util.ArrayList )
  stepList$add( stepConfig )
  request <- new( com.amazonaws.services.elasticmapreduce.model.AddJobFlowStepsRequest, jobFlowId, stepList)

  try(deleteS3Bucket(clusterObject$s3TempDirOut), silent=TRUE)
  
  #submit to EMR happens here
  print("SegueSC: Submitting Job")
  service$addJobFlowSteps(request)
  print("SegueSC: Job Submitted")
  
  Sys.sleep(15)
  
  print(paste0("SegueSC: cluster status = ",checkStatus(jobFlowId)))

  Sys.sleep(15)

  currentStatus <- checkStatus(jobFlowId)
  #currentStatus <- "COMPLETED"
  
  while (currentStatus  %in% c("COMPLETED", "FAILED", "TERMINATED", "WAITING", "CANCELLED")  == FALSE) {
    Sys.sleep(30)
    currentStatus <- checkStatus(jobFlowId)
    message(paste(currentStatus, " - ", Sys.time(), sep="" ))
  }
  if (stopClusterOnComplete==TRUE) {
    stopCluster(clusterObject)
  }
  return(currentStatus)
  
}

