
##' @import rJava
.onLoad <- function(lib, pkg) {
    pathToSdk <- paste(system.file(package = "segueSC") , "/aws-java-sdk/", sep="")

    jarPaths <- c(paste(pathToSdk, "lib/aws-java-sdk-1.10.69.jar", sep=""),
                  paste(pathToSdk, "third-party/lib/commons-logging-1.1.3.jar", sep=""),
                  paste(pathToSdk, "third-party/lib/commons-codec-1.6.jar", sep=""),
                  paste(pathToSdk, "third-party/lib/httpclient-4.3.6.jar", sep=""),
                  paste(pathToSdk, "third-party/lib/httpcore-4.3.3.jar", sep=""),
                  paste(pathToSdk, "third-party/lib/jackson-core-2.5.3.jar", sep=""),
                  paste(pathToSdk, "third-party/lib/jackson-databind-2.5.3.jar", sep=""),
                  paste(pathToSdk, "third-party/lib/jackson-dataformat-cbor-2.5.3.jar", sep=""),
                  paste(pathToSdk, "third-party/lib/jackson-annotations-2.5.0.jar", sep=""),
                  paste(pathToSdk, "third-party/lib/joda-time-2.8.1.jar", sep=""),
                  paste(pathToSdk, "third-party/", sep="")
                  )
    .jpackage(pkg, morePaths=jarPaths)
    attach( javaImport( c("java.lang", "java.io")))
    
    if (Sys.getenv("AWSACCESSKEY") != "" && Sys.getenv("AWSSECRETKEY") != ""){
      awsCreds <- new(com.amazonaws.auth.BasicAWSCredentials, Sys.getenv("AWSACCESSKEY"), Sys.getenv("AWSSECRETKEY"))
      assign("awsCreds", awsCreds, envir = .GlobalEnv)
      packageStartupMessage( "segueSC has loaded your AWS Credentials." )
    } else {
       packageStartupMessage( "segueSC did not find your AWS credentials. Please run the setCredentials() function." )
    }
}
