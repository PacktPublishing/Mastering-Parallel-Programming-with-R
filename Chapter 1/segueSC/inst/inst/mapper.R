#! /usr/bin/env Rscript

trim <- function(line) gsub("(^ +)|( +$)", "", line)


con <- file("stdin", open = "r")

## the libPath is set to include the process ID so that
## multiple instances of R on the same node don't have a
## library locking conflict
pid <- as.character(Sys.getpid())
libPath <- paste("/tmp/R", pid, "/", sep='')


## if you don't want to use the main CRAN site, you should
## change this to a mirror
options(repos=c(CRAN="http://cran.r-project.org/"))
dir.create(libPath)
.libPaths(libPath)

load("./emrData.RData") #contains:
                           # cranPackages - list of packages
                           # myFun - Function to apply
                           # funArgs - the arguments passed
                           # rObjectsOnNodes - a NAMED list of R objects the users wants
                           #                   on each node

attach(rObjectsOnNodes)

install.packages("bitops")
install.packages("caTools")
library(bitops)
library(caTools)

for (myPackage in cranPackages){
  try(install.packages(myPackage) )
  try(library(myPackage,  character.only = TRUE))
  cat("finished installing CRAN packages")
}

getPackageName <- function(sourceFile){
    packPath <- paste("/tmp/Rpackage", pid, "/", sep='')
    untar(sourceFile,  compressed="gzip", exdir=packPath)
    file <- paste(packPath, dir(packPath)[1], "/DESCRIPTION", sep="")
    dcf <- read.dcf(file = file)
    if (NROW(dcf) < 1L) 
          stop(gettextf("DESCRIPTION file of package '%s' is corrupt", 
          pkg), domain = NA)
    packageName <- as.list(dcf[1, ])$Package
    unlink(paste(packPath, dir(packPath)[1], sep=""), recursive=TRUE)
    return(packageName)
}

for (myPackage in customPackages){
  myPackage <- paste("/tmp/segue-upload/", myPackage, sep="")
  try(install.packages(myPackage, repos = NULL, type="source") )
  packageName <- getPackageName(myPackage)
  try(library(packageName,  character.only = TRUE)) 
  cat("finished installing custom packages")
}


## Feb 2010: JDL: moved this down in order to load CRAN packages
##   before the files from the nodes. Objects dependent on
##   CRAN packages would fail if loaded before CRAN packs
## files from filesOnNodes are uploaded to a tmp directory
## this copies the files to the current working directory
try( fileList <- list.files("/tmp/segue-upload/", full.names=TRUE ), silent=TRUE )
try( file.copy(fileList, getwd(), overwrite = TRUE), silent=TRUE)

## try to load the saved workplace image file. This will silently
## fail if there is no workspace file to load
try( load(file="local-workspace-image.RData"), silent=TRUE )

while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
  t <- try( { 
    cat("started readlines \n")
    key <-  as.numeric(trim(strsplit(line, split=",")[[1]][[1]]))
    value <- unserialize(base64decode(strsplit(line, split=",")[[1]][[2]], "raw"))
    value <- list(value)
    value <- c(value, funArgs)
    result <- do.call(myFun, value) # can you believe this one short line does
  } )                                # all the work?!?
  if (inherits(t, "try-error")) result <- paste( "error caught by Segue:", geterrmessage() )
  #serialize and encode the result
  sresult <- paste("<result>,", key, ",", base64encode(serialize(result, NULL, ascii=T)), "\n", sep="")
  cat(sresult, "|\n", sep = "")
}
close(con)
