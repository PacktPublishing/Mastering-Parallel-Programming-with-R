##' Converts an R list into a text CSV.
##' Serializes each element of an R list into ASCII characters then encodes then
##' for use as the input to a Hadoop Streaming job.
##' 
##' @param inList list to convert to CSV
##' @param outFileName file name of resulting CSV file
##' @return creates a CSV to file but returns nothing. 
##' @author James "JD" Long
##' @seealso csvToList
##' 
##' 
##' @export
listToCsv <-
function(inList, outFileName){
  #require(caTools)
  if (is.list(inList) == FALSE) 
        stop("listToCsv: The input list fails the is.list() check.")
  fileName <- outFileName
  cat("", file=fileName, append=FALSE)
  
  i <- 1
  for (item in inList) {
    myLine <- paste(i, ",", base64encode(serialize(item, NULL, ascii=TRUE)), "\n", sep="")
    cat(myLine, file=fileName, append=TRUE) 
    i <- i+1
  }
}

##' The inverse of the listToCsv() function
##' Takes a csv of serialized objects created by the listToCSV function and
##' turns it into a proper R list.
##' 
##' 
##' @param inFileName String pointing to the full path of the input CSV.
##' @return Returns a list object. Or an error. Hopefully a list.
##' @author James "JD" Long
##' @seealso listToCsv()
##' @examples
##'   myList <- NULL
##'   set.seed(1)
##'   for (i in 1:10){
##'     a <- c(rnorm(999), NA)
##'     myList[[i]] <- a
##'   }
##' 
##'   require(caTools)
##'   listToCsv(myList, "tst.csv")
##'   all.equal(myList,  csvToList("tst.csv" ))
##' 
##' @export
csvToList <- function(inFileName){
  #require(caTools)
  linesIn <- readLines(inFileName, n=-1)
  outList <- NULL
  
  i <- 1
  for (line in linesIn){
    outList[[i]] <- unserialize(base64decode(strsplit(linesIn[[i]], split=",")[[1]][[2]], "raw"))
    i <- i+1
  }
  return(outList)
}

