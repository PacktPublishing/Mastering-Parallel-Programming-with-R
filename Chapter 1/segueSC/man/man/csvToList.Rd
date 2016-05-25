\name{csvToList}
\alias{csvToList}
\title{The inverse of the listToCsv() function...}
\usage{csvToList(inFileName)}
\description{The inverse of the listToCsv() function
Takes a csv of serialized objects created by the listToCSV function and
turns it into a proper R list.}
\value{Returns a list object. Or an error. Hopefully a list.}
\author{James "JD" Long}
\seealso{listToCsv()}
\arguments{\item{inFileName}{String pointing to the full path of the input CSV.}}
\examples{myList <- NULL
set.seed(1)
for (i in 1:10){
a <- c(rnorm(999), NA)
myList[[i]] <- a
}

require(caTools)
listToCsv(myList, "tst.csv")
all.equal(myList,  csvToList("tst.csv" ))}
