\name{listToCsv}
\alias{listToCsv}
\title{Converts an R list into a text CSV.}
\usage{listToCsv(inList, outFileName)}
\description{Converts an R list into a text CSV.
Serializes each element of an R list into ASCII characters then encodes then
for use as the input to a Hadoop Streaming job.}
\value{creates a CSV to file but returns nothing.}
\author{James "JD" Long}
\seealso{csvToList}
\arguments{\item{inList}{list to convert to CSV}
\item{outFileName}{file name of resulting CSV file}}
