\name{emrlapply}
\alias{emrlapply}
\title{Parallel lapply() function using Amazon's EMR service.}
\usage{emrlapply(clusterObject, X, FUN, taskTimeout=10, ...)}
\description{Parallel lapply() function using Amazon's EMR service.}
\details{Parallel lapply() function for applying a function to every item in a list
using Amazon's EMR service.}
\value{Output as a list}
\arguments{\item{X}{list to which the function will be applied}
\item{FUN}{function to apply}
\item{clusterObject}{cluster on which to run the process}
\item{taskTimeout}{maximum time a single unit of work can run (in minutes)}
\item{\dots}{other params to pass to FUN}}
