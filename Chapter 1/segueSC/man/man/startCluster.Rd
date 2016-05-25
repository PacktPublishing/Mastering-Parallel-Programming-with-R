\name{startCluster}
\alias{startCluster}
\title{Starts a cluster on Amazon's EMR service...}
\usage{startCluster(clusterObject)}
\description{Starts a cluster on Amazon's EMR service}
\details{After a cluster has been defined with createCluster() this function actually
starts the machines running. Currently exported, but soon will be internal only.}
\value{a Job Flow ID}
\arguments{\item{clusterObject}{cluster object to start}}
