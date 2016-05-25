\name{submitJob}
\alias{submitJob}
\title{Submits a job to a running cluster...}
\usage{submitJob(clusterObject, stopClusterOnComplete=FALSE, taskTimeout=10)}
\description{Submits a job to a running cluster}
\details{After a cluster has been started this function submits jobs to that cluster.
If a job is submitted with enableDebugging=TRUE, all jobs submitted to that
cluster will also have debugging enabled. To turn debugging off, the cluster
must be stopped and restarted.}
\value{Execution status of this job}
\arguments{\item{clusterObject}{a cluster object to submit to}
\item{stopClusterOnComplete}{set to true if you want the cluster to be shut down
after job completes}
\item{taskTimeout}{maximum time a single unit of work can run (in minutes)}}
