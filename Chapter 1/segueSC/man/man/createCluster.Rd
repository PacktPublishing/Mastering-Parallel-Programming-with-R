\name{createCluster}
\alias{createCluster}
\title{Creates the configuration object, uploads needed files, and starts...}
\usage{createCluster(numInstances=2, cranPackages, customPackages,
    filesOnNodes, rObjectsOnNodes, enableDebugging=FALSE,
    instancesPerNode, masterInstanceType="m1.large",
    slaveInstanceType="m1.large", location="us-east-1c", ec2KeyName,
    copy.image=FALSE, otherBootstrapActions, sourcePackagesToInstall,
    masterBidPrice, slaveBidPrice)}
\description{Creates the configuration object, uploads needed files, and starts
a Segue Hadoop cluster on Elastic Map Reduce.}
\details{The the needed files are uploaded to S3 and the EMR nodes are started.}
\value{an emrlapply() cluster object with appropriate fields
populated. Keep in mind that this creates the cluster and starts the cluster running.}
\author{James "JD" Long}
\arguments{\item{numInstances}{number of nodes (EC2 instances)}
\item{cranPackages}{vector of string names of CRAN packages to load on each cluster node}
\item{customPackages}{vector of string file names of custom packages to load on each cluster node}
\item{filesOnNodes}{vector of string names of full path of files to be loaded on each node.
Files will be loaded into the local
path (i.e. ./file) on each node.}
\item{rObjectsOnNodes}{a named list of R objects which will be passed to the R
session on the worker nodes. Be sure the list has names. The list will be attached
on the remote nodes using attach(rObjectsOnNodes). If you list does not have names,
this will fail.}
\item{enableDebugging}{T/F whether EMR debugging should be enabled}
\item{instancesPerNode}{Number of R instances per node. Default of NULL uses AWS defaults.}
\item{masterInstanceType}{EC2 instance type for the master node}
\item{slaveInstanceType}{EC2 instance type for the slave nodes}
\item{location}{EC2 location name for the cluster}
\item{ec2KeyName}{EC2 Key used for logging into the main node. Use the user name 'hadoop'}
\item{copy.image}{T/F whether to copy the entire local environment to the nodes. If this feels
fast and loose... you're right! It's nuts. Use it with caution. Very handy when you really need it.}
\item{otherBootstrapActions}{a list-of-lists of other bootstrap actions to run; chlid list members}
\item{sourcePackagesToInstall}{vector of full paths to source packages to be installed on each node}
\item{masterBidPrice}{Bid price for master server}
\item{slaveBidPrice}{Bid price for slave (task) server}}
\examples{\dontrun{
myCluster   <- createCluster(numInstances=2,
cranPackages=c("Hmisc", "plyr"))
}}
