\name{downloadS3File}
\alias{downloadS3File}
\title{AWS Support Function: Downloads a key from an S3 Bucket into a local file.}
\usage{downloadS3File(bucketName, keyName, localFile)}
\description{AWS Support Function: Downloads a key from an S3 Bucket into a local file.}
\details{Pulls a key (file) from a bucket into a localFile. If the keyName = ".all" then
all files from the bucket are pulled and localFile should be a
directory name. Ignores "sub directories" in buckets.}
\author{James "JD" Long}
\arguments{\item{bucketName}{destination bucket}
\item{keyName}{key to download. ".all" to pull all keys}
\item{localFile}{local file name or path if ".all" is called for keyName}}
