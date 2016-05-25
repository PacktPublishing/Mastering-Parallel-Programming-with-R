\name{setCredentials}
\alias{setCredentials}
\title{##' AWS Support Function: set up credentials...}
\usage{setCredentials(awsAccessKeyText, awsSecretKeyText,
    setEnvironmentVariables=TRUE)}
\description{##' AWS Support Function: set up credentials}
\details{sets up the credentials needed to access AWS and optionally sets environment
variables for auto loading of credentials in the future}
\author{James "JD" Long}
\arguments{\item{awsAccessKeyText}{your AWS Access Key as a string}
\item{awsSecretKeyText}{your AWS Secret Key as a string}
\item{setEnvironmentVariables}{T/F would you like environment variables to be set so
Segue will read the credentials on load}}
