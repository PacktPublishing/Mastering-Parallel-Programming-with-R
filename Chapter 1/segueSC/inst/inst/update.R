#!/usr/bin/env r
#
# a simple example to update packages in /usr/lib/R/library
# from the littler examples (thanks Dirk!)

## adjust as necessary, see help('download.packages')
repos <- "http://cran.r-project.org"

lib.loc <- "/usr/lib/R/library"

## simply unrolling of all unlink over all files 'repos*' in $TMP
clearCache <- function() {
  sapply(list.files(path=tempdir(), pattern="libloc*", full.names=TRUE), unlink)
}

## Always clear caches of remote and local packages
clearCache()

## r use requires non-interactive use
update.packages(repos=repos, ask=FALSE, lib.loc=lib.loc)

## Always clear caches of remote and local packages
clearCache()
