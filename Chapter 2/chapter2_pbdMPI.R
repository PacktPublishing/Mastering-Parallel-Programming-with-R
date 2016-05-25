# File: chapter2_pbdMPI.R
library(pbdMPI, quietly=TRUE)
init()
pbdmpi_lastsend <- function() {
  myrank <- comm.rank()
  sender <- comm.size() - 1
  receiver <- comm.size() - 2
  if (myrank == sender) {
    msg <- paste("Hi from:", sender)
    send(msg, rank.dest=receiver)
  } else if (myrank == receiver) {
    buf <- recv(rank.source=sender)
  }
  comm.print(buf, rank.print=receiver)
}
pbdmpi_lastsend() # This is SPMD so all processes execute the same
finalize()
