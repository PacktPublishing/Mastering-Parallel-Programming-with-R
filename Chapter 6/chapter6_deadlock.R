##
## Copyright 2015 Simon Chapple
##
## Packt: "Mastering Parallelism with R"
## Chapter 6 - Deadlock
##

# run this script in a command shell with mpiexec:
# > mpiexec -np 2 Rscript chapter7_deadlock.R
suppressMessages(library(pbdMPI, quietly = TRUE))
init()
.comm.size <- comm.size()
.comm.rank <- comm.rank()

r <- .comm.rank
succ <- (r + 1) %% .comm.size
pred <- (r - 1) %% .comm.size
comm.print(succ, all.rank = TRUE)
comm.print(pred, all.rank = TRUE)

# vector length of 1000 works (no deadlock) - 10000 does not (deadlock)
v <- 1:10000 # send from this vector
v[1] <- r # set first element to my rank
w <- 1:10000 # receive into this vector

for (i in 1:length(v))
{
  v[i] = r
}

comm.print("sending data (no deadlock)", all.rank = TRUE)
if (r %% 2 == 0) { # even
  send(v,rank.dest = succ)
  w <- recv(w,rank.source = pred)
} else { # odd
  w <- recv(w,rank.source = pred)
  send(v,rank.dest = succ)
}
comm.print(sprintf("%d received message from %d",r,w[1]), all.rank = TRUE)

comm.print("sending data non-blocking", all.rank = TRUE)
isend(v,rank.dest = succ,request=1)
w <- recv(w,rank.source = pred)
wait(request=1)
comm.print(sprintf("%d received message from %d",r,w[1]), all.rank = TRUE)

comm.print("combined sendrecv in one call", all.rank = TRUE)
z <- sendrecv(v,x.buffer = w,rank.dest = succ,rank.source = pred)
comm.print(sprintf("%d received message from %d",r,w[1]), all.rank = TRUE)

comm.print("sending data (may deadlock)", all.rank = TRUE)
send(v,rank.dest = succ)
recv(w,rank.source = pred)
comm.print(sprintf("%d received message from %d",r,w[1]), all.rank = TRUE)

finalize()

##
## Packt: "Mastering Parallelism with R"
## Chapter 7 - Deadlock
##
## Copyright 2015 Simon Chapple
##
