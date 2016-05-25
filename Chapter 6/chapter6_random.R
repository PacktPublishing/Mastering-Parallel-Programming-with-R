##
## Copyright 2015 Simon Chapple
##
## Packt: "Mastering Parallelism with R"
## Chapter 6 - Random
##

# run this script in a command shell with mpiexec:
# > mpiexec -np 2 Rscript chapter7_random.R
suppressMessages(library(pbdMPI, quietly = TRUE))
init()

# Create unique seed for each MPI process
seed <- as.numeric(format(Sys.time(),"%OS6")) * as.numeric(Sys.getpid())
comm.print(seed, all.rank = TRUE)

finalize()

##
## Packt: "Mastering Parallelism with R"
## Chapter 7 - Random
##
## Copyright 2015 Simon Chapple
##
