#include <mpi.h>
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include "../../../sprint.h"

int hello(int n, ...)
{
// ignore input args.We don't need them in this example.
  int rank, size, result;

// MPI is initialised and finalised in SPRINT
//  MPI_Init_thread(NULL, NULL); 

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  DEBUG("MPI is initiated in phello rank %d \n", rank);
  Rprintf("Hello from rank %d out of %d\n", rank, size);

//  MPI_Finalize();
  MPI_Barrier(MPI_COMM_WORLD);
  result = 0; // successful execution

  return result;
}