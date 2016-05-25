#include <mpi.h>
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
SEXP hello(void);

SEXP hello(void)
{
  int rank, size;

  MPI_Init(NULL, NULL);

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  Rprintf("Hello from rank %d out of %d\n", rank, size);

  MPI_Finalize();

  SEXP result =PROTECT(result = NEW_INTEGER(1));
  INTEGER(result)[0] = 0;
  UNPROTECT(1);
  return result;
}
