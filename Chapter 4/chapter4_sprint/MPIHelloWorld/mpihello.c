#include <stdio.h>
#include <mpi.h>
int hello(void);

int main(void)
{
  return hello();
}

int hello(void)
{
  int rank, size;

  MPI_Init(NULL,NULL);

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  printf("Hello from rank %d out of %d\n", rank, size);

  MPI_Finalize();

  return 0;
}