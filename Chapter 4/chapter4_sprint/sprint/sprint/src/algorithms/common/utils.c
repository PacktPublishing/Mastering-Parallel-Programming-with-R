/**************************************************************************
 *                                                                        *
 *  SPRINT: Simple Parallel R INTerface                                   *
 *  Copyright © 2010 The University of Edinburgh                          *
 *                                                                        *
 *  This program is free software: you can redistribute it and/or modify  *
 *  it under the terms of the GNU General Public License as published by  *
 *  the Free Software Foundation, either version 3 of the License, or     *
 *  any later version.                                                    *
 *                                                                        *
 *  This program is distributed in the hope that it will be useful,       *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the          *
 *  GNU General Public License for more details.                          *
 *                                                                        *
 *  You should have received a copy of the GNU General Public License     *
 *  along with this program. If not, see <http://www.gnu.or/licenses/>.   *
 *                                                                        *
 **************************************************************************/
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Parse.h>

#include "../../sprint.h"
#include "utils.h"

int getMaximum(int *nmax, double *ammax, int world_size) {

  int maximum_indx = -1;
  double maximum = 0;

  for(int i=0; i<world_size; i++) {

    if (ammax[i] >= maximum) {
      maximum = ammax[i];
      maximum_indx = nmax[i];
    }
  }
  return maximum_indx; 
}

int getMinimumIndx(double *array, int array_length) {
  
  int minimum_indx = 0;
  double minimum = array[0];

  for(int i=1; i<array_length; i++) {

    if (array[i] <= minimum) {
      minimum = array[i];
      minimum_indx = i;
    }
  }
  return minimum_indx; 
}

void getMedoidIDs(int *nrepr, int *medoid_ids, int n_rows, int n_clusters) {

  int k = 0;

  for (int i=0; i<n_rows; i++) {
    if (nrepr[i] == 1) {
      medoid_ids[k] = i;
      k++;
      if (k == n_clusters)
        return;
    }
  }
}

// TODO: Replace with MPI_Gatherv, add a type parameter and implement using
// switch statement

void mMPI_AllgatherINT(int* send_buf, int world_size, int my_rank,
                    int my_start, int my_end, int n_rows, MPI_Comm comm) {

  int my_buf_size, recbuf_indx, i;
  int *buf_size;
  MPI_Status status;
    
  buf_size = Calloc(world_size, int);
  my_buf_size = my_end - my_start;
  recbuf_indx = 0;

  MPI_Gather(&my_buf_size, 1, MPI_INT, buf_size, 1, MPI_INT, MASTER_PROCESS, comm);

  if (my_rank == MASTER_PROCESS) {
    
    recbuf_indx = buf_size[MASTER_PROCESS];
    
    for (i=1; i<world_size; i++) {
     
      MPI_Recv(&send_buf[recbuf_indx], buf_size[i], MPI_INT, i, 0, comm, &status);
      recbuf_indx += buf_size[i];
    }
    
  } else {
   
    MPI_Send(&send_buf[my_start], my_buf_size, MPI_INT, MASTER_PROCESS, 0, comm);
  }
    
  MPI_Bcast(send_buf, n_rows, MPI_INT, 0, MPI_COMM_WORLD);
 
  Free(buf_size);
}

void mMPI_AllgatherDOUBLE(double* send_buf, int world_size, int my_rank,
                    int my_start, int my_end, int n_rows, MPI_Comm comm) {

  int my_buf_size, recbuf_indx, i;
  int *buf_size;
  MPI_Status status;

  buf_size = Calloc(world_size, int);
  my_buf_size = my_end - my_start;
  recbuf_indx = 0;

  MPI_Gather(&my_buf_size, 1, MPI_INT, buf_size, 1, MPI_INT, MASTER_PROCESS, comm);

  if (my_rank == MASTER_PROCESS) {
    
    recbuf_indx = buf_size[MASTER_PROCESS];
    
    for (i=1; i<world_size; i++) {
      MPI_Recv(&send_buf[recbuf_indx], buf_size[i], MPI_DOUBLE, i, 0, comm, &status);
      recbuf_indx += buf_size[i];
    }
    
  } else {
    MPI_Send(&send_buf[my_start], my_buf_size, MPI_DOUBLE, MASTER_PROCESS, 0, comm);
  }
    
  MPI_Bcast(send_buf, n_rows, MPI_DOUBLE, 0, MPI_COMM_WORLD);
 
  Free(buf_size);
}


/*
  Search for the largest distance 
*/
double getMaxDistance(double *distance_matrix, int my_start, int my_end, int nn) {

  int i, j;
  double maxD, maxD_reduced;
  
  /* s := max( dys[.] ), the largest distance */
  for (i = my_start, maxD = 0.; i < my_end; i++) {
    for(j = i; j < nn; j++) {
      
      //read only upper triangular of the distance matrix
      if (maxD < distance_matrix[DIST_INDEX(i,j,nn)]) {
        maxD = distance_matrix[DIST_INDEX(i,j,nn)];
      }
    }
  }

  MPI_Allreduce(&maxD, &maxD_reduced, 1, MPI_DOUBLE, MPI_MAX,MPI_COMM_WORLD);
  
  return maxD_reduced;
}

/*
  Medoids Initialisation

*/
void initMedoids(int *nrepr, int *med, Rboolean med_given, int nn, int n_clusters) {

  for (int i = 0; i < nn; ++i)
    nrepr[i] = 0;
  if(med_given) { /* if true, med[] contain initial medoids */

	/* for the moment, translate these to nrepr[] 0/1 :
	 * not assuming that the med[] indices are sorted */
	for (int k = 0; k < n_clusters; k++)
    nrepr[med[k] - 1] = 1;
  }
  return;
}

  

/*
  Distributes array (matrix) among processes
*/ 

void loopDistribute(int myid, int num_of_proc, int N,
                    int *my_start, int *my_end)
{

  int tail;

  if (num_of_proc == 1) {

    *my_start = 0;
    *my_end = N;
    
  } else {

    tail = N % (num_of_proc);

    *my_start = myid * (N/num_of_proc);
    *my_end   = *my_start + (N/num_of_proc);

    if(myid < tail && tail > 0) {
      if(myid == 0) {
        *my_end = *my_end + 1;
      } else {
        *my_end = *my_end + myid + 1;
        *my_start = *my_start + myid;
      }
    } else {
      *my_start = *my_start + tail;
      *my_end = *my_end + tail;
    }
  }
}

/* R_ParseVector is essentially the code used to implement parse(text=) at R level. */
SEXP parseExpression(SEXP expressionSexp) {

  SEXP  parsedCmd = R_NilValue;
  ParseStatus status; 
   
  parsedCmd = PROTECT(R_ParseVector(expressionSexp, -1, &status, R_NilValue));
	if (status != PARSE_OK) {
    UNPROTECT(1);
    error("invalid expression in parseExpression!");
  }

  UNPROTECT(1);
  return parsedCmd;
}

