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
 *                             2                                          *
 *  You should have received a copy of the GNU General Public License     *
 *  along with this program. If not, see <http://www.gnu.or/licenses/>.   *
 *                                                                        *
 **************************************************************************/

#include <R.h>
#include <Rinternals.h>

#include "../../../sprint.h"
#include "../../common/utils.h"
#include "apply.h"
#include "comms.h"

#define NAME_LENGTH 256

SEXP matrixApply(SEXP result, SEXP data, SEXP margin, SEXP function,
                 int worldRank, int worldSize) {

  SEXP ans, data_size;
  
  MPI_Datatype row_type, column_type;
  MPI_Status status;
  
  int my_start, my_end, N, function_nlines,
    nvectors, offset;
  int local_check = 0, global_check = 0;
  int dimensions[2];

  if (worldRank == MASTER_PROCESS) { 
    data_size = GET_DIM(data);
    
    dimensions[0] = INTEGER_POINTER(data_size)[0];
    dimensions[1] = INTEGER_POINTER(data_size)[1];

    /* function SEXP object is a vector of strings, each element contains
       a single line of the function definition */

    function_nlines = length(function);   
  }

  MPI_Bcast(dimensions, 2, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&function_nlines, 1, MPI_INT, 0, MPI_COMM_WORLD);

  /* margin provides the subscripts which the function will be
     applied over.  "1" indicates rows, "2" indicates columns,
     c(1,2)" indicates rows and columns */

  if(worldRank != MASTER_PROCESS)
    PROTECT(margin = allocVector(INTSXP, 1));

  MPI_Bcast(INTEGER(margin), 1, MPI_INT, 0, MPI_COMM_WORLD);

  /* Matrix dimensions in R are interpreted differen than in C.
     We will refer to R rows and columns ordering, so rows are not alligned
     in memory */

  if (INTEGER(margin)[0] == 1) {
    N = dimensions[0];

    /* define vector type type to handle R rows exchange
       (count, blocklength, stride)*/
    MPI_Type_vector (dimensions[1], 1, dimensions[0], MPI_DOUBLE, &row_type);
    MPI_Type_commit (&row_type);

  } else if (INTEGER(margin)[0] == 2) {
    N = dimensions[1];

    /* define contiguous type to handle R columns exchange */
    MPI_Type_contiguous(dimensions[0], MPI_DOUBLE, &column_type);
    MPI_Type_commit(&column_type);
    
  } else if (INTEGER(margin)[0] == 3) {
    // TODO
    DEBUG("Margin number 3 not yet implemented\n");
    return R_NilValue;
  } else {
    DEBUG("Don't know how to deal with margin number %d\n",
          INTEGER(margin)[0]);
    return R_NilValue;
  }
  
  if(worldRank != MASTER_PROCESS) {
  
    /* Allocate memory for SEXP objects on worker nodes.
       alloc... functions do their own error-checking and
       return if the allocation process will fail. */
    loopDistribute(worldRank, worldSize, N, &my_start, &my_end);
    
    if (INTEGER(margin)[0] == 1)
      PROTECT(data = allocMatrix(REALSXP, my_end-my_start, dimensions[1]));
    if (INTEGER(margin)[0] == 2)
      PROTECT(data = allocMatrix(REALSXP, dimensions[0], my_end-my_start));
              
    PROTECT(function = allocVector(STRSXP, function_nlines));
  }

  if ( (data == NULL) ||  (function == NULL) ) {
      local_check = 1;
  } else {
      local_check = 0;
  }

  /* Check whether memory was successfully allocated on all worker nodes */
  MPI_Allreduce(&local_check, &global_check, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  /*  No need to free memory if allocation fails on one of the workers
      R_alloc will release it after .Call returns to R */
  if ( global_check != 0 ) {
    /* Remove all references from the stack, I'm not sure if this is necessary */
    if(worldRank != MASTER_PROCESS)
      UNPROTECT(3);

    return ScalarInteger(-1);
  }

  /* Distribute data between processes */

  for (int worker_id=1; worker_id<worldSize; worker_id++) {

    if (worldRank == MASTER_PROCESS) {

      /* Calculate expected message length for each worker */
      loopDistribute(worker_id, worldSize, N, &my_start, &my_end);
      nvectors = my_end - my_start;

      /* If we applying over rows, as defined in R, we need to use the MPI vector type
         sending each row as a separate message */
      if (INTEGER(margin)[0] == 1) {
        for(int k=0; k<nvectors; k++) {
          offset = my_start+k;
          MPI_Send(&REAL(data)[offset], 1, row_type, worker_id, 0, MPI_COMM_WORLD);
        }
      }

      /* R defined columns are alligned in memory, single message of build from contiguous
         column_type elemensts is send */
      else if (INTEGER(margin)[0] == 2) {
        offset = my_start*dimensions[0];
        MPI_Send(&REAL(data)[offset], nvectors, column_type, worker_id, 0, MPI_COMM_WORLD);
      }
    }
    else if (worldRank == worker_id) {

      nvectors = my_end - my_start;

      if (INTEGER(margin)[0] == 1) {
        
        for(int k=0; k<nvectors; k++) {
          offset = k*dimensions[1];
          MPI_Recv(&REAL(data)[offset], dimensions[1], MPI_DOUBLE, MASTER_PROCESS, 0, MPI_COMM_WORLD, &status);
        }
      }
      else if (INTEGER(margin)[0] == 2) {
        MPI_Recv(REAL(data), nvectors, column_type, MASTER_PROCESS, 0, MPI_COMM_WORLD, &status);
      }
    }
  }

  /* Redo loop distribution for the Master process */
  if (worldRank == MASTER_PROCESS) {
      loopDistribute(worldRank, worldSize, N, &my_start, &my_end);
  }
    
  /* Bcast function name or definition, cover case when definition is split into
     several lines and stored as a SEXP string vector */
   bcastRFunction(function, function_nlines, worldRank);
  
  /* Response container, Vector of SEXPs, margin determines vector length */
  PROTECT(ans = allocVector(VECSXP, N));

  do_matrixApply(ans, data, margin, function, my_start, my_end, dimensions, worldRank);

  gatherData(result, ans, N, my_start, my_end, worldRank);
  
  if(worldRank != MASTER_PROCESS) {
    UNPROTECT(4);
  } else {
    UNPROTECT(1);
  }

  return result;

}

void do_matrixApply(SEXP ans,
                    SEXP data,
                    SEXP margin,
                    SEXP function,
                    int my_start,
                    int my_end,
                    int *dimensions,
                    int worldRank)
{

  SEXP data_chunk, R_fcall, parsedCmd = R_NilValue;
  int i, j;
  double *rdata, *rchunk;

  rdata = REAL(data);

  PROTECT(parsedCmd = parseExpression(function));

  /* Create R LANGSXP Vector, R function holder
     length of the vector is 1 + number of arguments */
  PROTECT(R_fcall = lang2(VECTOR_ELT(parsedCmd, 0), R_NilValue));
 
  if (INTEGER(margin)[0] == 1) {

    PROTECT(data_chunk = allocVector(REALSXP, dimensions[1]));
    rchunk = REAL(data_chunk);

    for(i=0; i<my_end-my_start; i++) {

      /* Master process won't have data aligned in memory in this case */
      if (worldRank != MASTER_PROCESS) {

        for(j=0; j<dimensions[1]; j++) {
          rchunk[j] = rdata[i*dimensions[1]+j];
        }
        
      } else {

        for(j=0; j<dimensions[1]; j++) {
          rchunk[j] = rdata[j*dimensions[0]+i];
        }
      }
      
      SETCADR(R_fcall, data_chunk);
      SET_VECTOR_ELT(ans, i, eval(R_fcall, R_GlobalEnv));
      
      }    
  } else if (INTEGER(margin)[0] == 2) {

    PROTECT(data_chunk = allocVector(REALSXP, dimensions[0]));
    rchunk = REAL(data_chunk);
   
    for(i=0; i<my_end-my_start; i++) {
      
      for(j=0; j<dimensions[0]; j++) {
        rchunk[j] = rdata[j+i*dimensions[0]];
      }

      SETCADR(R_fcall, data_chunk);
      SET_VECTOR_ELT(ans, i, eval(R_fcall, R_GlobalEnv));

    }
  }
  
  UNPROTECT(3);
  return;
  
}
