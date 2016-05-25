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
#include <R_ext/Parse.h>

#include "../../../sprint.h"
#include "../../common/utils.h"
#include "lapply.h"
#include "comms.h"


SEXP listApply(SEXP result, SEXP data, SEXP function,
               int worldRank, int worldSize) {

  SEXP ans, list_element, data_size;

  MPI_Status status;
  
  int my_start, my_end, function_nlines,
    list_length, nelements,
    worker_start, worker_end;
  int dimensions[2];

  if (worldRank == MASTER_PROCESS) {

    list_length = length(data);
    
    /* Function SEXP object is a vector of strings, each element contains
       a single line of a function definition */
    function_nlines = length(function);
  }
  /* Bcast number of elements in the list and length of function def */
  MPI_Bcast(&list_length, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&function_nlines, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (worldRank != MASTER_PROCESS) {
      /* Allocate a vector to hold multiline function definitions */
      PROTECT(function = allocVector(STRSXP, function_nlines));
  }
   /* Bcast function name or definition, cover cases when definition is split into
     several lines and stored as a SEXP string vector */
  bcastRFunction(function, function_nlines, worldRank); 

   /* Distribute data between processes */
   for (int worker_id=1; worker_id<worldSize; worker_id++) {

     if (worldRank == MASTER_PROCESS) {

       /* Calculate position and number of elements in the list to be sent
          to the worker process */
       loopDistribute(worker_id, worldSize, list_length, &worker_start, &worker_end);
       nelements = worker_end - worker_start;
       
       for (int i=worker_start; i<worker_end; i++) {
         list_element = VECTOR_ELT(data,i);
         data_size = GET_DIM(list_element);
         
         dimensions[0] = INTEGER_POINTER(data_size)[0];
         dimensions[1] = INTEGER_POINTER(data_size)[1];
         
         /* Send dimensions of each data element to the worker before sending
            the data */
         MPI_Send(dimensions, 2, MPI_INT, worker_id, 0, MPI_COMM_WORLD);
         MPI_Send(REAL(list_element), dimensions[0]*dimensions[1], MPI_DOUBLE, worker_id, 0, MPI_COMM_WORLD);
         
       }
     } else if (worldRank == worker_id) {
       /* loopDistribute calculates distribution of list_length elements based on process id
          and total number of processes */
       loopDistribute(worldRank, worldSize, list_length, &my_start, &my_end);
       nelements = my_end - my_start;

       /* Create a vector to hold designated part of the input list */
       PROTECT(data = allocVector(VECSXP, nelements));

       for (int i=0; i<nelements; i++) {

         /* Get dimensions of the data element and allocate memory accordingly */
         MPI_Recv(dimensions, 2, MPI_INT, MASTER_PROCESS, 0, MPI_COMM_WORLD, &status);
         PROTECT(list_element = allocMatrix(REALSXP, dimensions[0], dimensions[1]));

         /* Receive data element and write to newly created SEXP matrix */
         MPI_Recv(REAL(list_element), dimensions[0]*dimensions[1], MPI_DOUBLE, MASTER_PROCESS, 0, MPI_COMM_WORLD, &status);

         /* Set a new data element, dereference list_element object*/
         SET_VECTOR_ELT(data, i, list_element);
         UNPROTECT(1);

       }
     }
   }

  /* Calc loop distribution for the Master process */
   if (worldRank == MASTER_PROCESS) {
     loopDistribute(worldRank, worldSize, list_length, &my_start, &my_end);
   }
  
    /* Response container, Vector of SEXPs, list determines vector length */
   PROTECT(ans = allocVector(VECSXP, list_length));

   do_listApply(ans, data, function, my_end - my_start);

   gatherData(result, ans, list_length, my_start, my_end, worldRank);
   
   if (worldRank != MASTER_PROCESS) {
     UNPROTECT(3);     
   } else {
     UNPROTECT(1);
   }

   return result;
}

void do_listApply(SEXP ans,
                  SEXP data,
                  SEXP function,
                  int chunk_length)
{

  SEXP R_fcall, parsedCmd = R_NilValue;
  
  /* Parse the expression, creates function object that can be passed to
     R eval*/
  PROTECT(parsedCmd = parseExpression(function));

  /* Create R LANGSXP Vector, R function holder
     length of the vector is 1 + number of arguments */
  PROTECT(R_fcall = lang2(VECTOR_ELT(parsedCmd, 0), R_NilValue));

  for (int i=0;  i<chunk_length; i++) {
    SETCADR(R_fcall, VECTOR_ELT(data, i));
    SET_VECTOR_ELT(ans, i, eval(R_fcall, R_GlobalEnv));
  }
    
  UNPROTECT(2);
  return;
  
}


