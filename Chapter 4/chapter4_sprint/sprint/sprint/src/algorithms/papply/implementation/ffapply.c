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
#include "../../common/mmap.h"
#include "ffapply.h"
#include "comms.h"

SEXP ffApply(SEXP result, SEXP data, SEXP margin, SEXP function,
             SEXP nrows, SEXP ncols, int worldRank,
             SEXP out_filename, int worldSize) {

	  DEBUG("In ffapply\n");
  SEXP ans;
  
  int my_start, my_end, N, function_nlines;
  int local_check = 0, global_check = 0;

  char *filename, *file_out;
  int  *filesize;
  double *mapped_data_matrix;
  
  filesize = (int *) R_alloc(1, sizeof(int));

  if(worldRank == MASTER_PROCESS) {
    /* data argument is actually a path to a binary file where data is stored */
    filename = (char *)CHAR((STRING_ELT(data,0)));
    file_out = (char *)CHAR(STRING_ELT(out_filename,0));

    /* function SEXP object is a vector of strings, each element contains
       a single line of the function definition */
    function_nlines = length(function);
    
  } else {
    filename = (char *) R_alloc(FILENAME_LENGTH, sizeof(char));
    file_out = (char *) R_alloc(FILENAME_LENGTH, sizeof(char));

    DEBUG("\n**PROTECTING 3 : Worker process %3d \n", worldRank);
    PROTECT(nrows = allocVector(INTSXP, 1));
    PROTECT(ncols = allocVector(INTSXP, 1));
    PROTECT(margin = allocVector(INTSXP, 1));
  }

  MPI_Bcast(filename, FILENAME_LENGTH, MPI_CHAR, 0, MPI_COMM_WORLD);
  MPI_Bcast(file_out, FILENAME_LENGTH, MPI_CHAR, 0, MPI_COMM_WORLD); 
  MPI_Bcast(INTEGER(nrows), 1, MPI_INT, 0, MPI_COMM_WORLD); 
  MPI_Bcast(INTEGER(ncols), 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(INTEGER(margin), 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&function_nlines, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if(worldRank != MASTER_PROCESS) {

	    DEBUG("\n**PROTECTING 1 : Worker process %3d \n", worldRank);
    PROTECT(function = allocVector(STRSXP, function_nlines));
  }

  if((mapped_data_matrix = map_file(filename, filesize)) == NULL) {
       local_check = -1;
  }
  
  /* Check if all processes have successfully mapped the file to memory */
  MPI_Allreduce(&local_check, &global_check, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if ( global_check != 0 ) return ScalarInteger(-1);

  /* Matrix dimensions in R are interpreted differently than in C.
     We will refer to R rows and columns ordering, so rows are not alligned
     in memory */

  if (INTEGER(margin)[0] == 1) {
    N = INTEGER(nrows)[0];
  } else if (INTEGER(margin)[0] == 2) {
    N = INTEGER(ncols)[0];
  } else {
    DEBUG("Do not know how to distribute margin number %d\n",
          INTEGER(margin)[0]);
    return R_NilValue;
  }

  loopDistribute(worldRank, worldSize, N, &my_start, &my_end);

  /* Bcast function name or definition, cover case when definition is split into
     several lines and stored as a SEXP string vector */
  bcastRFunction(function, function_nlines, worldRank);

  DEBUG("\n**PROTECTING 1 : Worker process %3d \n", worldRank);
  /* Response container, Vector of SEXPs, margin determines vector length */
  PROTECT(ans = allocVector(VECSXP, 1));

  do_ffApply(ans, mapped_data_matrix, margin, function, my_start, my_end,
             INTEGER(nrows)[0], INTEGER(ncols)[0], worldRank, file_out);

  if(worldRank != MASTER_PROCESS) {

	    DEBUG("\n**UNPROTECTING 5 : Worker process %3d \n", worldRank);
    UNPROTECT(5);
  } else {
	    DEBUG("\n**UNPROTECTING 1 : Worker process %3d \n", worldRank);
    UNPROTECT(1);
  }

  return result;
}

void do_ffApply(SEXP ans,
                double *data,
                SEXP margin,
                SEXP function,
                int my_start,
                int my_end,
                int nrows,
                int ncols,
                int worldRank,
                char *out_filename) {

  MPI_Status stat;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_File fh;

  SEXP data_chunk, R_fcall, parsedCmd = R_NilValue;  
  double *rchunk;
  int i,k, offset=0, count=0;
  int no_of_protects = 0;
  
  /* Open the file handler */
  MPI_File_open(comm, out_filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

  /* The MPI_FILE_SET_VIEW routine changes the process's view of the data in the file */
  MPI_File_set_view(fh, 0, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);

  /* Parse the command, returns a function object */
  DEBUG("\n**PROTECTING 2 : Worker process %3d \n", worldRank);
  PROTECT(parsedCmd = parseExpression(function));
  no_of_protects++;

   /* Create R LANGSXP Vector, R function holder
     length of the vector is 1 + number of arguments */
  PROTECT(R_fcall = lang2(VECTOR_ELT(parsedCmd, 0), R_NilValue));
  no_of_protects++;


  if (INTEGER(margin)[0] == 1) {

	  DEBUG("\n**PROTECTING 1 : Worker process %3d \n", worldRank);
    PROTECT(data_chunk = allocVector(REALSXP, ncols));
    no_of_protects++;
    rchunk = REAL(data_chunk);

    for(i=my_start, k=0; i<my_end; i++, k++) {
      
      for(int j=0; j<ncols; j++) {
        rchunk[j] = data[j*nrows+i];
      }

      SETCADR(R_fcall, data_chunk);
      SET_VECTOR_ELT(ans, 0, eval(R_fcall, R_GlobalEnv));

      count = length(VECTOR_ELT(ans, 0));
      offset = i*count;
      DEBUG("\n**MPI_File_write_at** : Worker process %3d \n", worldRank);
      MPI_File_write_at(fh, offset, REAL(VECTOR_ELT(ans, 0)), count, MPI_DOUBLE, &stat); 
    }
  }

  if (INTEGER(margin)[0] == 2) {

	  DEBUG("\n**PROTECTING 1 : Worker process %3d \n", worldRank);
    PROTECT(data_chunk = allocVector(REALSXP, nrows));
    no_of_protects++;
    rchunk = REAL(data_chunk);

    for(i=my_start, k=0; i<my_end; i++, k++) {
      
      for(int j=0; j<nrows; j++) {
        rchunk[j] = data[j+nrows*i];
      }

      SETCADR(R_fcall, data_chunk);
      SET_VECTOR_ELT(ans, 0, eval(R_fcall, R_GlobalEnv));

      count = length(VECTOR_ELT(ans, 0));
      offset = i*count;

      DEBUG("\n**MPI_File_write_at** : Worker process %3d \n", worldRank);
      MPI_File_write_at(fh, offset, REAL(VECTOR_ELT(ans, 0)), count, MPI_DOUBLE, &stat);

    }

  }

  DEBUG("\n**count is %3d \n", count);
  DEBUG("\n**length(data_chunk) is %3d \n", length(data_chunk)) ;
  //TODO ET
  /* If count and length(data_chunk) are the same, then a matrix is written to the ff file,
   * if count is length ==1, then a list is written to the ff file (the current ff read code works for this).
   * Need to return a list/matix flag to the R code, so that it knows how to open the ff file.*/

  MPI_File_sync( fh ) ; 			// Causes all previous writes to be transferred to the storage device
  MPI_Barrier( MPI_COMM_WORLD ) ; 	// Blocks until all processes in the communicator have reached this routine.

  DEBUG("\n**No of things that are protected %d : Worker process %3d \n", no_of_protects, worldRank);
  DEBUG("\n**UNPROTECTING 3 : Worker process %3d \n", worldRank);
  UNPROTECT(3);
  
  /* Close file handler */
  MPI_File_close(&fh);

  MPI_Barrier( MPI_COMM_WORLD ) ; 	// Blocks until all processes in the communicator have reached this routine.
  return;
}

