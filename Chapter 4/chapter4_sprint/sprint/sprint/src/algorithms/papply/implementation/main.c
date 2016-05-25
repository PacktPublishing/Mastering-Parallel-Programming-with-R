/**************************************************************************
 *                                                                        *
 *  SPRINT: Simple Parallel R INTerface                                   *
 *  Copyright © 2008,2009 The University of Edinburgh                     *
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

#include <stdio.h>
#include <stdarg.h>
#include <mpi.h>

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <Rdefines.h>

#include "../../../sprint.h"
#include "main.h"
#include "apply.h"
#include "lapply.h"
#include "ffapply.h"
#include "comms.h"

#define NAME_LENGTH 256

/* Declares an enumeration data type called R_dataFormat */
enum R_dataFormat { 
  R_matrix,
  R_list,
  R_ff
};

/**
 ** 
 **/

int apply(int n,...) {

  int worldSize;
  int worldRank, dataFormatInt;
  va_list ap;

  SEXP data, function, margin, result,
    nrows, ncols, out_filename;

  enum R_dataFormat dataFormat;

  MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
  MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);

  if (worldRank == 0) {
    if (n == 6) {
      
      // Get input variables
      va_start(ap,n);
      data = va_arg(ap, SEXP);
      margin = va_arg(ap, SEXP);
      function = va_arg(ap, SEXP);
      result = va_arg(ap, SEXP);
      nrows = va_arg(ap,SEXP);
      ncols = va_arg(ap,SEXP);
      out_filename = va_arg(ap,SEXP);
      
      va_end(ap);

      if(isMatrix(data)) {
        dataFormat = R_matrix;
      } else if (isNewList(data)) {
        dataFormat = R_list;
      } else if (isString(data)) {
        dataFormat = R_ff;
      }
      else {          
        DEBUG("rank 0 passed incorrect arguments into correlation!\n");
      }
    } else {
        return -1;
    }
  } else {
      out_filename = NULL;
      ncols = NULL;
      nrows = NULL;
      result = NULL;
      margin = NULL;
      function = NULL;
      data = NULL;
  }

  dataFormatInt = (int)dataFormat;

  MPI_Bcast(&dataFormatInt, 1, MPI_INT, 0, MPI_COMM_WORLD);

  dataFormat = dataFormatInt;

  /* papply accept data in three different formats... */

  switch(dataFormat) {
    
  case R_matrix:

    matrixApply(result, data, margin, function, worldRank, worldSize);
    break;

  case R_list:

    listApply(result, data, function, worldRank, worldSize);
    break;

  case R_ff:

	  DEBUG("About to call ffApply\n");
    ffApply(result, data, margin, function, nrows,
            ncols, worldRank, out_filename, worldSize);
    break;

  default:
    ERR("Invalid data type passed to apply.c\n");
    break;
  }
   return 0;                 
}
