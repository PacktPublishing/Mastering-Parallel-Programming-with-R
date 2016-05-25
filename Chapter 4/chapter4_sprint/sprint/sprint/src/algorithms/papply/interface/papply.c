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

#include <Rdefines.h>
#include "../../../sprint.h"
#include "../../../functions.h"

extern int apply(int n,...);

/* ******************************************************** *
 *                    *
 *                    *
 * ******************************************************** */

SEXP papply(SEXP data,
            SEXP margin,
            SEXP function,
            SEXP nrows, /* ff matrix dimensions */
            SEXP ncols, /* if no ff object was passed this equals to 0 */
            SEXP out_filename)
{
  SEXP result = NULL;
  int i, response, worldSize, worldRank, intCode;
  enum commandCodes commandCode;
  
  // Check that MPI is initialized 
  MPI_Initialized(&response);
  if (response) {
    DEBUG("MPI is init'ed in papply\n");
  } else {
    
    DEBUG("MPI is NOT init'ed in papply\n");
    PROTECT(result = NEW_INTEGER(1));
    INTEGER(result)[0] = -1;
    UNPROTECT(1);
    
    return result;
  }

   // Get size and rank from communicator
  MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
  MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);
  
  /* Check the format of the input data */
  if(!isString(data)) {
    
    if(isMatrix(data)) {
      if (!IS_NUMERIC(data) && !IS_INTEGER(data)) {
        ERR("\npapply accepts only numeric data\n");
        result = NEW_INTEGER(1);
        INTEGER(result)[0] = -1;
        return result;
      }
    } else if (isNewList(data)) {
      
      for(i = 0; i < length(data); i++) {
        if (!IS_NUMERIC(VECTOR_ELT(data, i)) && !IS_INTEGER(VECTOR_ELT(data, i))) {
          ERR("\nlist element [%d] error: papply accepts only numeric data\n", i+1);
          result = NEW_INTEGER(1);
          INTEGER(result)[0] = -1;
          return result;
        }
      }
    } else {
      
      ERR("\ndata format not supported papply accepts only matrices and lists\n");
      result = NEW_INTEGER(1);
      INTEGER(result)[0] = -1;
      return result;      
    }
  }

  /* Check if function is serialised to char */
  if (!IS_CHARACTER(function) && !GET_LENGTH(function) > 0) {
    ERR("\nfun argument is invalid\n");
    result = NEW_INTEGER(1);
    INTEGER(result)[0] = -1;
    return result;
  }
  
  // broadcast command to other processors
  commandCode = PAPPLY;
  intCode = (int)commandCode;
  MPI_Bcast(&intCode, 1, MPI_INT, 0, MPI_COMM_WORLD);

  PROTECT(result = allocVector(VECSXP, 1));

  DEBUG("About to call apply\n");
  response  = apply(6,                
                    data,
                    margin,
                    function,
                    result,
                    nrows,
                    ncols,
                    out_filename);
  
  UNPROTECT(1);
  return result;
}
