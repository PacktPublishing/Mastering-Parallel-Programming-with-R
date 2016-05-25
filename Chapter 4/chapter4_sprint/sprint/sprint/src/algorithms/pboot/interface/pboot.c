/**************************************************************************
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

#include <stdarg.h>
#include <Rdefines.h>
#include "../../../sprint.h"
#include "../../../functions.h"

extern int boot(int n, ...);
  
SEXP pboot(SEXP list, SEXP fn)
{
  SEXP result;
  int i, response, worldSize, worldRank, intCode;
  enum commandCodes commandCode;
  
  // Check that MPI is initialized 
  MPI_Initialized(&response);
  if (response) {
    DEBUG("MPI is init'ed in ptest\n");
  } else {
    
    DEBUG("MPI is NOT init'ed in ptest\n");
    PROTECT(result = NEW_INTEGER(1));
    INTEGER(result)[0] = -1;
    UNPROTECT(1);
    
    return result;
  }

   // Get size and rank from communicator
  MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
  MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);

  // broadcast command to other processors
  commandCode = PBOOT;
  intCode = (int)commandCode;
  MPI_Bcast(&intCode, 1, MPI_INT, 0, MPI_COMM_WORLD);

  response = boot(3, list, fn, &result);

  return result;
  
}

