/**************************************************************************
 *                                                                        *
 *  SPRINT: Simple Parallel R INTerface                                   *
 *  Copyright © 2012 The University of Edinburgh                          *
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
#include <string.h>
#include "../../../sprint.h"
#include "../../../functions.h"

extern int stringDist(int n,...);

/* **************************************************************** *
 *  Accepts information from R gets response and returns it.        *
 *                                                                  *
 * **************************************************************** */
void pstringDist(char **x,
             char **out_filename,
             int *sample_width,
             int *number_of_samples)
{

  int response,intCode;
  int worldSize, worldRank;
  
  enum commandCodes commandCode;
  
  // Check that MPI is initialized
  MPI_Initialized(&response);
  if (response) {
    DEBUG("\nMPI is init'ed in pstringDist\n");
  } else {
    DEBUG("\nMPI is NOT init'ed in pstringDist\n");

    // This is called from .C and so can't return a value.
    // The number_of_samples will be checked in the R code as an indication that the code has worked.
    // The code could be refactored to use .Call instead of .C to return a value.
    *number_of_samples = -1;
    return;
  }
  
  // Get size and rank from communicatorx
  MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
  MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);
  
  // Broadcast command to other processors
  commandCode = PSTRINGDIST;
  intCode = (int)commandCode;
  //MPI_BCAST (buffer(address of data to be sent), count (no. of elements in buffer),
  // datatype, root (rank of root process, comm)
  MPI_Bcast(&intCode, 1, MPI_INT, 0, MPI_COMM_WORLD);
  
  // Call the stringDist function
  // The 4 is for the number of arguments
  response = stringDist(4, *x, *sample_width, *number_of_samples, *out_filename);
}
