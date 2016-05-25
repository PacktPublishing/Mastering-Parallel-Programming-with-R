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

#include "stringDist_kernel.h"
#include "../../../sprint.h"
#include "../../common/utils.h"

#define FILENAME_SIZE 256

int stringDist(int n, ...) {

  va_list ap; /*will point to each unnamed argument in turn*/
  int worldSize, worldRank, response;
  int sample_width, number_of_samples;
  
  int local_check=0, global_check=0;
  int my_start=0, my_end=0, chunk_size=0;
  
  char *DNAStringSet = NULL;
  char *out_filename = NULL;
  
  int *stringDist;

  // Get size and rank from communicator
  MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
  MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);

  if (worldRank == 0) {
    if (n == 4) {
      
      // Get input variables
      va_start(ap,n);
      DNAStringSet = va_arg(ap,char*);
      sample_width = va_arg(ap,int);
      number_of_samples = va_arg(ap,int);
      out_filename = va_arg(ap,char*);
      va_end(ap);
    } else {
      DEBUG("rank 0 passed incorrect arguments into correlation!");
    }
    
    // Sent the dimensions to the slave processes
    MPI_Bcast(&sample_width, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&number_of_samples, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    // Sent the number of arguments
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    // Master is always OK
    local_check = 0;
  } else {
    
    // Get the dimensions from the master
    MPI_Bcast(&sample_width, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&number_of_samples, 1, MPI_INT, 0, MPI_COMM_WORLD);
    DEBUG("Broadcasting width and length on %i. Got %i, %i\n", worldRank, sample_width, number_of_samples);
    
    // Get the number of arguments
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    DEBUG("Broadcasting number of arguments on %i. Got %i\n", worldRank, n);
    
    // Allocate memory for the input data array and out_filename
    DNAStringSet = (char *)malloc(sizeof(char) * sample_width * number_of_samples);
    out_filename = (char *)malloc(sizeof(char) * FILENAME_SIZE);
    
    // Check memory and make sure all slave processes have allocated successfully.
    // Perform an all-reduce operation to make sure everything is ok and then
    // move on to broadcast the data
    if ( DNAStringSet == NULL || out_filename == NULL) {
      local_check = 1;
      ERR("**ERROR** : Input data array memory allocation failed on slave process %d. Aborting.\n", worldRank);
    }
    else {
      local_check = 0;
    }
  }

  MPI_Allreduce(&local_check, &global_check, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if ( global_check != 0 && worldRank != 0 ) {
    if ( out_filename != NULL ) free(out_filename);
    if ( DNAStringSet != NULL ) free(DNAStringSet);
    return -1;
  }
  else if ( global_check != 0 && worldRank == 0 )
    return -1;

  // Broadcast the output filename
  MPI_Bcast(out_filename, FILENAME_SIZE, MPI_CHAR, 0, MPI_COMM_WORLD);
  DEBUG("Broadcasting output filename on %i. Got %s\n", worldRank, out_filename);
  
  MPI_Bcast(DNAStringSet, number_of_samples*sample_width, MPI_CHAR, 0, MPI_COMM_WORLD);
  DEBUG("Broadcasting DNAStringSet on %i.\n", worldRank);
  
  loopDistribute(worldRank, worldSize, number_of_samples, &my_start, &my_end);
  DEBUG("loopDistribute results on %i, %i %i\n", worldRank, my_start, my_end);

  chunk_size = allocateMaxChunk(worldRank, my_start, my_end, stringDist, number_of_samples);
  stringDist = (int *)malloc(sizeof(int) * number_of_samples * chunk_size);

  response = computeStringDist(worldRank, worldSize, DNAStringSet, stringDist, out_filename,
                            sample_width, number_of_samples, my_start, my_end, chunk_size);

  DEBUG("Done running stringDist kernel on %i\n", worldRank);

  // Free memory allocated for output filename and data array on slave processes
  if ( worldRank != 0 ) {
    
    free(out_filename);
    free(DNAStringSet);
    free(stringDist);
  }
  
  return (0);
}
