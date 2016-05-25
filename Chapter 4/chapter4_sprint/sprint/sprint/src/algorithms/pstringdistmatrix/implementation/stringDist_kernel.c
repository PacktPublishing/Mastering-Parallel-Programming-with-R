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

#include "../../../sprint.h"
#include "../../common/utils.h"

int computeStringDist(int worldRank, int worldSize, char* DNAStringSet, int *stringDist, char *out_filename,
                   int sample_width, int number_of_samples, int my_start, int my_end, int chunk_size) {

  MPI_Status stat;
  MPI_File fh;
  int offset=0, count=0;
  int diss = my_end-my_start;

  int i,j,k,c,diff;

  /* Open the file handler */
  MPI_File_open(MPI_COMM_WORLD, out_filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

  /* The MPI_FILE_SET_VIEW routine changes the process's view of the data in the file */
  MPI_File_set_view(fh, 0, MPI_INT, MPI_INT, "native", MPI_INFO_NULL);
  
  for (i=my_start,c=0; i<my_end; i++,c++) {
    for(j=0; j<number_of_samples; j++) {

      diff = 0;
      for(k=0; k<sample_width; k++) {
        
        if(DNAStringSet[i*sample_width+k] != DNAStringSet[j*sample_width+k])
          diff++;

      }
      stringDist[(c*number_of_samples)+j] = diff;
    }

    if(c==chunk_size) {

      MPI_File_write_at(fh, (my_start*number_of_samples)+(chunk_size*number_of_samples*count), stringDist, number_of_samples*chunk_size, MPI_INT, &stat);
      count++;
      c=0;
    }
    
  }
  /*
  Input Parameters MPI_File_write_all

  fh
  File handle (handle).
  buf
  Initial address of buffer (choice).
  count
  Number of elements in buffer (integer).
  datatype
  Data type of each buffer element (handle). */

 // MPI_File_write_all(fh, &cor[0], 0, MPI_DOUBLE, &stat); from pcor
 // MPI_File_write_all(fh, &cor[0], coeff_count, MPI_DOUBLE, &stat);

//  MPI_File_write_all(fh, stringDist, number_of_samples*c, MPI_INT, &stat);
  MPI_File_write_at(fh, (my_start*number_of_samples)+(chunk_size*number_of_samples*count), stringDist, number_of_samples*c, MPI_INT, &stat);

 /*
  * Input Parameters MPI_File_write_at

  fh
  File handle (handle).
  offset
  File offset (integer).
  buf
  Initial address of buffer (choice).
  count
  Number of elements in buffer (integer).
  datatype
  Data type of each buffer element (handle).*/
  MPI_File_sync( fh ) ; 			// Causes all previous writes to be transferred to the storage device
  MPI_Barrier( MPI_COMM_WORLD ) ; 	// Blocks until all processes in the communicator have reached this routine.

  /* Close file handler */
  MPI_File_close(&fh);

  MPI_Barrier( MPI_COMM_WORLD ) ; 	// Blocks until all processes in the communicator have reached this routine.
  return 0;
}


int allocateMaxChunk(int worldRank, int my_start, int my_end, int *stringDist,
                     int number_of_samples) {

  int chunk_size = (my_end-my_start);

/* First attempt at 1GB chunk sizes. Does not produce correct results though.
  int size_on_disk = sizeof(int)* number_of_samples * chunk_size;

  int ONE_GB = 100;//1073741824;
  	  	  	    //484 000 000
  Rprintf("chunk_size at beginning: %d\n", chunk_size);

  Rprintf("size_on_disk before check: %d\n", size_on_disk);

  while(size_on_disk > ONE_GB){
	  Rprintf("size_on_disk: %d\n", size_on_disk);
      chunk_size = chunk_size/2;
	  Rprintf("halving chunk size to: %d\n", chunk_size);
      size_on_disk = sizeof(int)* number_of_samples * chunk_size;
	  Rprintf("new size_on_disk: %d\n", size_on_disk);
  }
*/
  int malloc_failed=1;
    
  // This malloc fail handler is probably not used, but leaving it in just in case.
  if(NULL == (stringDist = (int *)malloc(sizeof(int)
                                              * number_of_samples * chunk_size)))
  {
    while(malloc_failed) { // This does not work malloc will allocate more virtual memory than is available.
      malloc_failed=0;
      chunk_size = chunk_size/2;
      
      stringDist = (int *)malloc(sizeof(int) * number_of_samples * chunk_size);
      
      if(stringDist == NULL) {
        malloc_failed=1;
      }
      free(stringDist);
    }
  }
  free(stringDist);

  return chunk_size;

}
