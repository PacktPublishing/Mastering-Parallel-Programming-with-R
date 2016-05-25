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

#include <stdarg.h>
#include "kernel.h"
#include "../../../sprint.h"

#define WIDTH 0
#define HEIGHT 1
#define FILENAME_SZ 256
#define MATRIX_Y_PRESENT 6


int correlation(int n,...) {
    int result = 0;
    va_list ap; /*will point to each unnamed argument in turn*/
    int worldSize, worldRank;
    double *dataMatrixX = NULL;
    double *dataMatrixY = NULL;
    int dimensions[2];
    char *out_filename = NULL;
    int local_check = 0, global_check = 0;
    int distance_flag=0;
    double start_time, end_time;

    // Get size and rank from communicator
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);
    
    DEBUG("Hello from %i/%i\n", worldRank, worldSize);
    
    if (worldRank == 0) {
        if (n == 5) {
          
          // Get input variables
          va_start(ap,n);
          dataMatrixX = va_arg(ap,double*);
          dimensions[0] = va_arg(ap,int);
          dimensions[1] = va_arg(ap,int);
          out_filename = va_arg(ap,char*);
          distance_flag = va_arg(ap,int);
          va_end(ap);
        } else if (n == 6) {
          
          // Get input variables
          va_start(ap,n);
          dataMatrixX = va_arg(ap,double*);
          dataMatrixY = va_arg(ap,double*);
          dimensions[0] = va_arg(ap,int);
          dimensions[1] = va_arg(ap,int);
          out_filename = va_arg(ap,char*);
          distance_flag = va_arg(ap,int);
          va_end(ap);
        } else {

          DEBUG("rank 0 passed incorrect arguments into correlation!");
        }
        
        // Sent the dimensions to the slave processes
        MPI_Bcast(dimensions, 2, MPI_INT, 0, MPI_COMM_WORLD);

        // Sent the number of arguments
        MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

        // Master is always OK
        local_check = 0;
    } else {

        // Get the dimensions from the master
        MPI_Bcast(dimensions, 2, MPI_INT, 0, MPI_COMM_WORLD);
        DEBUG("Broadcasting dims on %i. Got %i, %i\n", worldRank, dimensions[WIDTH], dimensions[HEIGHT]);

        // Get the number of arguments
        MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
        DEBUG("Broadcasting number of arguments on %i. Got %i\n", worldRank, n);

        // Allocate space for output filename
        out_filename = (char *)malloc(sizeof(char) * FILENAME_SZ);

        // Allocate memory for the input data array
        dataMatrixX = (double *)malloc(sizeof(double) * dimensions[WIDTH] * dimensions[HEIGHT]);
        if (n == MATRIX_Y_PRESENT) {
          dataMatrixY = (double *)malloc(sizeof(double) * dimensions[WIDTH] * dimensions[HEIGHT]);
        }
        
        // Check memory and make sure all slave processes have allocated successfully.
        // Perform an all-reduce operation to make sure everything is ok and then
        // move on to broadcast the data
        if ( (out_filename == NULL) || (dataMatrixX == NULL) || (n == MATRIX_Y_PRESENT && dataMatrixY == NULL)) {
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
        if ( dataMatrixX != NULL ) free(dataMatrixX);
        if ( dataMatrixY != NULL ) free(dataMatrixY);
        return -1;
    }
    else if ( global_check != 0 && worldRank == 0 )
        return -1;

    // Start timer
    start_time = MPI_Wtime();

    // Broadcast the output filename
    MPI_Bcast(out_filename, FILENAME_SZ, MPI_CHAR, 0, MPI_COMM_WORLD);
    DEBUG("Broadcasting output filename on %i. Got %s\n", worldRank, out_filename);

    DEBUG("Broadcasting distance flag on %i.\n", worldRank);
    MPI_Bcast(&distance_flag, 1, MPI_INT, 0, MPI_COMM_WORLD);

    DEBUG("Broadcasting data X on %i.\n", worldRank);
    MPI_Bcast(dataMatrixX, dimensions[WIDTH] * dimensions[HEIGHT], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (n == MATRIX_Y_PRESENT) {
      DEBUG("Broadcasting data Y on %i.\n", worldRank);
      MPI_Bcast(dataMatrixY, dimensions[WIDTH] * dimensions[HEIGHT], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
    
    // Stop timer
    end_time = MPI_Wtime();

    PROF(worldRank, "\nPROF_bcast (slaves=%d) : Time taken for bcast (filename, dims, data) : %g\n",
         worldSize-1, end_time - start_time);

    DEBUG("Running correlation kernel on %i\n", worldRank);

    result = correlationKernel(worldRank, worldSize, dataMatrixX, dataMatrixY, dimensions[WIDTH], dimensions[HEIGHT], out_filename, distance_flag);

    DEBUG("Done running correlation kernel on %i\n", worldRank);

    // Free memory allocated for output filename and data array on slave processes
    if ( worldRank != 0 ) {
        free(dataMatrixX);
        free(dataMatrixY);
        free(out_filename);
    }
    
    return result;
}

