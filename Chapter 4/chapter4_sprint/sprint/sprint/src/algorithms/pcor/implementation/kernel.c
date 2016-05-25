/**************************************************************************
 *                                                                        *
 *  SPRINT: Simple Parallel R INTerface                                   *
 *  Copyright ? 2008,2009 The University of Edinburgh                     *
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

#include <math.h>
#include "kernel.h"
#include "../../../sprint.h"
#include <limits.h>

// Percent of extra memory to be allocated for the correlation coefficients vector.
// This extra memory is added to the approximate memory needed by a fair division of
// all work assignments accross the workers. Small work imbalancing is expexted to
// and this extra space is going to cover for it.
#define MEM_PERC 0.2

// Pointer to allocate the local vector for the correlation
// coefficients. Had some problems when declaring this pointer inside
// "correlationKernel"... better here, everything works fine
double *cor = NULL;

// Vectors for storing Pearson correlation parameters in order not to re-compute
// them every time
double *Sxx_vector = NULL;
double *mean_value_vectorX = NULL;
double *Syy_vector = NULL;
double *mean_value_vectorY = NULL;

/* *************************************** *
 *  Free all dynamically allocated memory  *
 * *************************************** */
void free_all(double *cor,
              int *blocklens,
              int *indices,
              double *mean_value_vectorX,
              double *Sxx_vector,
              double *mean_value_vectorY,
              double *Syy_vector);

/* ************************************************************* *
 *  Compute parameters used by Pearson method                    *
 *  It will transform the input array to more meaningful values  *
 *  and also create vectors for the mean values and variances    *
 * ************************************************************* */
void compute_parameters(double *dataMatrixX, double *dataMatrixY, int rows, int columns) {

    int i, j;
    
      // Compute the X mean values of each row
      for(i=0; i < rows; i++) {
        mean_value_vectorX[i] = 0.0;
        for(j=0; j < columns; j++) {
          mean_value_vectorX[i] += dataMatrixX[(i*columns)+j];
        }
        mean_value_vectorX[i] /= columns;
      }
      
      // Transform the data_X array
      for(i=0; i < rows; i++) {
        Sxx_vector[i] = 0.0;
        for(j=0; j < columns; j++) {
          dataMatrixX[(i*columns)+j] -= mean_value_vectorX[i];
          Sxx_vector[i] += (dataMatrixX[(i*columns)+j] * dataMatrixX[(i*columns)+j]);
        }
        Sxx_vector[i] /= (columns-1);
      }
      
    if (dataMatrixY != NULL) {
      // Compute the Y mean values of each row
      for(i=0; i < rows; i++) {
        mean_value_vectorY[i] = 0.0;
        for(j=0; j < columns; j++) {
          mean_value_vectorY[i] += dataMatrixY[(i*columns)+j];
        }
        mean_value_vectorY[i] /= columns;
      }
      
      // Transform the data_Y array
      for(i=0; i < rows; i++) {
        Syy_vector[i] = 0.0;
        for(j=0; j < columns; j++) {
          dataMatrixY[(i*columns)+j] -= mean_value_vectorY[i];
          Syy_vector[i] += (dataMatrixY[(i*columns)+j] * dataMatrixY[(i*columns)+j]);
        }
        Syy_vector[i] /= (columns-1);
      }
    }
}

/* ****************************** *
 *  Pearson correlation function  *
 * ****************************** */
double pearson(double *data, int row_x, int row_y, int size) {

    int i;
    double sxy = 0.0;
    double r;
    double *x, *y;

    x = &data[row_x*size];
    y = &data[row_y*size];

    // Correlation
    for (i=0; i<size; i++) {
        sxy += x[i] * y[i];   
    }

    sxy /= (size-1);

    // This *should* set r=NAN if sxx*syy = 0 *and*
    // sxy = 0 - this replicates R's behaviour
    r = sxy / sqrt(Sxx_vector[row_x]*Sxx_vector[row_y]);

    return r;
}

/* **************************************************** *
 *  Pearson correlation function for two input matrices *
 * **************************************************** */
double pearson_XY(double *dataMatrixX, double *dataMatrixY, int row_x, int row_y, int size) {

    int i;
    double sxy = 0.0;
    double r;
    double *x, *y;

    x = &dataMatrixX[row_x*size];
    y = &dataMatrixY[row_y*size];

    // correlation
    for (i=0; i<size; i++) {
        sxy += x[i] * y[i];
    }

    sxy /= (size-1);

    // This *should* set r=NAN if sxx*syy = 0 *and*
    // sxy = 0 - this replicates R's behaviour
    r = sxy / sqrt(Sxx_vector[row_x]*Syy_vector[row_y]);

    return r;
}

/* *************************** *
 *  Main computational kernel  *
 * *************************** */
int correlationKernel(int rank,
                      int size,
                      double* dataMatrixX,
                      double* dataMatrixY,
                      int columns,
                      int rows,
                      char *out_filename,
                      int distance_flag) {

    int local_check = 0, global_check = 0;
    int i = 0, j, taskNo;
    int err, count = 0;
    unsigned long long fair_chunk = 0, coeff_count = 0;
    unsigned int init_and_cleanup_loop_iter=0;
    unsigned long long cor_cur_size = 0;
    
    double start_time, end_time;

    // Variables needed by the Indexed Datatype
    MPI_Datatype coeff_index_dt;
    MPI_File fh;
    int *blocklens, *indices;

    MPI_Status stat;
    MPI_Comm comm = MPI_COMM_WORLD;

    // Master processor keeps track of tasks
    if (rank == 0) {

        // Make sure everything will work fine even if there are
        // less genes than available workers (there are size-1 workers
        // master does not count)
        if ( (size-1) > rows )
            init_and_cleanup_loop_iter = rows+1;
        else
            init_and_cleanup_loop_iter = size;

        // Start timer
        start_time = MPI_Wtime();

        // Send out initial tasks (remember you have size-1 workers, master does not count)
        for (i=1; i<init_and_cleanup_loop_iter; i++) {
            taskNo = i-1;
            err = MPI_Send(&taskNo, 1, MPI_INT, i, 0, comm);
        }        

        // Terminate any processes that were not working due to the fact
        // that the number of rows where less than the actual available workers
        for(i=init_and_cleanup_loop_iter; i < size; i++) {
            PROF(rank, "\nPROF_idle : Worker %d terminated due to insufficient work load", i);
            err = -1;
            err = MPI_Send(&err, 1, MPI_INT, i, 0, comm);
        }

        // Wait for workers to finish their work assignment and ask for more
        for (i=init_and_cleanup_loop_iter-1; i<rows; i++) {
            err = MPI_Recv(&taskNo, 1, MPI_INT, MPI_ANY_SOURCE, 0, comm, &stat);

            // Check taskNo to make sure everything is ok. Negative means there is problem
            // thus terminate gracefully all remaining working workers
            if ( taskNo < 0 ) {
                // Reduce by one because one worker is already terminated
                init_and_cleanup_loop_iter--;
                // Break and cleanup
                break;
            }

            // The sending processor is ready to work:
            // It's ID is in stat.MPI_SOURCE
            // Send it the current task (i)
            err = MPI_Send(&i, 1, MPI_INT, stat.MPI_SOURCE, 0, comm);
        }

        // Clean up processors
        for (i=1; i<init_and_cleanup_loop_iter; i++) {
            // All tasks complete - shutdown workers
            err = MPI_Recv(&taskNo, 1, MPI_INT, MPI_ANY_SOURCE, 0, comm, &stat);
            // If process failed then it will not be waiting to receive anything
            // We have to ignore the send because it will deadlock
            if ( taskNo < 0 )
                continue;
            err = -1;
            err = MPI_Send(&err, 1, MPI_INT, stat.MPI_SOURCE, 0, comm);
        }

        // Master is *always* OK
        local_check = 0;
        MPI_Allreduce(&local_check, &global_check, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

        // Check failed, abort
        if ( global_check != 0 ) {
            return -1;
        }
        
        // Stop timer
        end_time = MPI_Wtime();
        PROF(rank, "\nPROF_comp (workers=%d) : Time taken by correlation coefficients computations : %g\n", size-1, end_time - start_time);

        // Start timer
        start_time = MPI_Wtime();

        // Master process must call MPI_File_set_view as well, it's a collective call
        // Open the file handler
        MPI_File_open(comm, out_filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

        // Create the file view
        MPI_File_set_view(fh, 0, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);

        // Write data to disk
        MPI_File_write_all(fh, &cor[0], 0, MPI_DOUBLE, &stat);

        // Stop timer
        end_time = MPI_Wtime();
        PROF(rank, "\nPROF_write (workers=%d) : Time taken for global write-file : %g\n",  size-1, end_time - start_time);

    } else {

        // Compute how many workers will share the work load
        // Two scenarios exist:
        // (1) more OR equal number of workers and rows exist
        // (2) more rows than workers
        if ( (size-1) > rows ) {
            // For this scenario each worker will get exaclty one work asssignment.
            // There is not going to be any other work so it only compute "rows" number
            // of coefficients
            fair_chunk = rows;
            cor_cur_size = fair_chunk;
        } else {
            // For this scenario we are going to allocate space equal to a fair
            // distribution of work assignments *plus* an extra amount of space to
            // cover any load imbalancing. This amount is expressed as a percentage
            // of the fair work distribution (see on top, 20% for now)

            // Plus 1 to round it up or just add some extra space, both are fine
            fair_chunk = (rows / (size-1)) + 1;
            DEBUG("fair_chunk %d \n", fair_chunk);

            // We can use "j" as temporary variable.
            // Plus 1 to avoid getting 0 from the multiplication.
            j = (fair_chunk * MEM_PERC) + 1;

            cor_cur_size = (fair_chunk + j) * rows;
            DEBUG("cor_cur_size %lld \n", cor_cur_size);
        }

        // Allocate memory
        DEBUG("cor_cur_size %lld \n", cor_cur_size);
        long long double_size = sizeof(double);
        DEBUG("malloc size %lld \n", (double_size * cor_cur_size));
        cor = (double *)malloc(double_size * cor_cur_size);

        blocklens = (int *)malloc(sizeof(int) * rows);
        indices = (int *)malloc(sizeof(int) * rows);

        mean_value_vectorX = (double *)malloc(sizeof(double) * rows);
        Sxx_vector = (double *)malloc(sizeof(double) * rows);
        mean_value_vectorY = (double *)malloc(sizeof(double) * rows);
        Syy_vector = (double *)malloc(sizeof(double) * rows);

        // Check that all memory is successfully allocated
        if ( ( cor == NULL ) || ( blocklens == NULL ) || ( indices == NULL ) || 
             ( mean_value_vectorX == NULL ) || ( Sxx_vector == NULL ) ||
             ( mean_value_vectorY == NULL ) || ( Syy_vector == NULL ) ) {
            ERR("**ERROR** : Memory allocation failed on worker process %d. Aborting.\n", rank);

            // Free allocated memory
            free_all(cor, blocklens, indices, mean_value_vectorX, Sxx_vector, mean_value_vectorY, Syy_vector);

            // Let the master process know its aborting in order to terminate
            // the rest of the working workers
            // We have to receive a work assignment first and then terminate
            // otherwise the master will deadlock trying to give work to this worker
            err = MPI_Recv(&taskNo, 1, MPI_INT, 0, 0, comm, &stat);
            taskNo = -1;
            err = MPI_Send(&taskNo, 1, MPI_INT, 0, 0, comm);

            // This worker failed
            local_check = 1;
            MPI_Allreduce(&local_check, &global_check, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

            return -1;
        }

        // Compute necessary parameters for Pearson method
        // (this will transform the values of the input array to more meaningful data
        //  and save us from a lot of redundant computations)
        compute_parameters(dataMatrixX, dataMatrixY, rows, columns);

        // Main loop for workers. They get work from master, compute coefficients,
        // save them to their *local* vector and ask for more work
        for(;;) {
            // Get work
            err = 0;
            err = MPI_Recv(&taskNo, 1, MPI_INT, 0, 0, comm, &stat);

            // If received task is -1, function is terminated
            if ( taskNo == -1 )  break;

            // Check if there is enough memory to store the new coefficients, if not reallocate
            // the current memory and expand it by MEM_PERC of the approximated size
            if ( cor_cur_size < (coeff_count + rows) ) {
                PROF(0, "\n**WARNING** : Worker process %3d run out of memory and reallocates. Potential work imbalancing\n", rank);
                DEBUG("\n**WARNING** : Worker process %3d run out of memory and reallocates. Potential work imbalancing\n", rank);

                // Use j as temporary again. Add two (or any other value) to avoid 0.
                // (two is just a random value, you can put any value really...)
                j = (fair_chunk * MEM_PERC) + 2;
                cor_cur_size += (j * rows);

                // Reallocate and check
                cor = (double *)realloc(cor, sizeof(double) * cor_cur_size);
                if ( cor == NULL ) {
                    ERR("**ERROR** : Memory re-allocation failed on worker process %d. Aborting.\n", rank);

                    // Let the master process know its aborting in order to terminate
                    // the rest of the working workers
                    taskNo = -1;
                    err = MPI_Send(&taskNo, 1, MPI_INT, 0, 0, comm);

                    // This worker failed
                    local_check = 1;
                    MPI_Allreduce(&local_check, &global_check, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

                    // Free all allocated memory
                    free_all(cor, blocklens, indices, mean_value_vectorX, Sxx_vector, mean_value_vectorY, Syy_vector);

                    return -1;
                }
            }

            // Compute the correlation coefficients
            if(dataMatrixY != NULL) {
              for (j=0; j < rows; j++) {
                cor[coeff_count] = pearson_XY(dataMatrixX, dataMatrixY, j, taskNo, columns);
                coeff_count++;
              }

            } else {
              for (j=0; j < rows; j++) {
                // Set main diagonal to 1
                if ( j == taskNo ) {
                  cor[coeff_count] = 1.0;
                  coeff_count++;
                  continue;
                }
                cor[coeff_count] = pearson(dataMatrixX, taskNo, j, columns);
                coeff_count++;
              }
            }

            // The value of blocklens[] represents the number of coefficients on each
            // row of the corellation array
            blocklens[count] = rows;

            // The value of indices[] represents the offset of each row in the data file
            indices[count] = (taskNo * rows);
            count++;

            // Give the master the taskID
            err = MPI_Send(&taskNo, 1, MPI_INT, 0, 0, comm);
        }

        // There are two possibilities
        //   (a) everything went well and all workers finished ok
        //   (b) some processes finished ok but one or more of the remaining working workers failed
        // To make sure all is well an all-reduce will be performed to sync all workers and guarantee success
        // before moving on to write the output file
        // This worker is OK
        local_check = 0;
        MPI_Allreduce(&local_check, &global_check, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

        // Check failed
        if ( global_check != 0 ) {
            // Free all allocated memory
          free_all(cor, blocklens, indices, mean_value_vectorX, Sxx_vector, mean_value_vectorY, Syy_vector);
            return -1;
        }

        PROF(0, "\nPROF_stats (thread %3d) : Fair chunk of work : %d \t\t Allocated : %d \t\t Computed : %d\n",
                rank, fair_chunk, cor_cur_size, coeff_count);

        // If the distance_flag is set, then transform all correlation coefficients to distances
        if ( distance_flag == 1 ) {
            for(j=0; j < coeff_count; j++) {
                cor[j] = 1 - cor[j];
            }
        }

        // Create and commit the Indexed datatype *ONLY* if there are data available
        if ( coeff_count != 0 ) {
            MPI_Type_indexed(count, blocklens, indices, MPI_DOUBLE, &coeff_index_dt);
            MPI_Type_commit(&coeff_index_dt);
        }

        // Open the file handler
        MPI_File_open(comm, out_filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

        // Create the file view
        if ( coeff_count != 0 ) {
            MPI_File_set_view(fh, 0, MPI_DOUBLE, coeff_index_dt, "native", MPI_INFO_NULL);
        } else {
            MPI_File_set_view(fh, 0, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);
        }

        // Write data to disk
        // TODO coeff_count cannot be greater than max int (for use in the MPI_File_write_all call). 
        // A better fix should be possible, for now throw error.
        
        DEBUG("\ncoeff_count is %lld\n", coeff_count);
        DEBUG("\INT_MAX is %d\n", INT_MAX);
        if(coeff_count>INT_MAX)
        {
            ERR("**ERROR** : Could not run as the chunks of data are too large. Try running again with more MPI processes.\n");

            // Free allocated memory
            free_all(cor, blocklens, indices, mean_value_vectorX, Sxx_vector, mean_value_vectorY, Syy_vector);

            // Let the master process know its aborting in order to terminate
            // the rest of the working workers
            // We have to receive a work assignment first and then terminate
            // otherwise the master will deadlock trying to give work to this worker
            err = MPI_Recv(&taskNo, 1, MPI_INT, 0, 0, comm, &stat);
            taskNo = -1;
            err = MPI_Send(&taskNo, 1, MPI_INT, 0, 0, comm);

            // This worker failed
            local_check = 1;
            MPI_Allreduce(&local_check, &global_check, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

            return -1;
        }

        
        
        DEBUG("\nWriting %d to disk\n", coeff_count);

        MPI_File_write_all(fh, &cor[0], coeff_count, MPI_DOUBLE, &stat);

        if (coeff_count != 0 )
            MPI_Type_free(&coeff_index_dt);

        // Free all allocated memory
        free_all(cor, blocklens, indices, mean_value_vectorX, Sxx_vector, mean_value_vectorY, Syy_vector);
    }

         DEBUG("\nAbout to write to disk %d\n", rank);
    MPI_File_sync( fh ) ;   		// Causes all previous writes to be transferred to the storage device
         DEBUG("\nWritten to disk %d\n",rank);
  //  MPI_Barrier( MPI_COMM_WORLD ) ; 	// Blocks until all processes in the communicator have reached this routine.
         DEBUG("\nAfter barrier \n", rank);

    // Close file handler
    MPI_File_close(&fh);
  DEBUG("\nAfter file closed /n");
   // MPI_Barrier( MPI_COMM_WORLD ) ; 	// Blocks until all processes in the communicator have reached this routine.
      DEBUG("\nAbout to return from kernel /n");
      return 0;
}

/* *************************************** *
 *  Free all dynamically allocated memory  *
 * *************************************** */
void free_all(double *cor,
              int *blocklens,
              int *indices,
              double *mean_value_vectorX,
              double *Sxx_vector,
              double *mean_value_vectorY,
              double *Syy_vector) {
    if ( cor != NULL ) free(cor);
    if ( blocklens != NULL ) free(blocklens);
    if ( indices != NULL ) free(indices);
    if ( mean_value_vectorX != NULL ) free(mean_value_vectorX);
    if ( Sxx_vector != NULL ) free(Sxx_vector);
    if ( mean_value_vectorY != NULL ) free(mean_value_vectorY);
    if ( Syy_vector != NULL ) free(Syy_vector);
}

