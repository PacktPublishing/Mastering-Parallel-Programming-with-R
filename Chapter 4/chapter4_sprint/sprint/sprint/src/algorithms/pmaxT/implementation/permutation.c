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
#include "../../../sprint.h"
#include "mt.h"

// Function to free all dynamically allocated memory
void perm_free(char **options, int *pnrow, int *pncol, int *pB, int *extra, int *L, double *pna,
                   double *data, double *T, int *index, int *generator_flag, int *res);

// Function to check all dynamically allocated memory
int perm_alloc_and_check(int *sizes, char ***options, int **pnrow, int **pncol, int **pB, int **extra,
                         int **L, double **pna, double **data, double **T, int **index, int **generator_flag,
                         int **total1, int **total2, int **count1, int **count2);

// Allocate and initialize vectors for storing the information necessary for computing
// the raw and ajusted values
void alloc_init_info_vectors(int vector_size, int **total1, int **total2, int **count1, int **count2);

// Function to return the number of iterations the process should perform
int split_iterations(int rank, int size, int pB);

/* ************************************************************************************************************************** *
 *  double  *d                  -->  all input data                                                                           *
 *  int     *pnrow              -->  scalar integer value                                                                     *
 *  int     *pncol              -->  scalar integer value                                                                     *
 *  int     *L                  -->  vector of integer values of size: (1) pncol **OR** (2) pncol/2 --> if(test == "pairt")   *
 *  double  *pna                -->  scalar double value                                                                      *
 *  double  *T                  -->  allocated memory for *all* T values         |                                            *
 *  double  *P                  -->  allocated memory for *all* P values         | Probably I need to                         *
 *  double  *adjP               -->  allocated memory for *all* adjusted values  | remove these initializations               *
 *  int     *pB                 -->  scalar integer value                                                                     *
 *  int     *index              -->  allocated memory for the index                                                           *
 *  char    **options           -->  c(test, side, fixed.seed.sampling)                                                       *
 *                                   ----------------------------------                                                       *
 *                                   test == c("t", "f", "blockf", "pairt", "wilcoxon", "t.equalvar")                         *
 *                                   side == c("upper","abs","lower")                                                         *
 *                                   fixed.seed.sampling == c("y","n")                                                        *
 *  int     *extra              -->  scalar integer value                                                                     *
 *  int     *generator_flag     -->  scalar integer value                                                                     *
 * ************************************************************************************************************************** */


int permutation(int n,...) {

    int worldSize, worldRank;
    int local_check = 0, global_check;
    va_list ap;

    // Variables for input arguments on master process
    double *d = NULL, *pna = NULL;
    int *pnrow = NULL, *pncol = NULL, *L = NULL, *generator_flag = NULL;
    double *T = NULL, *P = NULL, *adjP = NULL;
    int *pB = NULL, *index = NULL, *extra = NULL;
    char **options = NULL;
    int new_pB = 0, waste_cycles = 0;
    int skip_first = -1;

    // A buffer to store and send the options to slave processes
    // Also used as temp space after that
    int buf[8];

    // Variables for p-values information
    int i;
    // * **************************************************** *
    // *                   **WARNING**                        *
    // *                   -----------                        *
    // *  I allocate these vectors to be contiguous in memory *
    // *  in order to be able to perform a global reduction   *
    // *  to *all* of them in a single go                     *
    // * **************************************************** *
    int *total1 = NULL, *total2 = NULL;
    int *count1 = NULL, *count2 = NULL;

    // Get size and rank from communicator
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);

    if (worldRank == 0) {
        if (n != 13) {
            DEBUG("Rank 0 passed incorrect arguments into correlation!");
        }

        // Get input variables
        va_start(ap, n);
        d = va_arg(ap, double*);
        pnrow = va_arg(ap, int*);
        pncol = va_arg(ap, int*);
        L = va_arg(ap, int*);
        pna = va_arg(ap, double*);
        T = va_arg(ap, double*);
        P = va_arg(ap, double*);
        adjP = va_arg(ap, double*);
        pB = va_arg(ap, int*);
        index = va_arg(ap, int*);
        options = va_arg(ap, char**);
        extra = va_arg(ap, int*);
        generator_flag = va_arg(ap, int*);
        va_end(ap);

        // Send options first
        buf[0] = strlen(options[0]) + 1;
        buf[1] = strlen(options[1]) + 1;
        buf[2] = strlen(options[2]) + 1;
        buf[3] = *generator_flag;
        buf[4] = *pnrow;
        buf[5] = *pncol;
        buf[6] = *pB;
        buf[7] = *extra;

        // Broadcast options
        MPI_Bcast(buf, 8, MPI_INT, 0, MPI_COMM_WORLD);

        // Allocate vectors to store info about p-values
        alloc_init_info_vectors(buf[4], &total1, &total2, &count1, &count2);

        // Check memory allocated
        if ( total1 == NULL )
            local_check = 1;

        MPI_Allreduce(&local_check, &global_check, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        if ( global_check != 0 ) {
            // Free all allocated memory and exit
            if ( total1 != NULL ) free(total1);

            return 1;
        }

        /* ================================================================== *
         * All memory checks are done, memory is allocated in all processes,  *
         * proceed to broadcast, scatter                                      *
         * ================================================================== */

        // Broadcast the L vector
        MPI_Bcast(L, *pncol, MPI_INT, 0, MPI_COMM_WORLD);

        // Broadcast the "pna" value
        MPI_Bcast(pna, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // Broadcast string vector "options"
        MPI_Bcast(options[0], buf[0], MPI_CHAR, 0, MPI_COMM_WORLD);
        MPI_Bcast(options[1], buf[1], MPI_CHAR, 0, MPI_COMM_WORLD);
        MPI_Bcast(options[2], buf[2], MPI_CHAR, 0, MPI_COMM_WORLD);

        // Broadcast entire data set
        MPI_Bcast(d, buf[4] * buf[5], MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // Compute how many iterations the current process needs to perform
        new_pB = split_iterations(worldRank, worldSize, *pB);

        // Set flag to 0. Process with rank==0 will not skip the first
        // iteration. This process will be the *only* process to execute
        // the first permutation.
        skip_first = 0;
        waste_cycles = 0;

    } else {

        options = (char **)malloc(sizeof(char *) * 3);

        // Get some of the options (8 integer values)
        MPI_Bcast(buf, 8, MPI_INT, 0, MPI_COMM_WORLD);

        // Allocate and check *all* memory needed for arguments and results
        // Note : maxT internally will need more memory
        local_check = perm_alloc_and_check(buf, &options, &pnrow, &pncol, &pB, &extra, &L, &pna, &d, &T,
                                           &index, &generator_flag, &total1, &total2, &count1, &count2);

        // Check that all processes are OK with their memory allocations
        MPI_Allreduce(&local_check, &global_check, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        if ( global_check != 0 ) {
            perm_free(options, pnrow, pncol, pB, extra, L, pna, d, T, index, generator_flag, total1);
            return 1;
        }

        /* ================================================================== *
         * All memory checks are done, memory is allocated in all processes,  *
         * proceed to broadcasts.                                             *
         * ================================================================== */

        // Set values
        *pnrow = buf[4];   *pncol = buf[5];
        *pB = buf[6];      *extra = buf[7];
        *generator_flag = buf[3];

        // Get the L vector
        MPI_Bcast(L, *pncol, MPI_INT, 0, MPI_COMM_WORLD);

        // Get the pna value
        MPI_Bcast(pna, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // Get the "options" strings
        MPI_Bcast(options[0], buf[0], MPI_CHAR, 0, MPI_COMM_WORLD);
        MPI_Bcast(options[1], buf[1], MPI_CHAR, 0, MPI_COMM_WORLD);
        MPI_Bcast(options[2], buf[2], MPI_CHAR, 0, MPI_COMM_WORLD);

        // Get entire data set
        MPI_Bcast(d, buf[4] * buf[5], MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // * ******************************************************* *
        // *  Arguments received                                     *
        // *  Proceed to compute how many permutations each process  *
        // *  will perform and also initialize the generators to     *
        // *  correct values                                         *
        // * ******************************************************* *

        // Compute how many iterations the current process needs to perform
        new_pB = split_iterations(worldRank, worldSize, *pB);

        // If new_pB is zero then this process will not execute get_maxT so no need
        // to figure out the waste_cycles and skip_first values
        if( new_pB != 0 ) {
            // Before moving on with their permutations the processes need to
            // forward their generators to produce the exact permutations as the
            // serial version of the code. Depending on their rank, each process
            // will *waste* cycles from the random generators.
            waste_cycles = split_iterations(0, worldSize, *pB);
            waste_cycles = (waste_cycles * (worldRank-1) ) + (waste_cycles - 1);

            // The first iteration is *special* and it will be performed only
            // by the master process. The rest will skip it so we need to
            // increament the value of permutations
            new_pB++;
            skip_first = 1;
        }
    }

    /* ====================================================================== *
     *       Processes have all the data needed to start computations         *
     * ====================================================================== */

    // If permutations exist....execute
    if( new_pB != 0 ) {
        get_maxT(d, pnrow, pncol, L, pna, T, &new_pB, index, options, extra, total1, total2,
                 count1, count2, skip_first, *generator_flag, waste_cycles);
    }

    /* =================================================================== *
     *       Computations finished and the data must be gathered           *
     * =================================================================== */

    if ( worldRank == 0 ) {

        // A single reduce will sum up *all* 4 vectors because they are
        // contiguous in memory
        MPI_Reduce(MPI_IN_PLACE, total1, (*pnrow) * 4, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

        // Compute raw p-values
        for(i=0; i < *pnrow; i++) {
            if(total2[i] == 0) 
                P[i] = NA_FLOAT;
            else
                P[i] = count2[i]*1.0/total2[i];
        }

        // Compute adjusted p-values
        for(i=0; i < *pnrow; i++) {
            if(total1[i] == 0)
                adjP[i] = NA_FLOAT;
            else
                adjP[i] = count1[i]*1.0/total1[i];
        }

        // Enforce the montonicity
        for(i=1; i < *pnrow; i++)
            if( adjP[i] < adjP[i-1] )
                adjP[i] = adjP[i-1];


        // Free allocated memory
        // Note: all 4 vectors in one go
        if ( total1 != NULL )
            free(total1);

    } else {

        // A single reduce will sum up *all* 4 vectors because they are
        // contiguous in memory
        MPI_Reduce(total1, NULL, (*pnrow) * 4, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

        // Free all allocated memory
        perm_free(options, pnrow, pncol, pB, extra, L, pna, d, T, index, generator_flag, total1);
    }

    return 0;
}

/* =============================================================================== *
 *                              Helper functions                                   *
 *                              ----------------                                   *
 *      Mostly just to remove the confusing code from within the permutation       *
 *      function.                                                                  *
 * =============================================================================== */


/* ******************************************* *
 *  Splits the permutations into equal chunks  *
 * ******************************************* */
int split_iterations(int rank, int size, int pB) {

    int new_pB = -1;

    // Check if there are more processes than permtations requested
    // If it's true then assign one permutation to each processes and
    // zero to the remaining processes
    if ( pB < size ) {
        if ( rank < pB )
            return 1;
        else
            return 0;
    }

    // Compute how many iterations the current process needs to perform
    if ( (pB % size) != 0 ) {
        // Fair share
        new_pB = (pB/size) + 1;

        // Last process will get remaining iterations
        if ( rank == size-1 )
            new_pB = pB - (rank*new_pB);

        return new_pB;
    }
    else
        return (pB/size);
}

/* ********************************************************************************************** *
 *                                       perm_alloc_and_check                                     *
 *                                       --------------------                                     *
 *  Allocate all memory needed by the input arguments and also the memory needed for the vectors  *
 *  used to store information about the raw and adjusted values.                                  *
 * ********************************************************************************************** */
int perm_alloc_and_check(int *sizes, char ***options, int **pnrow, int **pncol, int **pB,
                         int **extra, int **L, double **pna,
                         double **data, double **T, int **index, int **generator_flag,
                         int **total1, int **total2, int **count1, int **count2)
{
    // First check if "options" variable is initialized ok
    if ( *options == NULL )
        return 1;

    // Allocate memory for the "options" strings
    (*options)[0] = (char *)malloc(sizeof(char) * sizes[0]);
    (*options)[1] = (char *)malloc(sizeof(char) * sizes[1]);
    (*options)[2] = (char *)malloc(sizeof(char) * sizes[2]);

    // Allocate and set values
    *pnrow = (int *)malloc(sizeof(int));
    *pncol = (int *)malloc(sizeof(int));
    *pB = (int *)malloc(sizeof(int));
    *extra = (int *)malloc(sizeof(int));
    *generator_flag = (int *)malloc(sizeof(int));

    // Allocate the L vector
    *L = (int *)malloc(sizeof(int) * sizes[5]);

    // Allocate the "pna" value
    *pna = (double *)malloc(sizeof(double));

    // Allocate memory for the new data produced
    *data = (double *)malloc(sizeof(double) * sizes[4] * sizes[5]);
    *T = (double *)malloc(sizeof(double) * sizes[4]);
    *index = (int *)malloc(sizeof(int) * sizes[4]);

    // Allocate vectors for p-values info
    alloc_init_info_vectors(sizes[4], total1, total2, count1, count2);

    // -----------------------------------
    //  Check the rest of the allocations
    // -----------------------------------

    // Check variables containing options
    if ( *pnrow == NULL || *pncol == NULL || *pna == NULL || *pB == NULL ) return 1;
    if ( *extra == NULL || *L == NULL || *generator_flag == NULL ) return 1;

    // Check data variables
    if ( *data == NULL || *index == NULL || *T == NULL) return 1;
    if ( *total1 == NULL ) return 1;

    // Success
    return 0;

}


/* ****************************************************************** *
 *                      alloc_init_info_vectors                       *
 *                      -----------------------                       *
 *  Allocate memory for the vectors that store information needed     *
 *  for computing the raw *and* adjusted p-values. (Read the warning  *
 *  notice inside the function).                                      *
 * ****************************************************************** */
void alloc_init_info_vectors(int vector_size, int **total1, int **total2, int **count1, int **count2)
{

    // * **************************************************** *
    // *                   **WARNING**                        *
    // *                   -----------                        *
    // *  I allocate these vectors to be contiguous in memory *
    // *  in order to be able to perform a global reduction   *
    // *  to *all* of them in a single go                     *
    // * **************************************************** *

    *total1 = (int *)calloc(vector_size * 4, sizeof(int));
    *total2 = *total1 + vector_size;
    *count1 = *total1 + (vector_size * 2);
    *count2 = *total1 + (vector_size * 3);
}


/* *************************** *
 *  Free all allocated memory  *
 * *************************** */
void perm_free(char **options, int *pnrow, int *pncol, int *pB, int *extra, int *L, double *pna,
               double *data, double *T, int *index, int *generator_flag, int *res)
{
    // Free variables containing options
    if ( options != NULL ) {
        if ( options[0] != NULL ) free(options[0]);
        if ( options[1] != NULL ) free(options[1]);
        if ( options[2] != NULL ) free(options[2]);
        free(options);
    }

    if ( pnrow != NULL ) free(pnrow);
    if ( pncol != NULL ) free(pncol);
    if ( pna != NULL ) free(pna);
    if ( pB != NULL ) free(pB);
    if ( extra != NULL ) free(extra);
    if ( L != NULL ) free(L);
    if ( generator_flag != NULL ) free(generator_flag);

    // Free data variables
    if ( data != NULL ) free(data); 
    if ( T != NULL ) free(T);
    if ( index != NULL ) free(index);
    if ( res != NULL ) free(res);
}


