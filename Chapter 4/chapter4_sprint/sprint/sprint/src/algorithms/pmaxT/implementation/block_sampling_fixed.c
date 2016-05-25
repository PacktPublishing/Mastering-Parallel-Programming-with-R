
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "string.h"
#include "mt.h"
#include "../../../sprint.h"

static int local_n=0;
static int local_perm_size=0;           /* The number of total simultaions      */
static int local_perm_count=1;          /* The number of permutations are done  */
static int local_random_flag=1;         /* The permuation is random or not      */
static int* local_L=NULL;
static int local_m=0;                   /* The number of treaments              */
static int* l_order_block=NULL;


/* ******************************************************** *
 *                  create_sampling_block                   *
 *                  ---------------------                   *
 *  Create the sampling structure and initialize all values *
 * ******************************************************** */
void create_sampling_block(int n, int *L, int B, int generator_flag, int initial_count)
{
    int i, m;

    m=0;
    for(i=0; i<n; i++)
        if(L[i]>m) {
            m++;
        }
    m++;

    local_L = (int *)R_alloc(n, sizeof(int));
    local_perm_size = B;
    local_perm_count = 1;
    local_n = n;
    local_m = m;

    if(generator_flag == 1) {

        // Doing complete permutation
        local_random_flag=0;

        // * ================================================================= *
        // *                     Forwarding logic                              *
        // *                     ----------------                              *
        // *  For parallel executions all processes, appart from the one with  *
        // *  rank == 0, should forward their random generators in order to    *
        // *  be able to reproduce the exact same random permutations as the   *
        // *  serial version of the code (or a version with only one thread)   *
        // *                                                                   *
        // *  The code below uses the initial_count variable to "burn" the     *
        // *  cycles from the random generator                                 *
        // * ================================================================= *

        // All processes apart from process with rank == 0 (translates to
        // initial_count == 0) will perform this forward
        if ( initial_count != 0 ) {

            // Initialize label
            init_label_block(local_L, local_n, local_m);

            // Burn the cycles
            // We *must* use L in order to initialize it in the proper permutation
            // L is given as the starting position for the next permutation
            for(i=0; i < initial_count; i++) {
                next_label_block(local_L, local_n, local_m);
            }
        }

    } else {
        // Doing random permutation
        local_random_flag=1;

        l_order_block = (int *)R_alloc(n, sizeof(int));
        init_label_block(l_order_block, n, m);


        // * ================================================================= *
        // *                     Forwarding logic                              *
        // *                     ----------------                              *
        // *  For parallel executions all processes, appart from the one with  *
        // *  rank == 0, should forward their random generators in order to    *
        // *  be able to reproduce the exact same random permutations as the   *
        // *  serial version of the code (or a version with only one thread)   *
        // *                                                                   *
        // *  The code below uses the initial_count variable to "burn" the     *
        // *  cycles from the random generator                                 *
        // * ================================================================= *

        // Set initial seed
        set_seed(g_random_seed);

        // Burn the cycles
        // ("local_L" is a safe scratch space for this)
        for(i=0; i < initial_count; i++) {
            sample(local_L, n);
        }

        // * ================================================================= *
        // * ================================================================= *

        // Initialize local_L with original L
        memcpy(local_L, L, sizeof(int)*n);
    }
}


/* ***************************************************************** *
 *                       first_sample_block                          *
 *                       ------------------                          *
 *  It initializes the first permutation vector depending on the     *
 *  user's options. For "random" permutations the first permutation  *
 *  is the input L vector. For "complete" permutations the function  *
 *  "init_label_block" is called to initialize the vector.           *
 * ***************************************************************** */
int first_sample_block(int *L)
{
    if(L == NULL)
        return local_perm_size;

    if(local_random_flag) {
        memcpy(L, local_L, sizeof(int)*local_n);
    } else {
        init_label_block(local_L, local_n, local_m);
        memcpy(L, local_L, sizeof(int) * local_n);
    }

    return 1;
}


/* ********************************************************************** *
 *                        next_sample_block                               *
 *                        -----------------                               *
 *  Get the next permutation. There are two permutation generators        *
 *  for doing this. The "random" generator and the "complete"             *
 *  permutations generator. Depending on the user options the             *
 *  function will choose between these two using the "local_random_flag"  *
 *  variable value.                                                       *
 * ********************************************************************** */
int next_sample_block(int *L)
{
    if(local_perm_count >= local_perm_size)
        return 0;

    if(local_random_flag){
        memcpy(L, l_order_block, sizeof(int)*local_n);
        sample_block(L, local_n, local_m);
    } else {
        next_label_block(local_L, local_n, local_m);
        memcpy(L, local_L, sizeof(int) * local_n);
    }

    local_perm_count++;

    return 1;
}


/* ****************************************************** *
 *                    sample_block                        *
 *                    ------------                        *
 *  Using the random generator function "sample" it will  *
 *  generate and return the next random permutation       *
 * ****************************************************** */
void sample_block(int *L, int n, int m)
{
    int block, b, s;

    s=0;
    block=n/m;

    for(b=0; b<block; b++) {
        sample(L+s, m);
        s+=m;
    }
}

/* **************************************************** *
 *               delete_sampling_block                  *
 *               ---------------------                  *
 *  Free dynamically allocated memory reserved for the  *
 *  permutation vectors and set pointer to NULL         *
 * **************************************************** */
void delete_sampling_block(void)
{
    local_L = NULL;

    if ( l_order_block != NULL ) {
        l_order_block = NULL;
    }
}

