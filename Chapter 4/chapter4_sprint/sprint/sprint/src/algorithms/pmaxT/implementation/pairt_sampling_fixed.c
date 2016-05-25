
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "mt.h"
#include "../../../sprint.h"

static int local_n=0;               /* The number of samples for permutations   */
static int local_perm_size=0;       /* The number of total simultaions          */
static int local_perm_count=1;      /* The number of permutations are done      */
static int* local_L=NULL;


/* ******************************************************** *
 *  Create the sampling structure and initialize all values *
 * ******************************************************** */
void create_sampling_pairt_fixed(int n, int *L, int B, int generator_flag, int initial_count)
{
    int i;

    local_n = n;
    local_perm_size = B;
    local_perm_count = 1;

    local_L = (int*)R_alloc(n,sizeof(int));

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

    // Copy original L in local_L
    memcpy(local_L, L, sizeof(int)*n);
}


/* ******************************* *
 *  Pop out the first permutation  *
 * ******************************* */
int first_sample_pairt_fixed(int *L)
{
    if(L == NULL)
        return local_perm_size;
    else{
        memcpy(L, local_L, sizeof(int)*local_n);
    }

    return 1;  
}


/* ******************************* *
 *  Generate the next permutation  *
 * ******************************* */
int next_sample_pairt_fixed(int *L)
{ 
    int n=local_n, i;
    double tmp;

    if(local_perm_count >= local_perm_size)
        return 0;

    for(i=0; i<n; i++) {
        tmp = get_rand();
        if(tmp>0.5)
            L[i]=1;
        else
            L[i]=0;       
    }
    
    // Increament the counter
    local_perm_count++;

    return 1;
}


/* ************************************************ *
 *  Free all memory in permutations array (struct)  *
 * ************************************************ */
void delete_sampling_pairt_fixed(void)
{
    if ( local_L != NULL ) {
        local_L=NULL;
    }
}


