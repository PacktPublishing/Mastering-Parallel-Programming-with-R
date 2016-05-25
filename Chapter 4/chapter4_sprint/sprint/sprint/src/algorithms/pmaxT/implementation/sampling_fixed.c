
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "string.h"
#include "mt.h"
#include "../../../sprint.h"

static int local_n=0;           /* The length of L                      */
static int local_k=0;           /* The number of groups                 */
static int* local_nk=NULL;      /* The number of objects in each groups */
static int* local_L=NULL;       /* The storrage of first label          */
static int local_b=1;           /* The number of permutations are done  */
static int local_B=0;           /* The number of all permutations       */
static int* local_permun=NULL;
static int* local_ordern=NULL;


/* ******************************************************** *
 *  Create the sampling structure and initialize all values *
 * ******************************************************** */
void create_sampling_fixed(int n, int *L, int B, int generator_flag, int initial_count)
{  
    int i,k;

    local_n = n;
    local_B = B;
    local_b = 1;

    // Allocate memory and copy class labels there
    local_L = (int *)R_alloc(n, sizeof(int));
    memcpy(local_L, L, sizeof(int)*n);

    // Count how many classes the set has
    k=0;
    for(i=0; i<n; i++)
        if(L[i]>k)
            k=L[i];
    k++;
    // Set the value of local_k to the class count
    local_k=k;

    local_nk = (int *)R_alloc(k, sizeof(int));
    memset(local_nk, 0, sizeof(int)*k);

    for(i=0; i<n; i++)
        local_nk[L[i]]++;

    local_permun = (int *)R_alloc(n, sizeof(int));
    local_ordern = (int *)R_alloc(n, sizeof(int));

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
    // ("local_permun" is a safe scratch space for this)
    for(i=0; i < initial_count; i++) {
        sample(local_permun, n);
    }

    // * ================================================================= *
    // * ================================================================= *

    for(i=0; i<n; i++){
        local_ordern[i]=i;
    }
}


/* ******************************* *
 *  Pop out the first permutation  *
 * ******************************* */
int first_sample_fixed(int *L)
{
    if(L == NULL)
        return local_B;
    else
        memcpy(L, local_L, sizeof(int)*local_n);

    return 1;
}


/* ******************************* *
 *  Generate the next permutation  *
 * ******************************* */
int next_sample_fixed(int *L)
{
    int n=local_n;

    if(local_b >= local_B)
        return 0;

    memcpy(local_permun, local_ordern, sizeof(int)*n);

    // Get next permutation
    sample(local_permun, n);

    // Change to labbeling
    sample2label(n, local_k, local_nk, local_permun, L);

    // Increase permutation count
    local_b++;

    return 1;
}


/* ************************************************ *
 *  Free all memory in permutations array (struct)  *
 * ************************************************ */
void delete_sampling_fixed(void)
{


    // Set pointers to NULL
    local_L=NULL;
    local_nk=NULL;
    local_permun=NULL;
    local_ordern=NULL;
}

