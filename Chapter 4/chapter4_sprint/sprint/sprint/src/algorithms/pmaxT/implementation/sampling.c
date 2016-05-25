
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "string.h"
#include "mt.h"
#include "../../../sprint.h"

typedef struct tagPERMU_ARRAY {
    int n;              /* The number of original observations (samples) needs to permute       */
    int k;              /* The number of classes, labelled from 0..(k-1)                        */
                        /* which functions as the base of the integar representation            */
    int* nk;            /* The number of groups in class 0..(k-1)                               */
    int B;              /* The number of permutations samples                                   */
    int len;            /* len = floor(log(imax,k)), where imax the maximum of integers         */
    int sz;             /* The number of integars for each permutation needed sz=ceil(n/len))   */
    unsigned int *v;    /* The array, which has size of B*sz)                                   */
                        /* unsigned integers                                                    */
} PERMU_ARRAY;

static int init_permu_array(PERMU_ARRAY *pa, int *L, int n, int B);
static int get_permu(PERMU_ARRAY *pa, int h, int *L);

static int set_permu(PERMU_ARRAY* pa, int h,int *L);
static void delete_permu_array(PERMU_ARRAY* pa);

static int local_perm_count = 1;      /* The number of permutations are done  */
static int local_perm_size = 0;       /* The number of all permutations       */

static PERMU_ARRAY local_pa;

static int* local_L = NULL;

/* ******************************************************** *
 *  Create the sampling structure and initialize all values *
 * ******************************************************** */
void create_sampling(int n, int *L, int B, int generator_flag, int initial_count)
{  
    int *ordern, *permun, *myL;
    int i;

    // Set the local size of permutations
    local_perm_size = B;
    local_perm_count = 1;

    // Allocate and initialize local L to original L
    local_L = (int *)R_alloc(n, sizeof(int));

    // To check random or complete
    if( generator_flag == 1 ) {

        // Initialize the permutation array (struct)
        // Sets the B in the struct equal to 0 to be
        // able to identify the generator later on
        init_permu_array(&local_pa, L, n, 0);

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
            init_label(local_pa.n, local_pa.k, local_pa.nk, local_L);

            // Burn the cycles
            // We *must* use L in order to initialize it in the proper permutation
            // L is given as the starting position for the next permutation
            for(i=0; i < initial_count; i++) {
                next_label(local_pa.n, local_pa.k, local_pa.nk, local_L);
            }
        }

        // * ================================================================= *
        // * ================================================================= *

    } else {

        // Intiailize the permu_array (struct)
        init_permu_array(&local_pa, L, n, B);

        permun = (int *)R_alloc(local_pa.n, sizeof(int));
        ordern = (int*)R_alloc(local_pa.n, sizeof(int));
        myL = (int *)R_alloc(local_pa.n, sizeof(int));

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
        // ("permun" is a safe scratch space for this)
        for(i=0; i < initial_count; i++) {
            sample(permun, n);
        }

        // * ================================================================= *
        // * ================================================================= *

        for(i=0; i<n; i++){
            ordern[i]=i;
        }

        // Allocate and assign the values for l_first_sample
        set_permu(&local_pa, 0, L);

        for(i=1; i<B; i++) {
            memcpy(permun, ordern, sizeof(int)*n);
            sample(permun, n);

            // Change to labbeling
            sample2label(n, local_pa.k, local_pa.nk, permun, myL);
            set_permu(&local_pa, i, myL);
        }


    }
}


/* ******************************* *
 *  Pop out the first permutation  *
 * ******************************* */
int first_sample(int *L)
{
    if(L == NULL)
        return local_perm_size;

    // If is random, we'll choose the original L as the first sample
    // Note: local_pa.B == 0 => complete permutations
    //       local_pa.B >  0 => random permutations
    if(local_pa.B > 0) {
        get_permu(&local_pa, 0, L);
    } else {
        init_label(local_pa.n, local_pa.k, local_pa.nk, L);
        memcpy(local_L, L, sizeof(int) * local_pa.n);
    }

    return 1;
}


/* ******************************* *
 *  Generate the next permutation  *
 * ******************************* */
int next_sample(int *L)
{
    if(local_perm_count >= local_perm_size)
        return 0;

    // Note: local_pa.B == 0 => complete permutations
    //       local_pa.B >  0 => random permutations
    if(local_pa.B > 0) {
        get_permu(&local_pa, local_perm_count, L);
    } else {
        next_label(local_pa.n, local_pa.k, local_pa.nk, local_L);
        memcpy(L, local_L, sizeof(int) * local_pa.n);
    }

    local_perm_count++;

    return 1;
}


/* ******************************************* *
 *  Initialize the permutation array (struct)  *
 * ******************************************* */
static int init_permu_array(PERMU_ARRAY *pa, int *L, int n, int B)
{
    int i;
    unsigned imax;
    pa->n=n;
    pa->B=B;
    pa->nk=NULL;
    pa->v=NULL;

    // Compute the k (number of different classes)
    pa->k=0;
    for(i=0; i<n; i++)
        if(L[i] > pa->k)
            pa->k=L[i];
    (pa->k)++;

    // Compute nk (how many rows each class has)
    pa->nk = (int *)R_alloc(pa->k, sizeof(int));
    memset(pa->nk, 0, sizeof(int)*pa->k);
    for(i=0; i<n; i++)
        pa->nk[L[i]]++;

    // Computer imax, len
    // Get all bits are 1 for the integars
    imax=~0;
    pa->len = floor(log(imax+1.0)/log(pa->k)); 
    pa->sz = ceil(n/(pa->len*1.0));

    // Allocate the space for v
    pa->v = (unsigned int*)R_alloc(B*pa->sz, sizeof(int));

    return 1;
}


/* ********************************************************** *
 *  Pop next permutation from the permutation array (struct)  *
 * ********************************************************** */
static int get_permu(PERMU_ARRAY *pa, int h, int *L)
{
    int i,j;
    unsigned val;

    memset(L, 0, sizeof(unsigned int)*pa->n);

    if((h+1) > pa->B)
        return 0;

    for(j=0; j<pa->sz; j++) {

        // Starting from the last bit
        i=j*pa->len; 
        val=pa->v[h*pa->sz+j];

        while(val>0) {
            // This code maybe faster if necessary
            L[i]=val%(unsigned int)(pa->k);
            i++;

            // To move another bit
            val/=(unsigned int)(pa->k);
        }
    }

    return 1;
}


/* *********************************************** *
 *  Add permutation to permutation array (struct)  *
 * *********************************************** */
static int set_permu(PERMU_ARRAY *pa, int h, int *L)
{
    int i, j, nextbound;
    unsigned val, pow;

    if((h+1) > pa->B)
        return 0;

    // Starting from the last bit
    i=0;
    for(j=0; j<pa->sz; j++) {
        nextbound=(j+1)*pa->len;

        if(nextbound > (pa->n))
            nextbound=pa->n;

        pow=1;
        val=0;

        while(i < nextbound) {
            val+=(unsigned int)(L[i])*pow;
            pow*=(unsigned int)pa->k;
            i++;
        }
        pa->v[h*pa->sz+j]=val;
    }
    return 1;
}


/* ************************************************ *
 *  Free all memory in permutations array (struct)  *
 * ************************************************ */
void delete_sampling(void)
{
    delete_permu_array(&local_pa);

    if ( local_L != NULL ) {
        local_L = NULL;
    }
}


/* *********************** *
 *  Free allocated memory  *
 * *********************** */
static void delete_permu_array(PERMU_ARRAY *pa)
{
    if(pa->nk != NULL ) {
        pa->nk = NULL;
    }

    if(pa->B != 0) {
        pa->v = NULL;
    }
}


