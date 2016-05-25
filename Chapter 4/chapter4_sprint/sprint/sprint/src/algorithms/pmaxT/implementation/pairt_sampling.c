
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "string.h"
#include "mt.h"
#include "../../../sprint.h"

static int local_n=0;                           /* The number of samples for permutations                            */
static int local_perm_size=0;                   /* The number of total simultaions                                   */
static int local_perm_count=1;                  /* The number of permutations are done                               */
static int local_current_permutation=1;         /* The number of the current int to be converted into a permutation  */
static int local_random_flag=-1;                /* The permuation is random or not                                   */
static unsigned int* l_all_samples=NULL;

// Store all the samples in random case
static int local_size=0;                        /* The number of bytes for per permutation  */
static int local_length=0;
static int get_binpermu(int h, int n, int sz, int len, int *L, int hMax, unsigned int *V);
static int set_binpermu(int *L, int h, int n, int sz, int len, int hMax, unsigned int *V);


/* ******************************************************** *
 *  Create the sampling structure and initialize all values *
 * ******************************************************** */
void create_sampling_pairt(int n, int *L, int B, int generator_flag, int initial_count)
{
    int i, j;
    unsigned int imax;
    int *myL;
    double tmp;

    local_n=n;
    local_perm_count=1;

    imax=(unsigned int)(~0);
    local_length=floor(log(imax+1.0)/log(2));
    local_size=ceil(n/(local_length*1.0));

    local_perm_size=B;

    if(generator_flag == 1) {

        local_random_flag = 0;
        local_current_permutation = 1;

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
            // Burn the cycles
            // Set the number of current permutation int to the number
            // we want to start at this process
            local_current_permutation = initial_count + 1;
        }

        // * ================================================================= *
        // * ================================================================= *


    } else {

        myL = (int *)R_alloc(n, sizeof(int));
        local_random_flag=1;

        l_all_samples = (unsigned int *)R_alloc(local_perm_size*local_size, sizeof(int));

        // Setting the first sample as the original data
        set_binpermu(L, 0, n,local_size, local_length, local_perm_size, l_all_samples);

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
        // ("myL" is a safe scratch space for this)
        for(i=0; i < initial_count; i++) {
            sample(myL, n);
        }

        // * ================================================================= *
        // * ================================================================= *

        // The extra as a buffer
        for(i=1; i<local_perm_size; i++) {

            for(j=0; j<n; j++) {
                tmp=get_rand();
                if(tmp>0.5)
                    myL[j]=1;
                else
                    myL[j]=0; 
            }
            set_binpermu(myL, i, n, local_size, local_length, local_perm_size, l_all_samples);
        }
    }
}


/* ******************************* *
 *  Pop out the first permutation  *
 * ******************************* */
int first_sample_pairt(int *L)
{
    // Return the number of real samples
    if(L==NULL)
        return local_perm_size;

    // Call different sampling function
    if(local_random_flag){
        get_binpermu(0, local_n, local_size, local_length, L, local_perm_size, l_all_samples);
    }
    else
        int2bin(0, L, local_n);

    return 1;
}


/* ******************************* *
 *  Generate the next permutation  *
 * ******************************* */
int next_sample_pairt(int *L)
{ 
    if(local_perm_count >= local_perm_size)
        return 0;

    // Call different sampling function
    if(local_random_flag)
        get_binpermu(local_perm_count, local_n, local_size, local_length, L, local_perm_size, l_all_samples);
    else
        int2bin(local_current_permutation, L, local_n);

    // We need two counters. One to count how many
    // permutations were executed and one to choose which
    // permutation to pop out in case of complete permutations
    local_perm_count++;
    local_current_permutation++;

    return 1; 
}


/* ********************************************************** *
 *  Pop next permutation from the permutation array (struct)  *
 * ********************************************************** */
static int get_binpermu(int h, int n, int sz, int len, int *L, int hMax, unsigned int *V)
{
    int i, j;
    unsigned val;

    memset(L, 0, sizeof(unsigned int)*n);

    if((h+1) > hMax)
        return 0;
    for(j=0; j<sz; j++) {
        // Starting from the last bit
        i=j*len; 
        val=V[h*sz+j];

        while(val>0) {
            // This code maybe faster if necessary
            L[i]=val&1;
            i++;
            // To move another bit
            val>>=1;
        }
    }

    return 1;
}


/* *********************************************** *
 *  Add permutation to permutation array (struct)  *
 * *********************************************** */
static int set_binpermu(int *L, int h, int n, int sz, int len, int hMax, unsigned int *V)
{

    int i,j,nextbound;
    unsigned val,pow;

    if((h+1) > hMax)
        return 0;

    // Starting from the last bit
    i=0; 
    for(j=0; j<sz; j++) {
        nextbound=(j+1)*len;

        if(nextbound>n)
            nextbound=n;
        pow=1;
        val=0;

        while(i<nextbound) {
            val+=(unsigned int)(L[i])*pow;
            pow<<=1;
            i++;
        }
        V[h*sz+j]=val;
    }

    return 1;
}


/* ************************************************ *
 *  Free all memory in permutations array (struct)  *
 * ************************************************ */
void delete_sampling_pairt(void)
{
    if(local_random_flag) {
        if(local_perm_size != 0) {
            l_all_samples=NULL;
        }
    }
}

