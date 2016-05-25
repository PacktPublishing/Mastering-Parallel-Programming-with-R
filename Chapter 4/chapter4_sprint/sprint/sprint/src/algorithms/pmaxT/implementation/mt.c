/* ********************************* *
 *            Header files           *
 * ********************************* */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <ctype.h>
#include <stdarg.h>
#include "mt.h"
#include "../../../sprint.h"

/* ************************************************** *
 *                malloc_gene_data                    *
 *                ----------------                    *
 *  Allocate the necessary space for the big data,    *
 *  see the comments about the structrue GENE_DATA    *
 * ************************************************** */
int myDEBUG = 0;
long int g_random_seed = 3455660;

void malloc_gene_data(GENE_DATA *pdata)
{
    int i;
    int nrow=pdata->nrow;
    int ncol=pdata->ncol;

    pdata->id = (char **)R_alloc(nrow, sizeof(char*));
    pdata->d = (double **)R_alloc(nrow, sizeof(double*));
    pdata->L = (int *)R_alloc(ncol, sizeof(int));

    /*initialization*/
    memset(pdata->L, 0, sizeof(int)*ncol);
    for(i=0; i<ncol; i++) 
        pdata->L[i] = 0;

    for (i=0; i<nrow; i++) {
        pdata->id[i] = (char *)R_alloc(MAX_ID, sizeof(char));
        pdata->d[i] = (double *)R_alloc(ncol, sizeof(double));
    }
}

/* ************************************ *
 *           free_gene_data             *
 *           --------------             *
 *  Free the space allocated for pdata  *
 * ************************************ */
void free_gene_data(GENE_DATA *pdata)
{
    int i;

}

/* **************************************************************************** *
 *                           compute_test_stat                                  *
 *                           -----------------                                  *
 *  L is an array which contains 0,1,2 ... for specifying class label           *
 *  T is the test_test needs to return                                          *
 *  func_stat is a functin pointer with the following protocol                  *
 *      double func_stat(double *Y, int* L,int n, double na,const void* extra)  *
 * **************************************************************************** */
void compute_test_stat(GENE_DATA *pdata, int *L, double *T,
                       FUNC_STAT func_stat, const void *extra)
{
    int i;

    for(i=0; i<pdata->nrow; i++)
        T[i]=(*func_stat)(pdata->d[i], L, pdata->ncol, pdata->na, extra);
}
 
/* ********************************************************************************************* *
 *                                      get1pvalue                                               *
 *  We'll do complete resampling with the function next_sample, which is determined by           *
 *  the function next_sample when it returns false.                                              *
 *  in the next_sample, you need to decide to choose complet resampling or not or                *
 *  whatever you like                                                                            *
 *                                                                                               *
 *  L: is the labelling of each experiment                                                       *
 *  T: is that test statistics                                                                   *
 *  P: unadjtesed P-values                                                                       *
 *                                                                                               *
 *  To use the function first_sample, and next_sample, they're needed to write into a separate   *
 *  file, where it provides the create_sampling to do initialization such as allocate the space  *
 *  (before use the sampling) and delete_sampling after we've done the sampling in the main().   *
 *                                                                                               *
 *  int first_sample(int *L)                                                                     *
 *      get the first sample of the labelling.                                                   *
 *  int next_sample(int* L)                                                                      *
 *      get the next sample, if it's done all the sampling, then it returns 0,                   *
 *      otherwise it returns 1.                                                                  *
 *                                                                                               *
 *  input: pdata, L, B,next_sample,func_stat                                                     *
 *  output: T,P                                                                                  *
 *                                                                                               *
 * ********************************************************************************************* */
void  get1pvalue(GENE_DATA *pdata, int *L, double *T, double *P,
                 FUNC_STAT func_stat, FUNC_SAMPLE func_first_sample,
                 FUNC_SAMPLE func_next_sample, FUNC_CMP func_cmp, const void *extra)
{
    int b=0,*bL,i,is_next,*total;
    double *bT, *count;
    int ncol=pdata->ncol;
    int nrow=pdata->nrow;

    // Allocate the space and initialziation
    bT = (double*)R_alloc(nrow, sizeof(double));
    bL = (int*)R_alloc(ncol, sizeof(int));
    count = (double*)R_alloc(nrow, sizeof(double));
    memset(count, 0, sizeof(double)*nrow); 
    total = (int*)R_alloc(nrow, sizeof(int));
    memset(total, 0, sizeof(int)*nrow);

    // Comuter the original one first
    compute_test_stat(pdata, L, T, func_stat,extra);

    // Iteration for permutaion
    (*func_first_sample)(bL);
    is_next=1;
    b=0;
    while(is_next) {
        compute_test_stat(pdata, bL, bT, func_stat, extra);
        for(i=0; i<nrow; i++) {
            if(bT[i] == NA_FLOAT)
                continue;
            if(T[i] == NA_FLOAT)
                continue;
            // Right now I only implements the 3 cases, which are pretty common
            if((func_cmp==cmp_high) && (bT[i] >= T[i]-EPSILON)){
                count[i]+=1;
            } else if((func_cmp==cmp_low) && (bT[i] <= T[i]+EPSILON)){
                count[i]+=1;
            } else if ((func_cmp==cmp_abs) && (fabs(bT[i]) >= fabs(T[i])-EPSILON)){
                count[i]+=1;
            }	      
            total[i]++;
        }
        b++;

        is_next=(*func_next_sample)(bL);
    }

    // Summarize the results
    for(i=0; i<nrow; i++){ 
        if(total[i]==0) 
            P[i] = NA_FLOAT;
        else
            P[i] = count[i]*1.0/total[i];
    }

    // Free the spaces

}


/* ********************************************************************************** *
 *                                 sort_gene_data                                     *
 *                                 --------------                                     *
 *  Description:                                                                      *
 *   sort the rows of gene_data such that row R[i] of is the first row, i=0,...,m-1,  *
 *   wher R[0],...,R[m-1] is a permutation of (0,...,m-1)                             *
 * ********************************************************************************** */										 
void sort_gene_data(GENE_DATA *pdata,int *R)
{
    int i, nrow=pdata->nrow;
    char **old_id;      /* the old addresses of the gene id */
    double **old_d;     /* th old addresses of the gene data */

    old_d = (double **)R_alloc(nrow, sizeof(double*));
    old_id = (char **)R_alloc(nrow, sizeof(char*));

    /* store the original pointers from pdata */
    for(i=0; i<nrow; i++) {
        old_d[i] = pdata->d[i];
        old_id[i] = pdata->id[i];
    }

    /*rearrange the data so that it's ordered according to R*/
    for(i=0; i<nrow; i++) {
        pdata->d[i] = old_d[R[i]];
        pdata->id[i] = old_id[R[i]];
    }


}


/* ************************************************************* *
 *                        sort_vector                            *
 *                        -----------                            *
 *  Desciption                                                   *
 *   sort the vector V according to the order R with n elemnets  *
 *   where R[0],...,R[n-1] is a permutation of 0,...n-1          *
 * ************************************************************* */
void sort_vector(double *V, int *R,int n)
{
    double *old_V;
    int i;

    old_V = (double *)R_alloc(n, sizeof(double));

    for(i=0; i<n; i++)
        old_V[i] = V[i];
    for(i=0; i<n; i++)
        V[i] = old_V[R[i]];


}


/* **************************************************************************** *
 *                              get_all_samples_P                               *
 *                              -----------------                               *
 *  Descriptions: Try to get all the unadjusted p-values for a gene with        *
 *  experssion values at V with n experiemtns.                                  *
 *                                                                              *
 *  int first_sample(int *L)                                                    *
 *  get the first sample of the labelling.                                      *
 *  if L==NULL, then it returns all the possible simulations, which depends on  *
 *  the initial function create_sampling.                                       *
 *  int next_sample(int* L)                                                     *
 *  get the next sample, if it's done all the sampling, then it returns 0,      *
 *  otherwise it returns 1.                                                     *
 *  output is P                                                                 *
 * **************************************************************************** */
void get_all_samples_P(double *V, int n, double *P, double na, 
                       FUNC_STAT func_stat, FUNC_SAMPLE func_first_sample, 
                       FUNC_SAMPLE func_next_sample, FUNC_CMP func_cmp,
                       const void *extra)
{
    int  *L,*R,i,oldb,is_next,b=0,B_new,B;
    double* T=P,oldf;

    B=(*func_first_sample)(NULL);

    // Allocate the spaces
    L = (int*)R_alloc(n, sizeof(int));
    R = (int*)R_alloc(B, sizeof(int));

    // Compute all the test_stat
    (*func_first_sample)(L);
    is_next=1;
    B_new=0;

    while(is_next) {
        T[b] = func_stat(V, L, n, na, extra);
        if(T[b] != NA_FLOAT)
            B_new++;
        b++;
        is_next=(*func_next_sample)(L);
    }

    if(B!=b){
      error("Error we have b(%d)!=B(%d)\n",b,B);
      return;
    }

    // Order the test_stat
    order_data(T,R,B,func_cmp);

    // Note the last elements of B-B_new has NA T-value
    // Assign the probabilites
    oldb=0;
    oldf=T[R[0]];
    for(b=1; b<B_new; b++){
        if((func_cmp == cmp_high) & (T[R[b]] >= oldf-EPSILON))
            continue;
        else if ((func_cmp == cmp_low ) && (T[R[b]] <= oldf+EPSILON))
            continue;
        else if((func_cmp == cmp_abs )&& fabs(T[R[b]]) >= fabs(oldf)-EPSILON)
            continue;

        for(i=oldb; i<b; i++)
            P[R[i]]=(b+0.0)/B_new;
        oldb=b;
        if(b<B_new-1)
            oldf=T[R[b]];
    }

    for(i=oldb; i<b; i++)
        P[R[i]]=1.0;

    // For NA test_stat, assign NA probabilites
    for(b=B_new;b<B;b++)
        P[R[b]]=NA_FLOAT;

} 

/*get all the samples of T and they're also ordered.It's used only for diagonsis*/
void get_all_samples_T(double* V, int n,double* T,double na, 
        FUNC_STAT func_stat,FUNC_SAMPLE func_first_sample, 
        FUNC_SAMPLE func_next_sample,const void* extra)
{
    int  *L,*R,is_next,b=0,B;

    B = (*func_first_sample)(NULL);
    /*allocate the spaces*/
    L = (int*)R_alloc(n, sizeof(int));
    R = (int*)R_alloc(B, sizeof(int));

    /*compute all the test_stat*/
    (*func_first_sample)(L);
    is_next=1;

    while(is_next){
        T[b]=func_stat(V,L,n,na,extra);
        b++;
        is_next=(*func_next_sample)(L);
    }

    if(B!=b){
      error("Error we have b(%d)!=B(%d)\n",b,B);
        return;
    }


} 

void adj_pvalue_quick(GENE_DATA *pdata, double *T, double *P,
        double *Adj_P, double* Adj_Lower,
        FUNC_STAT func_stat, FUNC_STAT func_stat_T,
        FUNC_SAMPLE func_first_sample,
        FUNC_SAMPLE func_next_sample, FUNC_CMP func_cmp, const void *extra)
{

    int *L, b, B, B_new, i, *R, neq; /*b for simulation*, neq is for the number of equal signs*/
    double *all_P, *all_Q, count;
    int ncol=pdata->ncol, nrow=pdata->nrow;

    // Allocate the space*/
    B=(*func_first_sample)(NULL);
    L = (int*)R_alloc(ncol, sizeof(int)); 
    R = (int*)R_alloc(nrow, sizeof(int));
    all_P = (double*)R_alloc(B, sizeof(double));
    all_Q = (double*)R_alloc(B, sizeof(double));

    // Get the original unadjusted p-values first
    // we'll use the normalized t-statistics
    get1pvalue(pdata,pdata->L,T,P,func_stat_T,func_first_sample,func_next_sample,func_cmp,extra);

    // Sort the test_stat
    order_mult_data(R, nrow, 2, P, cmp_low, T, func_cmp);
    /*order_data(P,R,nrow,func_cmp);*/

    // Rearrange the data according the unadjusted p-values
    sort_gene_data(pdata, R);
    sort_vector(T, R, nrow);
    sort_vector(P, R, nrow);

    // Initialze all_Q[]=NA_FLOAT
    for(b=0;b<B;b++)
        all_Q[b]=NA_FLOAT;

    // Loop for each gene
    for(i=nrow-1;i>=0;i--){
        get_all_samples_P(pdata->d[i], ncol, all_P,pdata->na,
                func_stat, func_first_sample, func_next_sample, func_cmp, extra);

        // Update all_Q
        count=0;
        B_new=0;
        neq=0;

        for(b=0;b<B;b++){
            if (all_P[b] == NA_FLOAT)
                break;/*we don't need care about NA pvlaues*/
            // Update q* by the value p
            if(all_Q[b] > all_P[b])
                all_Q[b]=all_P[b];
            // Skip NA
            if(all_Q[b] == NA_FLOAT)
                continue;
            if(all_Q[b] < P[i]){
                count+=1;
            } else if (all_Q[b] <= P[i]+EPSILON)
                neq++;
            B_new++;
        }

        // Assign the Adj_P and Adj_Lower for gene i
        if(B_new!=0) {
            Adj_P[i]=(count+neq)/B_new;

            if(neq==0) 
                Adj_Lower[i]=count/B_new;
            else Adj_Lower[i]=(count+1)/B_new;
        }
        else {
            Adj_P[i]=NA_FLOAT;
            Adj_Lower[i]=NA_FLOAT; 
        }
    }

    /* to make monotone of Adj_P and Adj_Lower*/
    for(i=1;i<nrow;i++)
        if(Adj_P[i]<Adj_P[i-1])
            Adj_P[i]=Adj_P[i-1];

    for(i=1; i<nrow; i++)
        if(Adj_Lower[i] < Adj_Lower[i-1])
            Adj_Lower[i] = Adj_Lower[i-1];

    /*free the spaces*/

}
	
/* ************************************************************************************* *
 *                               adj-by_t                                                *
 *                               --------                                                *
 *  We'll do complete resampling with the function next_sample, which is determined by   *
 *  the function next_sample when it returns false.                                      *
 *   in the next_sample, you need to decide to choose complet resampling or not or       *
 *   whatever you like                                                                   *
 *                                                                                       *
 *  L: is the labelling of each experiment                                               *
 *  T: is that test statistics                                                           *
 *  P: unadjtesed P-values                                                               *
 *  Adj_P:ajusted p-values by using the max|T|                                           *
 *  To use the function first_sample, and next_sample, they're needed to write into a    *
 *  separate file, where it provides the create_sampling to do initialization such as    *
 *  allocate the space (before use the sampling) and delete_sampling after we've done    *
 *  the sampling in the main().                                                          *
 *                                                                                       *
 *  int first_sample(int *L) :  get the first sample of the labelling.                   *
 *  int next_sample(int* L)  :  get the next sample, if it's done all the sampling,      *
 *                              then it returns 0, otherwise it returns 1.               *
 *                                                                                       *
 *  input: pdata, L, B,next_sample,func_stat                                             *
 *  output: T,P                                                                          *
 *                                                                                       *
 * ************************************************************************************* */

void  adj_by_T(GENE_DATA *pdata, double *T, FUNC_STAT func_stat, FUNC_SAMPLE func_first_sample, 
               FUNC_SAMPLE func_next_sample, FUNC_CMP func_cmp, const void *extra,
               int *total1, int *total2, int *count1, int *count2, int skip_first)
{
    int b=0, *bL, i, is_next, *R;
    // qT is the successiv maxima
    double *bT, qT;
    int ncol=pdata->ncol;
    int nrow=pdata->nrow;

    // Allocate the space and initialziation
    bT = (double*)R_alloc(nrow, sizeof(double));
    bL = (int*)R_alloc(ncol, sizeof(int));
    R = (int*)R_alloc(nrow, sizeof(int));

    // Compute the original t-stat first
    compute_test_stat(pdata, pdata->L, T, func_stat, extra);

    // Sort the T
    order_data(T, R, nrow, func_cmp);
    sort_gene_data(pdata, R);
    sort_vector(T, R, nrow);

    // In parallel runs only the master process needs to perform
    // the first iteration. The rest of the processes *must* skip
    // to the next permutation
    if( skip_first )
        (*func_next_sample)(bL);
    else
        (*func_first_sample)(bL);

    // Changed to the orignal stat, which is monotone of t and centered
    is_next=1;
    b=0;
    while(is_next) {

        compute_test_stat(pdata, bL, bT, func_stat, extra);

        // Deal with unajdused value first
        for(i=0; i<nrow; i++){

            if(T[i]==NA_FLOAT)
                continue;

            if( bT[i] != NA_FLOAT ) {
                if( (func_cmp == cmp_high) && (bT[i]+EPSILON >= T[i]) )
                    count2[i]++;
                if( (func_cmp == cmp_low) && (bT[i] <= T[i]+EPSILON))
                    count2[i]++;
                if( (func_cmp == cmp_abs) && (fabs(bT[i]) >= fabs(T[i])-EPSILON))
                    count2[i]++;

                total2[i]++;
            }
        }

        // Deal with adjusted values
        // Intitalize the qT
        qT=NA_FLOAT;

        // Looping the row reversely
        for(i=nrow-1; i >= 0; i--) {

            if( T[i] == NA_FLOAT )
                continue;

            // Right now I only implements the 3 cases, which are pretty common
            if( func_cmp == cmp_high ) {

                if( (bT[i] != NA_FLOAT) && (qT != NA_FLOAT) && (bT[i] > qT) )
                    qT=bT[i];
                if( (bT[i] != NA_FLOAT) && (qT == NA_FLOAT) )
                    qT=bT[i];
                if( (qT != NA_FLOAT) && (qT >= T[i]-EPSILON) )
                    count1[i]+=1;

            } else if( func_cmp == cmp_low ) {

                if( (bT[i] != NA_FLOAT) && (qT != NA_FLOAT) && (bT[i] < qT) )
                    qT=bT[i];
                if( (bT[i] != NA_FLOAT) && (qT == NA_FLOAT) )
                    qT=bT[i];
                if( (qT != NA_FLOAT) && (qT <= T[i]+EPSILON) )
                    count1[i]+=1;

            } else if( func_cmp == cmp_abs ) {

                if( (bT[i] != NA_FLOAT) && (qT != NA_FLOAT) && (fabs(bT[i]) > qT) )
                    qT=fabs(bT[i]);
                if( (bT[i] != NA_FLOAT) && (qT == NA_FLOAT) )
                    qT=fabs(bT[i]);
                if( (qT != NA_FLOAT) && (qT >= fabs(T[i])-EPSILON) )
                    count1[i]+=1;
            }
      
            if( qT != NA_FLOAT )
                total1[i]++;
        }
        b++;

        is_next=(*func_next_sample)(bL);
    }


}      

void set_seed_sampling(long int seed) {
    g_random_seed=seed;
}

