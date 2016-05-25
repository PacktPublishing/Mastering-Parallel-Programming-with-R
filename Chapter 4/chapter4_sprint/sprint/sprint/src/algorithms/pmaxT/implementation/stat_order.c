
#include "stdio.h"
#include "stdlib.h"
#include "stdarg.h"
#include "math.h"
#include "mt.h"
#include "../../../sprint.h"

static double* gp_arr;   

typedef struct tagCMP_DATA{
    double *V;
    FUNC_CMP func_cmp;
}CMP_DATA;

// The lower are used for the order_mult_data
CMP_DATA *gp_cmp_data;
int g_ncmp;
static int cmp_mult(const void* v1, const void* v2);


/*

 */
static int cmp_mult(const void* v1, const void* v2)
{
  int i;
  int ret=-2;
  for(i=0;i<g_ncmp;i++){
    gp_arr=gp_cmp_data[i].V; /*used in the function func_cmp provided in this file*/
    ret=(gp_cmp_data[i].func_cmp)(v1,v2);
    if(ret!=0) return ret;/*return the result*/		      
  }
  return ret;
}


/******************************************************************************
 *                               order_data                                   *
 *                               ----------                                   *
 * n: the dimension of the data                                               *
 * V1, func_cmp1,...Vk,funct_cmpk                                             *
 * are the parameters to compare such that using the comparing function       *
 * func_cmpi to                                                               *
 * order Vi, every Vi needs to be double array                                *
 * the result is stored at int array R.                                       *
 ********************************************************************************/
void order_mult_data(int* R,int n,int k,...)
{
    CMP_DATA *cmp_data;
    va_list ap;
    int i;

    // Allocate, initialize and check memory allocation
    cmp_data = (CMP_DATA*)R_alloc(k, sizeof(CMP_DATA));

    // Get input arguments
    va_start(ap,k);
    for(i=0; i<k; i++) {
        cmp_data[i].V = va_arg(ap, double*);
        cmp_data[i].func_cmp = va_arg(ap, FUNC_CMP);
    }
    va_end(ap);

    // Both used in the function cmp_mult
    gp_cmp_data = cmp_data;
    g_ncmp = k; 
    for(i=0; i<n; i++)
        R[i] = i;
    qsort(R, n, sizeof(R[0]), cmp_mult);

}
    
void order_data(double* V,int*R,int n,FUNC_CMP func_cmp)
{
    int i;
    for(i=0;i<n;i++)
        R[i]=i;
    gp_arr=V;
    qsort(R,n,sizeof(R[0]),func_cmp);
}
  

/********************************************************************************
 *                          Comparing functions                                 *
 *                          -------------------                                 *
 *  1 cmp_abs: comparing the absolute values,                                   *
 *  2 cmp_low: comparing the lower tail                                         *
 *  3 cmp_high: comparing the higher tail                                       *
 *                                                                              *
 *  Note                                                                        *
 *  ----                                                                        *
 *  The gp_arr is the global pointer which has to be used in the program qsort  *
 *  After using the quick with those functions, it will be odered such that     *
 *  1 in cmp_abs: the bigger values in absolute values will in lower index      *
 *                similar to two sided-test                                     *
 *  2 in cmp_low: the smaller values will in lower index                        *
 *                similar to lower tail test                                    *
 *  3 in cmp_high:the bigger values will in lower index                         *
 *                similar to hight tail test.                                   *
 *                                                                              *
 *always put the absolute value at the end of the array                         *
 ********************************************************************************/
int cmp_abs(const void *v1, const void *v2) {

    double f1=fabs(*(gp_arr+*(int *)v1));
    double f2=fabs(*(gp_arr+*(int *)v2));

    if(f1==NA_FLOAT)
        return 1;
    if(f2==NA_FLOAT)
        return -1;
    if (f1<f2)
        return 1;
    if (f1>f2)
        return -1;
    else
        return 0; 
}

int cmp_low(const void *v1, const void *v2) {

    if((*(gp_arr+*(int *)v1))==NA_FLOAT)
        return 1;
    if((*(gp_arr+*(int *)v2))==NA_FLOAT)
        return -1;
    if ((*(gp_arr+*(int *)v1))<(*(gp_arr+*(int *)v2)))
        return -1;
    if ((*(gp_arr+*(int *)v1))>(*(gp_arr+*(int *)v2)))
        return 1;
    else
        return 0; 
} 

int cmp_high(const void *v1, const void *v2) {

    if((*(gp_arr+*(int *)v1))==NA_FLOAT)
        return -1;
    if((*(gp_arr+*(int *)v2))==NA_FLOAT)
        return 1;
    if ((*(gp_arr+*(int *)v1))<(*(gp_arr+*(int *)v2)))
        return 1;
    if ((*(gp_arr+*(int *)v1))>(*(gp_arr+*(int *)v2)))
        return -1;
    else
        return 0; 
}

