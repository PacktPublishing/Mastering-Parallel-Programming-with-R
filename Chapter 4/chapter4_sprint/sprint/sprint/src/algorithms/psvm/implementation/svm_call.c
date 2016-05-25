#include <stdio.h>
#include "svm.h"
#include "Rsvm.h"
#include "mpi.h"
#include "../interface/psvm.h"

#include "../../../functions.h"
#include "../../../sprint.h"

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <Rdefines.h>
#include <stdarg.h>
void print_svmType(int n);
void print_svmKernel(int n);
void printPSVM(params* par, int* nPSV,int* nclasses);

int svm_call(int n, ...)
{
  int i;
  double *x, *y;
  int rows, columns, result;
  int nclass, cross;
  
  va_list ap;
  
  int worldSize, worldRank;

  int    *nclasses;
  int    *nr;
  int    *index;
  int    *labels;
  int    *nSV;
  params *pp;
  double *rho;
  double *coefs;
  double *sigma;
  double *probA;
  double *probB;
  
  double *cresults;
  double *ctotal1;
  double *ctotal2;
  char   **error;
  
  //make room for svm parameters
  pp = (params*)malloc(sizeof(params));

  MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
  MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);
  MPI_Status stat; 
  
  MPI_Datatype MPI_parType, oldtypes[2];
  int          blockcounts[2]; 
  MPI_Aint    offsets[2], extent; 

  // setup description of the 15 MPI_INT fields  
  offsets[0] = 0; 
  oldtypes[0] = MPI_INT; 
  blockcounts[0] = 15; 
  
  // setup description of the 8 MPI_DOUBLE fields n, type  
  // Need to first figure offset by getting size of MPI_INT
  MPI_Type_extent(MPI_INT, &extent); 
  offsets[1] = 15 * extent; 
  oldtypes[1] = MPI_DOUBLE; 
  blockcounts[1] = 9; 

  // Now define structured type and commit it 
  MPI_Type_struct(2, blockcounts, offsets, oldtypes, &MPI_parType); 
  MPI_Type_commit(&MPI_parType); 
  
  // Get input variables
  if (worldRank == 0) {

    // Get input variables
    va_start(ap,n);

    //matrix
    x       = va_arg(ap,double*);
    y       = va_arg(ap,double*);

    //svm parameters
    pp      = va_arg(ap,params*);

    //no more parameters expected
    va_end(ap);
  } 
  
  //broadcast parameters to all nodes
  MPI_Bcast(pp, 1, MPI_parType, 0, MPI_COMM_WORLD);

  //make room for out parameters
  nclasses = (int*)malloc(sizeof(int));
  nr       = (int*)malloc(sizeof(int));
  index    = (int*)malloc(pp->rows * sizeof(int));
  labels   = (int*)malloc(pp->nclasses * sizeof(int));
  nSV      = (int*)malloc(pp->nclasses * sizeof(int));

  rho      = (double*)malloc((pp->nclasses * (pp->nclasses - 1) / 2)* sizeof(double));
  coefs    = (double*)malloc((pp->rows   * (pp->nclasses - 1))    * sizeof(double));
  sigma    = (double*)malloc(sizeof(double));
  probA    = (double*)malloc((pp->nclasses * (pp->nclasses - 1) / 2)* sizeof(double));
  probB    = (double*)malloc((pp->nclasses * (pp->nclasses - 1) / 2)* sizeof(double)); 
  
  cresults = (double*)malloc(pp->cross * sizeof(double));
  ctotal1  = (double*)malloc(sizeof(double));
  ctotal2  = (double*)malloc(sizeof(double));

  error    = (char**)malloc(255*sizeof(char));
  for(i = 0; i< 255; i++)
    error[i] = (char*)malloc(255*sizeof(char));
  

  if(worldRank != 0){
    x = (double*)malloc(pp->rows * pp->columns * sizeof(double));  
    y = (double*)malloc(pp->rows * sizeof(double)); 
  }

  //broadcast input matrix 
  MPI_Bcast(x, pp->rows * pp->columns, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(y, pp->rows, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  
  //if( worldRank!=0 ){
  
  svmtrain (x, &(pp->rows), &(pp->columns),
	    y,
	    &(pp->sparseRows), 
	    &(pp->sparseCols),
	    &(pp->type),
	    &(pp->kernel),
	    &(pp->degree),
	    &(pp->gamma),
	    &(pp->coef0),
	    &(pp->cost),
	    &(pp->nu),
	    &(pp->whiteLabels),
	    &(pp->classWeights),
	    &(pp->lengthWhiteLabels),
	    &(pp->cacheSize),
	    &(pp->tolerance),
	    &(pp->epsilon),
	    &(pp->shrinking),
	    &(pp->cross),
	    &(pp->sparse),
	    &(pp->probability),
	    &(pp->seed),
	    worldRank,
	    worldSize,
	    
	    //Outputs
	    nclasses,
	    nr,
	    index,
	    labels,
	    nSV,
	    rho,
	    coefs,
	    sigma,
	    probA,
	    probB,
	    
	    cresults,
	    ctotal1,
	    ctotal2,
	    error);
  
  if( worldRank == 0 )  
    printPSVM(pp, nSV, nclasses);
  //  }
  result = 0;
  return result;
}


void printPSVM(params* par, int* nPSV, int* nclasses){

  int i, nTSV;
  Rprintf("\n\n");

  print_svmType(par->type);
  print_svmKernel(par->kernel);

  Rprintf("\n\n");
  Rprintf("Parameters :\n");
  
    
  if(par->type == 0 || par->type ==3 || par->type == 4)
    Rprintf("cost: %g\n", par->cost);
  
  if(par->kernel == 1)
    Rprintf("degree: %d\n", par->degree);
  
  Rprintf("gamma: %g \n", par->gamma);
  
  if(par->kernel == 1 || par->kernel == 3)
    Rprintf("coef 0: %g\n", par->coef0);
  
  if(par->type == 1 || par->type == 2 || par->type == 4)
    Rprintf("nu: %g \n", par->nu);
  
  if(par->type == 3)
    Rprintf("epsilon: %g\n", par->epsilon);
  
  nTSV = 0;
  for(i = 0; i < *nclasses; i++ )
    nTSV += nPSV[i];
  
  Rprintf("Number of Support Vectors: %d \n",nTSV);
  Rprintf("\n\n");
  
   
  if(par->type << 2)
    Rprintf("Number of Classes: %d\n", par->nclasses);
  
}


void print_svmType(int n){

  switch(n){
    
  case 0:
    Rprintf(" svm-Type: C-classification\n");
    break;
   
  case 1:
    Rprintf(" svm-Type: nu-classification\n");
    break;
    
  case 2:
    Rprintf(" svm-Type: one-classification\n");
    break;
    
  case 3:
    Rprintf(" svm-Type: eps-regression\n");
    break;
    
  case 4:
    Rprintf(" svm-Type: nu-regression\n");
    break;
    
  default:
    Rprintf(" svm-Type: error \n");
    break;
    
  }
}

void print_svmKernel(int n){

  switch(n){
    
  case 0:
    Rprintf(" svm-Kernel: Linear\n");
    break;
   
  case 1:
    Rprintf(" svm-Kernel: Polynomial\n");
    break;
    
  case 2:
    Rprintf(" svm-Kernel: Radial\n");
    break;
    
  case 3:
    Rprintf(" svm-Kernel: Sigmoid\n");
    break;
    
  default:
    Rprintf(" svm-Kernel: error \n");
    break;
    
  }
}
