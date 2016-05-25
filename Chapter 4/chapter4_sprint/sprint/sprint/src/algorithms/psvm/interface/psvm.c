#include <Rdefines.h>
#include "../../../functions.h"
#include "../../../sprint.h"
#include "psvm.h"

extern int svm_call(int n, ...);



SEXP psvm(SEXP x, 
	  SEXP rows, 
	  SEXP cols, 
	  SEXP y, 
	  SEXP nclass, 
	  SEXP cross, 
	  SEXP sparseRows, 
	  SEXP sparseCols,
	  SEXP type, 
	  SEXP kernel, 
	  SEXP degree, 
	  SEXP gamma, 
	  SEXP coef0, 
	  SEXP cost,
	  SEXP nu, 
	  SEXP weightlabels, 
	  SEXP classweight, 
	  SEXP lcw, 
	  SEXP cache, 
	  SEXP tolerance,
	  SEXP epsilon, 
	  SEXP shrinking, 
	  SEXP sparse,
	  SEXP probability, 
	  SEXP seed )
{
  
  SEXP result;
  int response,intCode;
  params* p = (params*)malloc(sizeof(params));
  
  // store x matrix and label vector
  double* x_matrix;
  double* y_vector;
  
 
  enum commandCodes commandCode;
  
  int message = 10;

  PROTECT( x      = AS_NUMERIC( x      ) );
  PROTECT( y      = AS_NUMERIC( y      ) );
  
  PROTECT( cross  = AS_INTEGER( cross  ) );
  PROTECT( rows   = AS_INTEGER( rows   ) );
  PROTECT( cols   = AS_INTEGER( cols   ) );
  PROTECT( nclass = AS_INTEGER( nclass ) );
  
     
  PROTECT( sparseRows  = AS_INTEGER( sparseRows  ) );
  PROTECT( sparseCols  = AS_INTEGER( sparseCols ) );
  
  PROTECT( type    = AS_INTEGER( type  ) );
  PROTECT( kernel  = AS_INTEGER( kernel  ) );
  PROTECT( degree  = AS_INTEGER( degree  ) );
  PROTECT( gamma   = AS_NUMERIC( gamma  ) );
  
  PROTECT( coef0  = AS_NUMERIC( coef0  ) );
  PROTECT( cost  = AS_NUMERIC( cost  ) );
  PROTECT( nu  = AS_NUMERIC( nu  ) );
  
  PROTECT( weightlabels = AS_INTEGER(weightlabels  ) );
  PROTECT( classweight  = AS_NUMERIC(classweight ) );
  PROTECT( lcw  = AS_INTEGER(lcw  ) );

  PROTECT( cache  = AS_NUMERIC(cache  ) );
  PROTECT( tolerance  = AS_NUMERIC(tolerance  ) );
  PROTECT( epsilon  = AS_NUMERIC(epsilon  ) );
 

  PROTECT( shrinking  = AS_INTEGER( shrinking  ) );
  PROTECT( sparse  = AS_INTEGER( sparse  ) );
  PROTECT( probability  = AS_INTEGER( probability  ) );
  PROTECT( seed  = AS_INTEGER( seed  ) );
  
 
  x_matrix = NUMERIC_POINTER(x);
  y_vector = NUMERIC_POINTER(y);
  
    
  p->rows                = INTEGER_VALUE(rows);
  p->columns             = INTEGER_VALUE(cols);
  p->nclasses            = INTEGER_VALUE(nclass);
  p->cross               = INTEGER_VALUE(cross);
  p->type                = INTEGER_VALUE(type);
  p->sparseRows          = INTEGER_VALUE(sparseRows);
  p->sparseCols          = INTEGER_VALUE(sparseCols);
  p->degree              = INTEGER_VALUE(degree);
  p->kernel              = INTEGER_VALUE(kernel);
  p->whiteLabels         = INTEGER_VALUE(weightlabels);
  p->lengthWhiteLabels   = INTEGER_VALUE(lcw);
  p->shrinking           = INTEGER_VALUE(shrinking);
  p->sparse              = INTEGER_VALUE(sparse);
  p->probability         = INTEGER_VALUE(probability);
  p->seed                = INTEGER_VALUE(seed);
  
  
  p->gamma           = NUMERIC_VALUE(gamma);
  p->coef0           = NUMERIC_VALUE(coef0);
  p->cost            = NUMERIC_VALUE(cost);
  p->nu              = NUMERIC_VALUE(nu);
  p->classWeights    = NUMERIC_VALUE(classweight);
  p->cacheSize       = NUMERIC_VALUE(cache);
  p->tolerance       = NUMERIC_VALUE(tolerance);
  p->epsilon         = NUMERIC_VALUE(epsilon);
  
  

  //We don't need to keep data any longer
  //so we can unprotect it
  
  UNPROTECT(25);
  
      
  MPI_Initialized(&response);
  if (response) {
    DEBUG("MPI is init'ed in ptest\n");
  } else {
    DEBUG("MPI is NOT init'ed in ptest\n");
  }
  
 
 
  // broadcast command to other processors
  commandCode = PSVM;
  intCode = (int)commandCode;
  
  MPI_Bcast(&intCode, 1, MPI_INT, 0, MPI_COMM_WORLD);
  
  response = svm_call( 3, x_matrix, y_vector, p);
  
  //result = ScalarInteger(response);
  
  result = ScalarInteger(2);
  
  return result;
}
