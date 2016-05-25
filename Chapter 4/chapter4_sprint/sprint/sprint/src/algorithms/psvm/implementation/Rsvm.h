void foote(int a);

struct svm_node ** sparsify (double *x, int r, int c);
struct svm_node ** transsparse (double *x, int r, 
				int *rowindex, 
				int *colindex
				);

void do_cross_validation(struct svm_problem *prob,
			 struct svm_parameter *param,
			 int nr_fold,
			 double* cresults,
			 double* ctotal1,
			 double* ctotal2 );

void do_cross_validation_in_parallel(struct svm_problem *prob,
                                     struct svm_parameter *param,
                                     int nr_fold, int rank,
                                     double* cresults,
                                     double* ctotal1,
                                     double* ctotal2);

void svmtrain (double *x, int *r, int *c, 
	       double *y,
	       int    *rowindex, int *colindex,
	       int    *svm_type,
	       int    *kernel_type,
	       int    *degree,
	       double *gamma,
	       double *coef0,
	       double *cost,
	       double *nu,
	       int    *weightlabels,
	       double *weights,
	       int    *nweights,
	       double *cache,
	       double *tolerance,
	       double *epsilon,
	       int    *shrinking,
	       int    *cross,
	       int    *sparse,
	       int    *probability,
	       int    *seed,
	       int     rank,
	       int     proc,
	       
	       int    *nclasses,
	       int    *nr,
	       int    *index,
	       int    *labels,
	       int    *nSV,
	       double *rho,
	       double *coefs,
	       double *sigma,
	       double *probA,
	       double *probB,

	       double *cresults,
	       double *ctotal1,
	       double *ctotal2,
	       char   **error
	       );

void svmpredict (int    *decisionvalues,
		 int    *probability,
		 
		 double *v, int *r, int *c,
		 int    *rowindex,
		 int    *colindex,
		 double *coefs,
		 double *rho,
		 int    *compprob,
		 double *probA,
		 double *probB,
		 int    *nclasses,
		 int    *totnSV,
		 int    *labels,
		 int    *nSV,
		 int    *sparsemodel,
		 
		 int    *svm_type,
		 int    *kernel_type,
		 int    *degree,
		 double *gamma,
		 double *coef0,
		 
		 double *x, int *xr,
		 int    *xrowindex,
		 int    *xcolindex,
		 int    *sparsex,
		 
		 double *ret,
		 double *dec,
		 double *prob
		 );

void svmwrite (double *v, int *r, int *c,
	       int    *rowindex,
	       int    *colindex,
	       double *coefs,
	       double *rho,
	       double *probA,
	       double *probB,
	       int    *nclasses,
	       int    *totnSV,
	       int    *labels,
	       int    *nSV,
	       int    *sparsemodel,
	       
	       int    *svm_type,
	       int    *kernel_type,
	       int    *degree,
	       double *gamma,
	       double *coef0,
	       
	       char **filename
	       );
