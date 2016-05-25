#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <math.h>

SEXP
bootloop(SEXP fbody, SEXP X, SEXP W, SEXP p, SEXP n, SEXP B, SEXP samp)
{
  int B_len= INTEGER(B)[0], p_len=INTEGER(p)[0], num_samples, b, i, j;
  SEXP Xb, Wb, Sb, Tb, muboot;

  num_samples = INTEGER(n)[0];

  PROTECT(Xb=allocVector(REALSXP,num_samples));
  PROTECT(Wb=allocVector(REALSXP,num_samples));
  PROTECT(Sb=allocVector(INTSXP,num_samples));
  PROTECT(Tb=allocVector(REALSXP,3));
  PROTECT(muboot = allocVector(REALSXP, B_len * p_len));

  SEXP e, ptr;
  PROTECT(e=allocVector(LANGSXP,4)); // this includes the samp
  SETCAR(e, fbody);

  for(b=0;b<B_len;b++){

    if((b%100==0.0) & (b>0.0)) /* modulo 100 */
      Rprintf("%d ",b);

    for(j=0;j<p_len;j++){

      for(i=0;i<num_samples;i++){
	INTEGER(Sb)[i]=INTEGER(samp)[num_samples*b+i];
	REAL(Xb)[i]=REAL(X)[(INTEGER(samp)[num_samples*b+i]-1)*p_len+j];
	REAL(Wb)[i]=REAL(W)[(INTEGER(samp)[num_samples*b+i]-1)*p_len+j];
      }

      ptr = CDR(e);
      SETCAR(ptr, Xb);
      ptr = CDR(ptr);
      SETCAR(ptr, Wb);
      ptr = CDR(ptr);
      SETCAR(ptr, Sb);

      Tb = eval(e, R_GlobalEnv);
      REAL(muboot)[p_len*b+j] = REAL(Tb)[2]*REAL(Tb)[0]/REAL(Tb)[1];
    }
  }

  Rprintf("%d\n",B_len);

  UNPROTECT(6);

  return(muboot);
}

