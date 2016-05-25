/**************************************************************************
 *                                                                        *
 *  SPRINT: Simple Parallel R INTerface                                   *
 *  Copyright © 2008,2009 The University of Edinburgh                     *
 *                                                                        *
 *  This program is free software: you can redistribute it and/or modify  *
 *  it under the terms of the GNU General Public License as published by  *
 *  the Free Software Foundation, either version 3 of the License, or     *
 *  any later version.                                                    *
 *                                                                        *
 *  This program is distributed in the hope that it will be useful,       *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the          *
 *  GNU General Public License for more details.                          *
 *                                                                        *
 *  You should have received a copy of the GNU General Public License     *
 *  along with this program. If not, see <http://www.gnu.or/licenses/>.   *
 *                                                                        *
 **************************************************************************/

#include <Rdefines.h>
#include "../../../sprint.h"
#include "../../common/serialize.h"
#include "../../common/utils.h"

SEXP runBootstrapCall(SEXP list, SEXP fn, int worldRank, int worldSize);
SEXP combine_boot_results(SEXP a, SEXP b);

int boot(int n, ...)
{

  int worldSize, worldRank;
  int nbytes;
  SEXP list, fn, tmp;
  SEXP *ret;
  va_list ap;

  int objects_counter = 0;

  MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
  MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);

  if (worldRank == 0) {

      va_start(ap,n);
      list = va_arg(ap, SEXP);
      fn = va_arg(ap, SEXP);
      ret = va_arg(ap, SEXP *);
      va_end(ap);

  }

  /* Broadcast function and data using serialize_form and unserialize_form
    functions from coomon/serialize.c */
  
  if(worldRank == 0) {

    PROTECT(tmp = serialize_form(list)); objects_counter++;
    nbytes = length(tmp);
    
    MPI_Bcast(&nbytes, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(RAW(tmp), nbytes, MPI_BYTE, 0, MPI_COMM_WORLD);

    PROTECT(tmp = serialize_form(fn)); objects_counter++;
    nbytes = length(tmp);
    MPI_Bcast(&nbytes, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(RAW(tmp), nbytes, MPI_BYTE, 0, MPI_COMM_WORLD);
    
  } else {

    MPI_Bcast(&nbytes, 1, MPI_INT, 0, MPI_COMM_WORLD);
    PROTECT(tmp = allocVector(RAWSXP, nbytes)); objects_counter++;
    MPI_Bcast(RAW(tmp), nbytes, MPI_BYTE, 0, MPI_COMM_WORLD);
    PROTECT(list = unserialize_form(tmp)); objects_counter++;

    MPI_Bcast(&nbytes, 1, MPI_INT, 0, MPI_COMM_WORLD);
    tmp = allocVector(RAWSXP, nbytes);
    MPI_Bcast(RAW(tmp), nbytes, MPI_BYTE, 0, MPI_COMM_WORLD);
    PROTECT(fn = unserialize_form(tmp)); objects_counter++;

  }

  PROTECT(tmp = runBootstrapCall(list, fn, worldRank, worldSize));
  objects_counter++;

  MPI_Status status;
  SEXP chunk;

    /* Do a tree reduction across the processes */

  if(worldRank == 0) {

    reduce_combine(tmp, ret, &combine_boot_results, 0, MPI_COMM_WORLD);

  } else {

    reduce_combine(tmp, NULL, &combine_boot_results, 0, MPI_COMM_WORLD);

  }

  UNPROTECT(objects_counter);
  return 0;
 
}

SEXP runBootstrapCall(SEXP list, SEXP fn, int worldRank, int worldSize) {

  SEXP rho = R_GlobalEnv;
  R_len_t i, k, n;
  SEXP R_fcall, ans;
  int my_start, my_end;

  if(!isNewList(list)) error("’list’ must be a list");
  if(!isFunction(fn)) error("’fn’ must be a function");
  if(!isEnvironment(rho)) error("’rho’ should be an environment");
  
  loopDistribute(worldRank, worldSize, length(list), &my_start, &my_end);
  n = my_end - my_start;

  PROTECT(R_fcall = lang2(fn, R_NilValue));
  PROTECT(ans = allocVector(VECSXP, n));
  

  for(i = my_start, k=0; i < my_end; i++, k++) {
    SETCADR(R_fcall, VECTOR_ELT(list, i));
    SET_VECTOR_ELT(ans, k, eval(R_fcall, rho)); 
  }

  setAttrib(ans, R_NamesSymbol, getAttrib(list, R_NamesSymbol));

  
  UNPROTECT(2);

  return(ans);
}

/* Use the interpreter "c" function to combine results */
SEXP combine_boot_results(SEXP a, SEXP b)
{
    SEXP thunk;
    SEXP ret;
    PROTECT(thunk = allocVector(LANGSXP, 3));
    SETCAR(thunk, install("c"));
    SETCADR(thunk, a);
    SETCADDR(thunk, b);
    ret = eval(thunk, R_GlobalEnv);
    UNPROTECT(1);
    return ret;
}
