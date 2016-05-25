/**************************************************************************
 *                                                                        *
 *  SPRINT: Simple Parallel R INTerface                                   *
 *  Copyright © 2008-2011 The University of Edinburgh                     *
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
 *  along with this program. If not, see <http://www.gnu.org/licenses/>.  *
 *                                                                        *
 **************************************************************************/

#include <stdarg.h>
#include "../../../sprint.h"
/* Urg, but I can't think of a better way. */
static char *old_rng_kind[2];

int init_rng_worker(int n, ...)
{
    MPI_Comm comm;
    int rank;
    int size;
    SEXP thunk;
    SEXP package;
    SEXP ret;
    SEXP seed;
    va_list ap;

    comm = MPI_COMM_WORLD;

    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    PROTECT(seed = R_NilValue);
    if ( rank == 0 ) {
        va_start(ap, n);
        seed = va_arg(ap, SEXP);
        va_end(ap);
    } else {
        seed = allocVector(REALSXP, 1);
    }

    MPI_Bcast(REAL(seed), 1, MPI_DOUBLE, 0, comm);
    PROTECT(thunk = allocVector(LANGSXP, 4));

    PROTECT(package = eval(lang2(install("getNamespace"),
                                 ScalarString(mkChar("sprint"))),
                           R_GlobalEnv));
    SETCAR(thunk, findFun(install("init.rng.worker"), package));
    SETCADR(thunk, ScalarInteger(rank));
    SETCADDR(thunk, ScalarInteger(size));
    SETCADDDR(thunk, seed);

    PROTECT(ret = eval(thunk, R_GlobalEnv));

    old_rng_kind[0] = strdup(CHAR(STRING_ELT(ret, 0)));
    old_rng_kind[1] = strdup(CHAR(STRING_ELT(ret, 1)));

    UNPROTECT(4);
    return 0;
}

int reset_rng_worker(int n, ...)
{
    SEXP thunk;
    SEXP arg;
    SEXP package;
    int size;
    MPI_Comm comm;

    comm = MPI_COMM_WORLD;
    MPI_Comm_size(comm, &size);

    PROTECT(thunk = allocVector(LANGSXP, 3));
    PROTECT(package = eval(lang2(install("getNamespace"),
                                 ScalarString(mkChar("sprint"))),
                           R_GlobalEnv));

    SETCAR(thunk, findFun(install("reset.rng.worker"), package));

    PROTECT(arg = allocVector(STRSXP, 2));

    SET_STRING_ELT(arg, 0, mkChar(old_rng_kind[0]));
    SET_STRING_ELT(arg, 1, mkChar(old_rng_kind[1]));

    free(old_rng_kind[0]);
    free(old_rng_kind[1]);

    SETCADR(thunk, arg);

    SETCADDR(thunk, ScalarInteger(size));
    eval(thunk, R_GlobalEnv);

    UNPROTECT(3);

    return 0;
}
