/**************************************************************************
 *                                                                        *
 *  SPRINT: Simple Parallel R INTerface                                   *
 *  Copyright © 2008,2009,2010 The University of Edinburgh                *
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

#include <stdarg.h>
#include "../../../sprint.h"
#include "../../common/serialize.h"
#include <R_ext/Parse.h>
#include <R.h>

/* These match up with the R interface.  So that nthcdr(args,
 * arguments_t) gives the relevant SEXP value. */
typedef enum _arguments_t {
    X = 0,
    Y,
    XTEST,
    YTEST,
    NTREE,
    MTRY,
    REPLACE,
    CLASSWT,
    CUTOFF,
    STRATA,
    SAMPSIZE,
    NODESIZE,
    MAXNODES,
    IMPORTANCE,
    LOCALIMP,
    NPERM,
    PROXIMITY,
    OOB_PROX,
    NORM_VOTES,
    DO_TRACE,
    KEEP_FOREST,
    COR_BIAS,
    KEEP_INBAG
} arguments_t;

/* Call the serial randomForest call with the args we were passed in
 * parallel.  We use the size and rank arguments to calculate the
 * subset of the total trees we're going to compute. */
SEXP serial_randomForest(SEXP args, int rank, int size)
{
    SEXP thunk;
    SEXP ret;
    SEXP tmp;
    int length;
    int i;
    int ntree;
    int chunksize;
    int remainder;
#ifdef _PROF
    double t;
    t = MPI_Wtime();
#endif
    length = length(args) + 1;  /* Space for function name too */

    PROTECT(thunk = allocVector(LANGSXP, length));
    SETCAR(thunk, install("randomForest"));

    ntree = INTEGER(coerceVector(CAR(nthcdr(args, NTREE)), INTSXP))[0];

    if ( ntree < size ) {
        UNPROTECT(1);
        return NULL;
    }
    chunksize = ntree/size;
    remainder = ntree - (chunksize * size);
    if ( rank < remainder ) {
        ++chunksize;
    }
    /* We modify this value in place but it's come from R, so we'll need to
     * restore the value afterwards */
    REAL(CAR(nthcdr(args, NTREE)))[0] = (double)chunksize;

    /* Need to keep pointer to head of args list, so just modify a
     * pointer copy here. */
    PROTECT(tmp = args);
    for ( i = 1; i < length; i++ ) {
        /* Because args is a pairlist, but thunk is a linked list, Nil
         * values in thunk confuse the R->C call that the evaluation
         * below carries out.  Basically the Nil value is seen as the
         * end of the list, rather than just a null argument.  To fix
         * this, if we got any Nil arguments through the R call, we
         * replace them with "missing" values.  These get fixed up in
         * the serial randomForest call. */
        if ( CAR(tmp) == R_NilValue ) {
            SET_MISSING(CAR(tmp), 1);
            SETCAR(tmp, R_MissingArg);
        }
        SETCAR(nthcdr(thunk, i), CAR(tmp));
        tmp = CDR(tmp);
    }
    UNPROTECT(1);
    PROTECT(ret = eval(thunk, R_GlobalEnv));

    /* Don't need to send the data back over, so set the "call" part of
     * the list to NULL, we repopulate in the R interface. */
    setListElement(ret, "call", R_NilValue);
    /* The xlevels and ncat variables are the same for a given
     * dataset, so we don't need to transfer these over either if
     * we're not rank 0. */
    if ( rank != 0 ) {
        PROTECT(tmp = getListElement(ret, "forest"));
        setListElement(tmp, "xlevels", R_NilValue);
        setListElement(tmp, "ncat", R_NilValue);
        UNPROTECT(1);
    }
    /* Need to restore saved number of trees value */
    REAL(CAR(nthcdr(args, NTREE)))[0] = (double)ntree;
    UNPROTECT(2);
#ifdef _PROF
    Rprintf("[%d: subforest-gen] %g (%d trees)\n", rank, MPI_Wtime() - t, chunksize);
#endif
    return ret;
}

/* Load a package.  We need this so that we can load the randomForest
 * package on the slaves. */
static void require_package(char *package_name)
{
    SEXP thunk;
    PROTECT(thunk = allocVector(LANGSXP, 2));
    SETCAR(thunk, install("require"));
    SETCADR(thunk, install(package_name));
    eval(thunk, R_GlobalEnv);
    UNPROTECT(1);
}

/* Combining forests uses the R function pcombine we've just
 * defined which calls out to the random forest library function
 * combine to do most of the heavy lifting. */
SEXP combine_forests(SEXP a, SEXP b)
{
    SEXP thunk;
    SEXP ret;
    PROTECT(thunk = allocVector(LANGSXP, 3));
    SETCAR(thunk, install("pcombine"));
    SETCADR(thunk, a);
    SETCADDR(thunk, b);
    ret = eval(thunk, R_GlobalEnv);
    UNPROTECT(1);
    return ret;
}

/* Workhorse function.  This is where we jump into as a slave. */
int random_forest_driver(int n, ...)
{
    int size;
    int rank;
    MPI_Comm comm;
    va_list ap;
    SEXP args;               /* args to randomForest */
    SEXP tmp;                /* serialized data */
    SEXP *ret;               /* result passed back (on root) */
    int length;
#ifdef _PROF
    double t;
#endif
    comm = MPI_COMM_WORLD;

    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    if ( 0 == rank ) {
        va_start(ap, n);
        args = va_arg(ap, SEXP);
        ret = va_arg(ap, SEXP *);
        va_end(ap);

#ifdef _PROF
        t = MPI_Wtime();
#endif
        /* Throw away name of function call */
        PROTECT(args = CDR(args));
        /* Serialize input call.  The data come through here too. */
        PROTECT(tmp = serialize_form(args));

        length = length(tmp);

        /* Send this serialized chunk to everyone else. */
        MPI_Bcast(&length, 1, MPI_INT, 0, comm);
        MPI_Bcast(RAW(tmp), length, MPI_BYTE, 0, comm);
        UNPROTECT(1);           /* tmp */
#ifdef _PROF
        Rprintf("[0: setup-and-send] %g (%d bytes)\n", MPI_Wtime() - t, length);
#endif
        /* Do some work on master too */
        PROTECT(tmp = serial_randomForest(args, rank, size));
        if ( NULL == tmp ) {
            Rprintf("Can't call prandomForest with more processes than trees, returning NULL\n");
            *ret = R_NilValue;
            UNPROTECT(2);
            return 0;
        }

#ifdef _PROF
        t = MPI_Wtime();
#endif
        /* Get the results back */
        reduce_combine(tmp, ret, &combine_forests, 0, comm);
#ifdef _PROF
        Rprintf("[%d: combine-forests] %g\n", rank, MPI_Wtime() - t);
#endif
        UNPROTECT(2);           /* tmp and args */
    } else {
#ifdef _PROF
        t = MPI_Wtime();
#endif
        /* Make sure randomForest library is loaded on slave */
        require_package("randomForest");
        /* Get data */
        MPI_Bcast(&length, 1, MPI_INT, 0, comm);
        PROTECT(tmp = allocVector(RAWSXP, length));
        MPI_Bcast(RAW(tmp), length, MPI_BYTE, 0, comm);
        PROTECT(args = unserialize_form(tmp));
        UNPROTECT_PTR(tmp);
#ifdef _PROF
        Rprintf("[%d: setup-and-recv] %g (%d bytes)\n", rank, MPI_Wtime() - t, length);
#endif
        /* Do our bit of work */
        PROTECT(tmp = serial_randomForest(args, rank, size));
        if ( NULL == tmp ) {
            UNPROTECT(2);
            return 0;
        }

#ifdef _PROF
        t = MPI_Wtime();
#endif
        reduce_combine(tmp, NULL, &combine_forests, 0, comm);
#ifdef _PROF
        Rprintf("[%d: combine-forests] %g\n", rank, MPI_Wtime() - t);
#endif
        UNPROTECT(2);           /* tmp and args */
    }
    return 0;
}
