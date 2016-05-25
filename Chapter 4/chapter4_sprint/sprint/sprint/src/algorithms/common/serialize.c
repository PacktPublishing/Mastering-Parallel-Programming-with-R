/**************************************************************************
 *                                                                        *
 *  SPRINT: Simple Parallel R INTerface                                   *
 *  Copyright © 2012 The University of Edinburgh                          *
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

#include "../../sprint.h"

/* When the serialize.c interface is made available to third party
 * packages this could be pulled entirely into the C layer.  As it is
 * we have to go through the interpreter. */
SEXP serialize_form(SEXP form)
{
    SEXP thunk;
    SEXP ret;
    PROTECT(thunk = allocVector(LANGSXP, 3));
    SETCAR(thunk, install("serialize"));
    SETCADR(thunk, form);
    SETCADDR(thunk, R_NilValue);

    ret = eval(thunk, R_GlobalEnv);
    UNPROTECT(1);
    return ret;
}


/* As above. */
SEXP unserialize_form(SEXP form)
{
    SEXP thunk;
    SEXP ret;
    PROTECT(thunk = allocVector(LANGSXP, 2));
    SETCAR(thunk, install("unserialize"));
    SETCADR(thunk, form);

    ret = eval(thunk, R_GlobalEnv);
    UNPROTECT(1);
    return ret;
}

SEXP getListElement(SEXP list, char *str)
{
  SEXP elmt = R_NilValue;
  SEXP names = getAttrib(list, R_NamesSymbol);
  int i;

  for (i = 0; i < length(list); i++)
    if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
      elmt = VECTOR_ELT(list, i);
      break;
    }
  return elmt;
}

void setListElement(SEXP list, char *str, SEXP value)
{
  SEXP names = getAttrib(list, R_NamesSymbol);
  int i;

  for (i = 0; i < length(list); i++)
    if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
      SET_VECTOR_ELT(list, i, value);
      break;
    }
}

/* Do a tree reduction across the processes in COMM combining IN
 * using COMBINE_FN (an associate, commutative, side-effect-free
 * function), writing the result to OUT on ROOT process.  OUT need
 * only be defined on the root process (you can pass NULL on other
 * processes). */
void reduce_combine(SEXP in, SEXP *out, SEXP (*combine_fn)(SEXP, SEXP),
                    int root, MPI_Comm comm)
{
    MPI_Group g;
    MPI_Group doing_comms;
    MPI_Group sitting_out;
    MPI_Comm new;
    MPI_Status status;
    MPI_Request request;
    SEXP ret;
    SEXP tmp;
    int size;
    int rank;
    int neighbour;
    int nbytes;
    int *ranks;
    int i;

    if ( comm == MPI_COMM_NULL ) {
        return;
    }

    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);
    if ( root >= size ) {
        Rprintf("Invalid root rank specified in reduce, aborting\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    /* Rotate ranks so that root is at rank zero, this makes the logic
     * simpler. */
    MPI_Comm_split(comm, 0, (rank - root + size) % size, &new);
    MPI_Comm_dup(new, &comm);
    MPI_Comm_free(&new);
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_group(comm, &g);

    PROTECT(ret = in);
    /* We've reached the root of the tree, assign to the output
     * buffer, which should be defined since we're on the root
     * process, and return. */
    if ( size == 1 ) {
        *out = ret;
        UNPROTECT(1);
        return;
    }

    /* Create group containing the ranks that will participate in this
     * round of communication.  If size is even, this is all of them,
     * otherwise it misses the last one off. */
    ranks = calloc((size_t)size, sizeof(int));
    for ( i = 0; i < size - (size % 2); i++ ) {
        ranks[i] = i;
    }
    MPI_Group_incl(g, size - (size % 2), ranks, &doing_comms);
    /* Create group containing the ranks not participating in the
     * current communication round. */
    MPI_Group_difference(g, doing_comms, &sitting_out);

    MPI_Comm_create(comm, doing_comms, &new);
    /* If we're in the comms round */
    if ( new != MPI_COMM_NULL ) {
        MPI_Comm_rank(new, &rank);
        /* Even rank receives, odd rank sends. */
        if ( rank % 2 == 0 ) {
            neighbour = rank + 1;
            MPI_Recv(&nbytes, 1, MPI_INT, neighbour, 0, new, &status);
            PROTECT(tmp = allocVector(RAWSXP, nbytes));
            MPI_Recv(RAW(tmp), nbytes, MPI_BYTE, neighbour, 0, new, &status);
            ret = (*combine_fn)(in, unserialize_form(tmp));
            UNPROTECT(1);
        } else {
            neighbour = rank - 1; 
            PROTECT(tmp = serialize_form(in));
            nbytes = length(tmp);
            MPI_Isend(&nbytes, 1, MPI_INT, neighbour, 0, new, &request);
            MPI_Isend(RAW(tmp), nbytes, MPI_BYTE, neighbour, 0, new, &request);
            MPI_Wait(&request, &status);
            UNPROTECT(1);
        }
    }

    /* Now create a group containing all the ranks that will
     * participate in the next comms round.  All the even ones from
     * this round... */
    for ( i = 0; i < size/2; i++ ) {
        ranks[i] = i * 2;
    }

    MPI_Group_free(&doing_comms);
    MPI_Group_incl(g, size/2, ranks, &doing_comms);
    MPI_Group_free(&g);
    free(ranks);

    /* ...plus the one that was left out in this round. */
    MPI_Group_union(doing_comms, sitting_out, &g);
    MPI_Group_free(&doing_comms);
    MPI_Group_free(&sitting_out);

    if ( new != MPI_COMM_NULL ) {
        MPI_Comm_free(&new);
    }
    /* Create communicator */
    MPI_Comm_create(comm, g, &new);
    MPI_Group_free(&g);
    MPI_Comm_free(&comm);
    /* Do next level of reduction.  Now the root is zero, because we
     * rotated ranks in the top-level call. */
    reduce_combine(ret, out, combine_fn, 0, new);

    if ( new != MPI_COMM_NULL ) {
        MPI_Comm_free(&new);
    }

    UNPROTECT(1);
    return;
}
