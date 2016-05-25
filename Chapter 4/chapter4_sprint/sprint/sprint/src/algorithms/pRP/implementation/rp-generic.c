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

#define _GNU_SOURCE
#include <stdlib.h>
#include <stdbool.h>
#include <mpi.h>
#include "rank-product.h"
#include "../../../sprint.h"

/* Comparison function sorting a before b */
inline int compare_fcs_up(const void *a, const void *b)
{
    const double *da = &((const fold_change_t *) a)->val;
    const double *db = &((const fold_change_t *) b)->val;

    return (*da > *db) - (*da < *db);
}

/* Comparison function sorting b before a */
inline int compare_fcs_down(const void *a, const void *b)
{
    return -compare_fcs_up(a, b);
}

inline int compare_doubles(const void *a, const void *b)
{
    const double *da = (const double *) a;
    const double *db = (const double *) b;

    return (*da > *db) - (*da < *db);
}

/* Return the position in BUF of VAL.  If VAL doesn't match an entry
 * exactly, return the position of the largest entry in BUF smaller
 * than VAL. */
int bsearch_approx(double *buf, const double val, const size_t ngenes)
{
    double tmp;
    int high;
    int low;
    int mid;

    high = ngenes;
    low = 0;

    while ( low < high ) {
        mid = low + (high - low)/2;
        tmp = buf[mid];
        if ( val > tmp )
            low = mid + 1;
        else
            high = mid;
    }
    return low;
}

/* Calculate the fold change (pointwise division) v1/v2 and put the
 * values into BUF, which must be initialised. */
void fold_change(fold_change_t *buf,
                 const double *v1, const double *v2,
                 const size_t n)
{
    size_t i;
    
    if ( NULL == v2 ) {
        for ( i = 0; i < n; i++ ) {
            buf[i].idx = i;
            buf[i].val = v1[i];
        }
    } else {
        for ( i = 0; i < n; i++ ) {
            buf[i].idx = i;
            buf[i].val = v1[i] / v2[i];
        }
    }
}

/* Calculate the log fold change (pointwise subtraction) v1 - v2 and
 * put the values into BUF, which must be initialised. */
void log_fold_change(fold_change_t *buf,
                     const double *v1, const double *v2,
                     const size_t n)
{
    size_t i;

    if ( NULL == v2 ) {
        for ( i = 0; i < n; i++ ) {
            buf[i].idx = i;
            buf[i].val = v1[i];
        }
    } else {
        for ( i = 0; i < n; i++ ) {
            buf[i].idx = i;
            buf[i].val = v1[i] - v2[i];
        }
    }
}

/* Compute fold changes of v1 against v2 and rank them in increasing
 * (or decreasing if rev_sorting is true) order.  BUF must be
 * initialised beforehand and afterwards the entries of BUF contain
 * the fold-change value, the original gene this change came from,
 * and the rank of this fold-change. */
void rank_fold_changes(fold_change_t *buf,
                       const double *v1,
                       const double *v2,
                       const size_t n,
                       const bool logarithmic_data,
                       const bool rev_sorting)
{
    size_t i;
    size_t j;
    double tmp;
    double sum;
    size_t last;

    /* Fill buf with fold-change data */
    if ( logarithmic_data ) {
        log_fold_change(buf, v1, v2, n);
    } else {
        fold_change(buf, v1, v2, n);
    }
    
    qsort(buf, n, sizeof(fold_change_t),
          rev_sorting ? compare_fcs_down : compare_fcs_up);

    /* R's rank function breaks ties by averaging.  So do that here
     * too.  We walk over the sorted list of fold changes.  While the
     * fold change for neighbouring elements is the same, we continue
     * to increment the sum of the ranks.  When a different value is
     * encountered we go back and write in the average rank of all
     * the fold-changes we just skipped over. */
    for ( i = 0; i < n; ) {
        tmp = buf[i].val;
        last = i;
        sum = 0;
        do {                    /* Keep going until a different value
                                 * is encountered */
            ++i;
            sum += i;
        } while ( i < n         /* Don't read off the end of buf */
                  && tmp == buf[i].val )
            ;
        /* Go back and write the "rank" in to the slots we skipped
         * over. */
        for ( j = last; j < i; j++ ) {
            buf[j].rank = sum / (i - last);
        }
    }
}

/* In place shuffle (Fisher-Yates) of a single gene vector */
void shuffle(double *data, const size_t n)
{
    size_t i;
    size_t j;
    double tmp;
    for ( i = n - 1; i > 0; i-- ) {
        j = (size_t)(unif_rand() * i);
        tmp = data[j];
        data[j] = data[i];
        data[i] = tmp;
    }
}

/* Shuffle entire dataset (nclass gene vectors) */
void shuffle_data(double *data, const int nclass, const size_t ngenes)
{
    int i;
    GetRNGstate();
    for ( i = 0; i < nclass; i++ ) {
        shuffle(data + i*ngenes, ngenes);
    }
    PutRNGstate();
}

void shuffle_data_multi(double *data, const int *nclass, const int norigins,
                        const size_t ngenes)
{
    int i;
    size_t offset;
    GetRNGstate();
    offset = 0;
    for ( i = 0; i < norigins; i++ ) {
        shuffle_data(data + offset, nclass[i], ngenes);
        offset += ngenes * nclass[i];
    }
    PutRNGstate();
}

void *xmalloc(const size_t n, const MPI_Comm comm)
{
    register void *ret;

    ret = malloc(n);

    if ( NULL == ret ) {
        Rprintf("Unable to allocate memory\n");
        MPI_Abort(comm, EXIT_FAILURE);
    }

    return ret;
}

void *xcalloc(const size_t n, const size_t eltsize, const MPI_Comm comm)
{
    register void *ret;

    ret = calloc(n, eltsize);

    if ( NULL == ret ) {
        Rprintf("Unable to allocate memory\n");
        MPI_Abort(comm, EXIT_FAILURE);
    }

    return ret;
}
