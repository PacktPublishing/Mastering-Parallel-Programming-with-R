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
#include <stdarg.h>
#include "rank-product.h"
#include "../../../sprint.h"

/* Bootstrap a null distribution for the Rank Product statistic. */
void boot_rp(double *data1, const int nclass1,
             double *data2, const int nclass2,
             const double *experimental_rp,
             const size_t ngenes,
             const bool logarithmic_data, const bool rev_sorting,
             const bool rank_sum,
             const int nperms,
             double *ret)

{
    int i;
    int j;
    fold_change_t *buf;
    double *boot_rp;
    /* Pre-allocate to avoid memory churn. */
    buf = xmalloc(sizeof(fold_change_t) * ngenes, MPI_COMM_WORLD);
    boot_rp = xmalloc(sizeof(double) * ngenes, MPI_COMM_WORLD);

    for ( i = 0; i < ngenes; i++ ) {
        ret[i] = 0.0;
    }
    /* Generate random experiment and calculate the rank product.
     * Repeat many times to obtain a distribution */
    for ( i = 0; i < nperms; i++ ) {
        shuffle_data(data1, nclass1, ngenes);
        shuffle_data(data2, nclass2, ngenes);
        rank_product(buf, data1, nclass1, data2, nclass2, ngenes,
                     logarithmic_data, rev_sorting, rank_sum, boot_rp);
        /* For each experimental rank product, count the number of
         * generated rank products that are smaller.  I.e., where in
         * the distribution are we? This is O(ngenes log(ngenes))
         * since the search is log(ngenes) for each entry. */
        qsort(boot_rp, ngenes, sizeof(double), compare_doubles);
        for ( j = 0; j < ngenes; j++ ) {
            ret[j] += bsearch_approx(boot_rp, experimental_rp[j], ngenes);
        }
    }
    free(boot_rp);
    free(buf);
}

/* Bootstrap null distribution for data from multiple origins. */
void boot_rp_multi(double *data1, const int *nclass1,
                   double *data2, const int *nclass2,
                   const double *experimental_rp,
                   const size_t ngenes,
                   const int norigins,
                   const bool logarithmic_data,
                   const bool rev_sorting,
                   const bool rank_sum,
                   const int nperms,
                   double *ret)
{
    int i;
    int j;
    fold_change_t *buf;
    double *boot_rp;
    /* Pre-allocate to avoid memory churn. */
    buf = xmalloc(sizeof(fold_change_t) * ngenes, MPI_COMM_WORLD);
    boot_rp = xmalloc(sizeof(double) * ngenes, MPI_COMM_WORLD);

    for ( i = 0; i < ngenes; i++ ) {
        ret[i] = 0.0;
    }

    for ( i = 0; i < nperms; i++ ) {
        shuffle_data_multi(data1, nclass1, norigins, ngenes);
        shuffle_data_multi(data2, nclass2, norigins, ngenes);
        rank_product_multi(buf, data1, nclass1, data2, nclass2, ngenes,
                           norigins, logarithmic_data, rev_sorting,
                           rank_sum, boot_rp);
        /* For each experimental rank product, count the number of
         * generated rank products that are smaller.  I.e., where in
         * the distribution are we? This is O(ngenes log(ngenes))
         * since the search is log(ngenes) for each entry. */
        qsort(boot_rp, ngenes, sizeof(double), compare_doubles);
        for ( j = 0; j < ngenes; j++ ) {
            ret[j] += bsearch_approx(boot_rp, experimental_rp[j], ngenes);
        }
    }
    free(boot_rp);
    free(buf);
}

void gather_boot_data(double *data, const size_t genes,
                      double *ret, MPI_Comm comm)
{
    MPI_Reduce(data, ret, genes, MPI_DOUBLE, MPI_SUM, 0, comm);
}

/* Interface that dispatches to the one- or two-class implementation
 * as appropriate */
int boot_rank_product(int n, ...)
{
    double *data1_ = NULL;
    int nclass1;
    double *data2_ = NULL;
    int nclass2;
    double *experimental_rp = NULL;
    int ngenes;
    int logarithmic_data;
    int rev_sorting;
    int rank_sum;
    int chunk;
    int nperms;
    double *ret = NULL;
    double *local_distribution = NULL;
    double *data1 = NULL;
    double *data2 = NULL;
    int rank;
    int size;
    MPI_Comm comm;
    va_list ap;

    comm = MPI_COMM_WORLD;

    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    if ( 0 == rank ) {
        va_start(ap, n);

        data1_ = va_arg(ap, double *);
        nclass1 = va_arg(ap, int);
        data2_ = va_arg(ap, double *);
        nclass2 = va_arg(ap, int);
        experimental_rp = va_arg(ap, double *);
        ngenes = (int)va_arg(ap, size_t);
        logarithmic_data = va_arg(ap, int);
        rev_sorting = va_arg(ap, int);
        rank_sum = va_arg(ap, int);
        nperms = va_arg(ap, int);
        ret = va_arg(ap, double *);

        va_end(ap);
    }

    MPI_Bcast(&nclass1, 1, MPI_INT, 0, comm);
    MPI_Bcast(&nclass2, 1, MPI_INT, 0, comm);
    MPI_Bcast(&ngenes, 1, MPI_INT, 0, comm);
    MPI_Bcast(&logarithmic_data, 1, MPI_INT, 0, comm);
    MPI_Bcast(&rev_sorting, 1, MPI_INT, 0, comm);
    MPI_Bcast(&rank_sum, 1, MPI_INT, 0, comm);
    MPI_Bcast(&nperms, 1, MPI_INT, 0, comm);

    chunk = nperms/size;
    if ( rank < (nperms - chunk * size) ) {
        ++chunk;
    }

    local_distribution = xmalloc(ngenes * sizeof(double), comm);

    /* We're not allowed to touch the input data since it's just a
     * pointer to R data and so we'd be messing with the ordering
     * when bootstrapping, so copy things into buffers for further
     * use. */
    data1 = xmalloc(ngenes * nclass1 * sizeof(double), comm);

    if ( 0 == rank ) {
        memcpy(data1, data1_, ngenes * nclass1 * sizeof(double));
    }
    MPI_Bcast(data1, ngenes * nclass1, MPI_DOUBLE, 0, comm);

    if ( nclass2 != 0 ) {
        data2 = xmalloc(ngenes * nclass2 * sizeof(double), comm);
        if ( 0 == rank ) {
            memcpy(data2, data2_, ngenes * nclass2 * sizeof(double));
        }
        MPI_Bcast(data2, ngenes * nclass2, MPI_DOUBLE, 0, comm);
    }

    if ( rank != 0 ) {
        experimental_rp = xmalloc(ngenes * sizeof(double), comm);
    }

    MPI_Bcast(experimental_rp, ngenes, MPI_DOUBLE, 0, comm);

    boot_rp(data1, nclass1, data2, nclass2,
            experimental_rp, ngenes,
            logarithmic_data, rev_sorting,
            rank_sum,
            chunk, local_distribution);
    free(data1);
    free(data2);

    gather_boot_data(local_distribution, ngenes, ret, comm);

    free(local_distribution);

    return 0;
}

/* Interface that dispatches to the one- or two-class implementation
 * as appropriate */
int boot_rank_product_multi(int n, ...)
{
    double *data1_ = NULL;
    int *nclass1 = NULL;
    double *data2_ = NULL;
    int *nclass2 = NULL;
    double *experimental_rp = NULL;
    int ngenes;
    int norigins;
    int logarithmic_data;
    int rev_sorting;
    int rank_sum;
    int chunk;
    int nperms;
    double *ret = NULL;
    double *local_distribution = NULL;
    double *data1 = NULL;
    double *data2 = NULL;
    int rank;
    int size;
    MPI_Comm comm;
    va_list ap;
    size_t total_classes;
    int i;

    comm = MPI_COMM_WORLD;

    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    if ( 0 == rank ) {
        va_start(ap, n);

        data1_ = va_arg(ap, double *);
        nclass1 = va_arg(ap, int *);
        data2_ = va_arg(ap, double *);
        nclass2 = va_arg(ap, int *);
        experimental_rp = va_arg(ap, double *);
        ngenes = (int)va_arg(ap, size_t);
        norigins = va_arg(ap, int);
        logarithmic_data = va_arg(ap, int);
        rev_sorting = va_arg(ap, int);
        rank_sum = va_arg(ap, int);
        nperms = va_arg(ap, int);
        ret = va_arg(ap, double *);

        va_end(ap);
    }

    MPI_Bcast(&norigins, 1, MPI_INT, 0, comm);

    if ( 0 != rank ) {
        nclass1 = xmalloc(norigins * sizeof(int), comm);
        nclass2 = xmalloc(norigins * sizeof(int), comm);
    }
    MPI_Bcast(nclass1, norigins, MPI_INT, 0, comm);
    MPI_Bcast(nclass2, norigins, MPI_INT, 0, comm);
    MPI_Bcast(&ngenes, 1, MPI_INT, 0, comm);
    MPI_Bcast(&logarithmic_data, 1, MPI_INT, 0, comm);
    MPI_Bcast(&rev_sorting, 1, MPI_INT, 0, comm);
    MPI_Bcast(&rank_sum, 1, MPI_INT, 0, comm);
    MPI_Bcast(&nperms, 1, MPI_INT, 0, comm);

    chunk = nperms/size;
    if ( rank < (nperms - chunk * size) ) {
        ++chunk;
    }

    local_distribution = xmalloc(ngenes * sizeof(double), comm);

    /* We're not allowed to touch the input data since it's just a
     * pointer to R data and so we'd be messing with the ordering
     * when bootstrapping, so copy things into buffers for further
     * use. */
    total_classes = 0;
    for ( i = 0; i < norigins; i++ ) {
        total_classes += nclass1[i];
    }
    data1 = xmalloc(ngenes * total_classes * sizeof(double), comm);

    if ( 0 == rank ) {
        memcpy(data1, data1_, ngenes * total_classes * sizeof(double));
    }
    MPI_Bcast(data1, ngenes * total_classes, MPI_DOUBLE, 0, comm);

    total_classes = 0;
    for ( i = 0; i < norigins; i++ ) {
        total_classes += nclass2[i];
    }
    if ( total_classes != 0 ) {
        data2 = xmalloc(ngenes * total_classes * sizeof(double), comm);
        if ( 0 == rank ) {
            memcpy(data2, data2_, ngenes * total_classes * sizeof(double));
        }
        MPI_Bcast(data2, ngenes * total_classes, MPI_DOUBLE, 0, comm);
    }

    if ( rank != 0 ) {
        experimental_rp = xmalloc(ngenes * sizeof(double), comm);
    }

    MPI_Bcast(experimental_rp, ngenes, MPI_DOUBLE, 0, comm);

    boot_rp_multi(data1, nclass1, data2, nclass2,
                  experimental_rp, ngenes,
                  norigins,
                  logarithmic_data, rev_sorting,
                  rank_sum,
                  chunk, local_distribution);
    free(data1);
    free(data2);

    gather_boot_data(local_distribution, ngenes, ret, comm);

    free(local_distribution);
    if ( rank != 0 ) {
        free(nclass1);
        free(nclass2);
        free(experimental_rp);
    }
    return 0;
}
