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
#include <math.h>
#include "rank-product.h"
#include "../../../sprint.h"

/* Compute the Rank Product for data of a single class.  This will
 * typically be from a two-colour microarray where each entry has
 * already calculated the fold-change internally. */
void one_class_rp(fold_change_t *buf,
                  const double *data, const int nclass,
                  const size_t ngenes,
                  const bool logarithmic_data, const bool rev_sorting,
                  const bool rank_sum,
                  double *ret)
{
    size_t i;
    size_t k;
    /* Don't have to compute pairwise fold-changes, since the data
     * already have this built in.  So just rank the fold-changes from
     * each "gene" vector and compute the rank product. */
    for ( i = 0; i < nclass; i++ ) {
        rank_fold_changes(buf,
                          data + ngenes * i,
                          NULL,
                          ngenes,
                          logarithmic_data,
                          rev_sorting);
        /* The R implementation does a different calculation when the
         * number of samples and/or genes is very small, but both
         * options preserve the monotonic property of the rank
         * product, so we just use this case.  This is more likely to
         * avoid floating point overflow too. */
        if ( rank_sum ) {
            for ( k = 0; k < ngenes; k++ ) {
                ret[buf[k].idx] += buf[k].rank/nclass;
            }
        } else {
            for ( k = 0; k < ngenes; k++ ) {
                ret[buf[k].idx] *= pow(buf[k].rank, 1.0/nclass);
            }
        }
    }
}

/* Compute the Rank product for two class data.  This data will
 * typically be from single-channel microarrays, where some subset of
 * the samples are in class 1 and the rest are in class 2. */
void two_class_rp(fold_change_t *buf,
                  const double *data1, const int nclass1,
                  const double *data2, const int nclass2,
                  const size_t ngenes,
                  const bool logarithmic_data, const bool rev_sorting,
                  const bool rank_sum,
                  double *ret)
{
    size_t i;
    size_t j;
    size_t k;

    /* For each sample in the first class */
    for ( i = 0; i < nclass1; i++ ) {
        /* We compute the fold-change relative to each sample in the
         * second class */
        for ( j = 0; j < nclass2; j++ ) {
            rank_fold_changes(buf,
                              data1 + ngenes * i,
                              data2 + ngenes * j,
                              ngenes,
                              logarithmic_data,
                              rev_sorting);
            if ( rank_sum ) {
                for ( k = 0; k < ngenes; k++ ) {
                    ret[buf[k].idx] += buf[k].rank/(nclass1 * nclass2);
                }
            } else {
                /* And compute the rank product for each gene. */
                for ( k = 0; k < ngenes; k++ ) {
                    ret[buf[k].idx] *= pow(buf[k].rank, 1.0/(nclass1 * nclass2));
                }
            }
        }
    }
}

/* Compute Rank Product for some microarray data.  This just
 * dispatches to the relevant one- or two-class implementation as
 * necessary. */
void rank_product(fold_change_t *buf,
                  const double *data1, const int nclass1,
                  const double *data2, const int nclass2,
                  const size_t ngenes,
                  const bool logarithmic_data,
                  const bool rev_sorting,
                  const bool rank_sum,
                  double *ret)
{
    size_t i;

    for ( i = 0; i < ngenes; i++ ) {
        ret[i] = rank_sum ? 0.0 : 1.0;
    }

    if ( 0 == nclass2 ) {
        one_class_rp(buf, data1, nclass1, ngenes,
                     logarithmic_data, rev_sorting,
                     rank_sum, ret);
    } else {
        two_class_rp(buf, data1, nclass1, data2, nclass2, ngenes,
                     logarithmic_data, rev_sorting, rank_sum, ret);
    }

}

void rank_product_multi(fold_change_t *buf,
                        const double *data1, const int *nclass1,
                        const double *data2, const int *nclass2,
                        const size_t ngenes,
                        const int norigins,
                        const bool logarithmic_data,
                        const bool rev_sorting,
                        const bool rank_sum,
                        double *ret)
{
    size_t base_offset1;
    size_t base_offset2;
    double total_rank;
    double *tmp;
    int i;
    int j;
    for ( i = 0; i < ngenes; i++ ) {
        ret[i] = rank_sum ? 0.0 : 1.0;
    }
    tmp = xmalloc(sizeof(double) * ngenes, MPI_COMM_WORLD);

    total_rank = 0.0;
    for ( i = 0; i < norigins; i++ ) {
        if ( nclass2[i] == 0 ) {
            total_rank += (double)nclass1[i];
        } else {
            total_rank += (double)(nclass1[i] * nclass2[i]);
        }
    }
    base_offset1 = 0;
    base_offset2 = 0;
    for ( i = 0; i < norigins; i++ ) {
        rank_product(buf,
                     data1 + base_offset1,
                     nclass1[i],
                     data2 + base_offset2,
                     nclass2[i],
                     ngenes,
                     logarithmic_data,
                     rev_sorting,
                     rank_sum,
                     tmp);
        base_offset1 += nclass1[i] * ngenes;
        base_offset2 += nclass2[i] * ngenes;
        if ( rank_sum ) {
            for ( j = 0; j < ngenes; j++ ) {
                ret[j] += tmp[j] * (nclass2[i] == 0
                                    ? nclass1[i]
                                    : (nclass1[i] * nclass2[i])) / total_rank;
            }
        } else {
            for ( j = 0; j < ngenes; j++ ) {
                /* We undo the root taken on individual origin and then
                 * take root relative to total number of classes across
                 * all origins */
                ret[j] *= pow(tmp[j],
                              (nclass2[i] == 0
                               ? nclass1[i]
                               : (nclass1[i] * nclass2[i])) / total_rank);
            }
        }
    }
    free(tmp);
}
