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
#include "../implementation/rank-product.h"
#include "../../../sprint.h"

void rank_product_(const double *data1, const int *nclass1_,
                   const double *data2, const int *nclass2_,
                   const int *ngenes_, const int *logarithmic_data_,
                   const int *rev_sorting_,
                   const int *rank_sum_,
                   double *ret)
{
    int nclass1;
    int nclass2;
    const size_t ngenes = (const size_t)(*ngenes_);
    bool logarithmic_data;
    bool rev_sorting;
    bool rank_sum;
    fold_change_t *buf;
    logarithmic_data = (bool)(*logarithmic_data_);
    rev_sorting = (bool)(*rev_sorting_);
    rank_sum = (bool)(*rank_sum_);

    nclass1 = *nclass1_;
    if ( NULL == nclass2_ ) {
        nclass2 = 0;
    } else {
        nclass2 = *nclass2_;
    }
    /* Pre-allocate to avoid memory churn. */
    buf = xmalloc(sizeof(fold_change_t) * ngenes, MPI_COMM_WORLD);

    rank_product(buf, data1, nclass1, data2, nclass2, ngenes,
                 logarithmic_data, rev_sorting, rank_sum, ret);
    free(buf);
}

void rank_product_multi_(const double *data1, const int *nclass1_,
                         const double *data2, const int *nclass2_,
                         const int *ngenes_, const int *norigins_,
                         const int *logarithmic_data_,
                         const int *rev_sorting_,
                         const int *rank_sum_,
                         double *ret)
{
    int *nclass1;
    int *nclass2;
    int norigins;
    const size_t ngenes = (const size_t)(*ngenes_);
    bool logarithmic_data;
    bool rev_sorting;
    bool rank_sum;
    fold_change_t *buf;
    logarithmic_data = (bool)(*logarithmic_data_);
    rev_sorting = (bool)(*rev_sorting_);
    rank_sum = (bool)(*rank_sum_);
    norigins = *norigins_;

    nclass1 = (int *)nclass1_;
    if ( NULL == nclass2_ ) {
        nclass2 = xcalloc(norigins, sizeof(int), MPI_COMM_WORLD);
    } else {
        nclass2 = (int *)nclass2_;
    }
    /* Pre-allocate to avoid memory churn. */
    buf = xmalloc(sizeof(fold_change_t) * ngenes, MPI_COMM_WORLD);

    rank_product_multi(buf, data1, nclass1, data2, nclass2, ngenes,
                       norigins,
                       logarithmic_data, rev_sorting, rank_sum, ret);
    free(buf);
}
