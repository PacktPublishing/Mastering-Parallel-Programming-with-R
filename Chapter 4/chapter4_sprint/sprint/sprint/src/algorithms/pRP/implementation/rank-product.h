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

#ifndef _RANK_PRODUCT_H
#define _RANK_PRODUCT_H
#define _GNU_SOURCE
#include <stdlib.h>
#include <stdbool.h>
#include <stdarg.h>
#include <mpi.h>

typedef struct _tagged_fc_t {
    double val;                 /* fold-change value */
    double rank;                /* rank of this gene */
    int idx;                    /* which gene did the fold-change come
                                 * from? */
} fold_change_t;

typedef enum _class_t { CLASS_1 = 0, CLASS_2 } class_t;

int compare_fcs_up(const void *, const void *);
int compare_fcs_down(const void *, const void *);
int compare_doubles(const void *, const void *);

void fold_change(fold_change_t *,
                 const double *, const double *,
                 const size_t);

void log_fold_change(fold_change_t *,
                     const double *, const double *,
                     const size_t);

void rank_fold_changes(fold_change_t *,
                       const double *, const double *,
                       const size_t,
                       const bool, const bool);

void one_class_rp(fold_change_t *,
                  const double *, const int,
                  const size_t,
                  const bool, const bool,
                  const bool, double *);

void two_class_rp(fold_change_t *,
                  const double *, const int,
                  const double *, const int,
                  const size_t, const bool, const bool,
                  const bool, double *);

void rank_product(fold_change_t *,
                  const double *, const int,
                  const double *, const int,
                  const size_t,
                  const bool, const bool,
                  const bool,
                  double *);

void rank_product_multi(fold_change_t *,
                        const double *, const int *,
                        const double *, const int *,
                        const size_t,
                        const int,
                        const bool, const bool,
                        const bool,
                        double *);

void rank_product_(const double *, const int *,
                   const double *, const int *,
                   const int *,
                   const int *,
                   const int *,
                   const int *,
                   double *);

void rank_product_multi_(const double *, const int *,
                         const double *, const int *,
                         const int *,
                         const int *,
                         const int *,
                         const int *,
                         const int *,
                         double *);
void shuffle(double *, const size_t);

void shuffle_data(double *, const int, const size_t);

void shuffle_data_multi(double *, const int *, const int, const size_t);

void boot_rp(double *, const int,
             double *, const int,
             const double *,
             const size_t,
             const bool, const bool,
             const bool,
             const int,
             double *);

void boot_rp_multi(double *, const int *,
                   double *, const int *,
                   const double *,
                   const size_t,
                   const int,
                   const bool,
                   const bool,
                   const bool,
                   const int,
                   double *);

void gather_boot_data(double *, const size_t,
                      double *, MPI_Comm);

int boot_rank_product(int, ...);

int boot_rank_product_multi(int, ...);

void boot_rank_product_(const double *, const int *,
                        const double *, const int *,
                        const double *,
                        const int *,
                        const int *, const int *,
                        const int *,
                        const int *,
                        double *,
                        double *);

void boot_rank_product_multi_(const double *, const int *,
                              const double *, const int *,
                              const double *,
                              const int *,
                              const int *,
                              const int *,
                              const int *,
                              const int *,
                              const int *,
                              double *,
                              double *);

int bsearch_approx(double *, const double, const size_t);

void *xmalloc(const size_t, const MPI_Comm);

void *xcalloc(const size_t, const size_t, const MPI_Comm);
#endif
