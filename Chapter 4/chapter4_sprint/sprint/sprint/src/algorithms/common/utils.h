/**************************************************************************
 *                                                                        *
 *  SPRINT: Simple Parallel R INTerface                                   *
 *  Copyright © 2010 The University of Edinburgh                          *
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

#ifndef _UTILS_H
#define _UTILS_H

#define DIST_INDEX(i,j,n) (i*n + j)
#define MASTER_PROCESS 0

void getMedoidIDs(int *nrepr, int *medoid_ids, int n_rows, int n_clusters);
void mMPI_AllgatherDOUBLE(double* send_buf, int world_size, int my_rank, int my_start,
                    int my_end, int n_rows, MPI_Comm comm);
void mMPI_AllgatherINT(int* send_buf, int world_size, int my_rank, int my_start,
                       int my_end, int n_rows, MPI_Comm comm);
int getMaximum(int *nmax, double *ammax, int world_size);
int getMinimumIndx(double *array, int array_length);
double getMaxDistance(double *dys, int my_start, int my_end, int nn);
void loopDistribute(int myid, int num_of_proc, int N, int *my_start, int *my_end);
void initMedoids(int *nrepr, int *med, Rboolean med_given, int nn, int n_clusters);
SEXP parseExpression(SEXP expressionSexp);
#endif
