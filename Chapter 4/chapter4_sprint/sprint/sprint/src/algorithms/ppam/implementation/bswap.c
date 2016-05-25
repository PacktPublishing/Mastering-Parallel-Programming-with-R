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
#include <R.h>

#include "../../../sprint.h"
#include "../../common/utils.h"
#include "bswap.h"

// Include the time.h library
#include <time.h>


/*--------------------------------------------------------------------------
  
  dysma - array contains a distance between an object and the nearest medoid
  s - greater than max(dys[])

  --------------------------------------------------------------------------*/

/*
     bswap(): the clustering algorithm in 2 parts:  I. build,	II. swap
*/
void bswap(int my_rank, int world_size, int my_start, int my_end, int n_clusters,
           int n_rows, int *nrepr,
           Rboolean med_given, Rboolean do_swap, int trace_lev,
           /* nrepr[]: here is boolean (0/1): 1 = "is representative object"  */
           double *dysma, double *dysmb, double *beter,
           double *distance_matrix, double *sky, double s, double *obj)
{

  int i, j, k,h;
  double sky_local;

  int *medoid_ids;
  double *dysma_local;

  /* Parameter adjustments */
  
  /* [MP]: I cannot comprehend why this was done here
  --nrepr;
  --beter;
    
  --dysma; --dysmb; */

  /* Zeros the memory allocated */
  medoid_ids = (int *) S_alloc(n_clusters, sizeof(int));
  dysma_local = (double *) S_alloc(n_rows, sizeof(double));

  if(trace_lev) Rprintf("pam()'s bswap(*, s=%g): ", s);
    
  s = s * 1.1 + 1.;/* larger than all dys[];
                      replacing by DBL_MAX  changes result - why ? */

  // =========  CALCULATE dysma IF MEDOIDS ARE GIVEN ON INPUT ======== /*
  if(med_given) {

    for (i = my_start; i < my_end; i++)
      dysma_local[i] = s;

    if(trace_lev) Rprintf("medoids given\n");
    
    /* Get medoid ids and store it in an array */
    getMedoidIDs(nrepr, medoid_ids, n_rows, n_clusters);
      
    /* compute dysma[] : dysma[j] = D(j, nearest_representative) */
    for (i = my_start; i < my_end; i++) {
      for (j = 0; j < n_clusters; j++) {
        if (dysma_local[i] > distance_matrix[DIST_INDEX(i,medoid_ids[j],n_rows)])
          dysma_local[i] = distance_matrix[DIST_INDEX(i,medoid_ids[j],n_rows)];
      }
    }

    MPI_Allreduce(dysma_local, dysma, n_rows, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  // =========  END  ======== /*
    
  } else {
      
    /* ======== BUILD MEDOIDS ============= */
       
    /* find kk representatives aka medoids :  */

     // Record the start time bswap:
    int nmax_local = -1; 
    double ammax_local, cmd;
    
    int *nmax;
    double *ammax;

    nmax = (int *) R_alloc(world_size, sizeof(int));
    ammax = (double *) R_alloc(world_size, sizeof(double));
    
    /* dysma - array contains a distance between an object and the nearest medoid */
    
    /* set all dysma values to s > max_distance */
    for (i = 0; i < n_rows; i++)
      dysma_local[i] = s;
    
    for (k = 0; k < n_clusters; k++) {
      
      R_CheckUserInterrupt();
      
      /* compute beter[i] for all non-representatives:
       * also find ammax := max_{..} and nmax := argmax_i{beter[i]} ... */
      nmax_local = -1;       
      ammax_local = 0.;
      
      for (i = my_start; i < my_end; i++) {
        if (nrepr[i] == 0) {
          beter[i] = 0.;
          for (j = 0; j < n_rows; j++) {
            cmd = dysma_local[j] - distance_matrix[DIST_INDEX(i,j,n_rows)];
            if (cmd > 0.)
              beter[i] += cmd;
          }
          if (ammax_local <= beter[i]) {
            /*  does < (instead of <= ) work too? -- NO! */
            ammax_local = beter[i];
            nmax_local = i;
          }
        }
      }
      
      //Gather the local maximums
      MPI_Allgather(&nmax_local, 1, MPI_INT, nmax, 1, MPI_INT, MPI_COMM_WORLD);
      MPI_Allgather(&ammax_local, 1, MPI_DOUBLE, ammax, 1, MPI_DOUBLE, MPI_COMM_WORLD);
      mMPI_AllgatherDOUBLE(beter, world_size, my_rank, my_start, my_end,
                           n_rows, MPI_COMM_WORLD);
      
      nmax_local = getMaximum(nmax, ammax, world_size); 
      
      nrepr[nmax_local] = 1; /* = .true. : found new representative */
      if (trace_lev >= 2 && my_rank == 0)
        Rprintf("    new repr. %d\n", nmax_local);
      
          /* update dysma[] : dysma[j] = D(j, nearest_representative) */
      for (j = my_start; j < my_end; j++) {
        if (dysma_local[j] > distance_matrix[DIST_INDEX(nmax_local,j,n_rows)])
          dysma_local[j] = distance_matrix[DIST_INDEX(nmax_local,j,n_rows)];       
      }
      
      mMPI_AllgatherDOUBLE(dysma_local, world_size, my_rank, my_start, my_end,
                           n_rows, MPI_COMM_WORLD);
      
    }   
  }

  if(trace_lev && my_rank == 0) {
    Rprintf("  after build: medoids are");
    for (i = 1; i <= n_rows; ++i)
	    if(nrepr[i] == 1) Rprintf(" %2d", i);
    Rprintf(" \n\n", i);
  }

  // Calculate global sum
  *sky = 0.; sky_local = 0;
  for (j = my_start; j < my_end; j++)
    sky_local += dysma_local[j];

  MPI_Allreduce(&sky_local, sky, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  obj[0] = *sky / n_rows;
  
  if (do_swap && (n_clusters > 1 || med_given)) {

    double dzsky, dys_ij, dys_hj;
    int hbest = -1, nbest = -1, indx;/* init: -Wall*/

    double *dzsky_global;
    int *hbest_global, *nbest_global;

    hbest_global = (int *) R_alloc(world_size, sizeof(int));
    nbest_global = (int *) R_alloc(world_size, sizeof(int));
    dzsky_global = (double *) R_alloc(world_size, sizeof(double));

    /* ====== second algorithm: SWAP. ====== */

    /* Hmm: In the following, we RE-compute dysma[];
     *      don't need it first time; then only need *update* after swap */

/*--   Loop : */
L60:
    for (j = my_start; j < my_end; j++) {
      /*  dysma[j] := D_j  d(j, <closest medi>)  [KR p.102, 104]
       *  dysmb[j] := E_j  d(j, <2-nd cl.medi>)  [p.103] */
      
      dysma[j] = s;
      dysmb[j] = s;
      for (i = 0; i < n_rows; i++) {
        if (nrepr[i]) {
          dys_ij = distance_matrix[DIST_INDEX(j,i,n_rows)];
          
          if (dysma[j] > dys_ij) {
            dysmb[j] = dysma[j];
            dysma[j] = dys_ij;
            
          } else if (dysmb[j] > dys_ij) {
            dysmb[j] = dys_ij;
          }
        }
      }
    }

    mMPI_AllgatherDOUBLE(dysma, world_size, my_rank, my_start, my_end,
                   n_rows, MPI_COMM_WORLD);
    mMPI_AllgatherDOUBLE(dysmb, world_size, my_rank, my_start, my_end,
                   n_rows, MPI_COMM_WORLD);
    
    dzsky = 1.; /* 1 is arbitrary > 0; only dzsky < 0 matters in the end */

    for (h = my_start; h < my_end; h++) if (!nrepr[h]) {
         
        for (i = 0; i < n_rows; i++) if (nrepr[i]) {
            
            double dz = 0.;
            
            for (j = 0; j < n_rows; j++) { /* if (!nrepr[j]) { */
              
              dys_hj = distance_matrix[DIST_INDEX(h,j,n_rows)];
              dys_ij = distance_matrix[DIST_INDEX(i,j,n_rows)];
              
              if (dys_ij == dysma[j]) {
                
                double small = dysmb[j] > dys_hj ? dys_hj : dysmb[j];
                dz += (- dysma[j] + small);
                
              } else if (dys_hj < dysma[j]) /* 1c. */
                dz += (- dysma[j] + dys_hj);
              
            }
            
            if (dzsky > dz) {
              dzsky = dz; /* dzsky := min_{i,h} T_{i,h} */
              hbest = h;
              nbest = i;
              
            }
        }
    }

    MPI_Allgather(&hbest, 1, MPI_INT, hbest_global, 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Allgather(&nbest, 1, MPI_INT, nbest_global, 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Allgather(&dzsky, 1, MPI_DOUBLE, dzsky_global, 1, MPI_DOUBLE, MPI_COMM_WORLD);

    indx = -1;
    indx = getMinimumIndx(dzsky_global, world_size);
      
    dzsky = dzsky_global[indx];
    hbest = hbest_global[indx];
    nbest = nbest_global[indx];

    if (dzsky < 0.) { /* found an improving swap */

      if(trace_lev >= 2)
        Rprintf( "   swp new %d <-> %d old; decreasing diss. by %g\n", hbest, nbest, dzsky);
      
      nrepr[hbest] = 1;
      nrepr[nbest] = 0;
      *sky += dzsky;
      
      goto L60;
    }
    
  }
  obj[1] = *sky / n_rows;

}

