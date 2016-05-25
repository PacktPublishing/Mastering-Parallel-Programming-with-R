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
#include "cstat.h"


/* -----------------------------------------------------------
 cstat(): Compute STATistics (numerical output) concerning each partition
*/
void cstat(int my_rank, int world_size, int my_start, int my_end, int n_clusters,
           int n_rows, int *nsend, int *nrepr, Rboolean all_stats,
           double *radus, double *damer, double *avsyl, double *separ, double *s,
           double *distance_matrix, int *ncluv, int *nelem, int *med, int *nisol)
{
  int j, k, ja, jk, nplac, ksmal = -1/* -Wall */;
  double ss = *s * 1.1 + 1.;

     int my_kk_start, my_kk_end;

     loopDistribute(my_rank, world_size, n_clusters, &my_kk_start, &my_kk_end);
       
    /* Zeros the memory allocated */
    //nsend_local = (int *) S_alloc(n_rows, sizeof(int));
    //dysma_local = (double *) S_alloc(n_rows, sizeof(double));

    /* Parameter adjustments [MP]: why he was doing that?uflen
    --nisol;
    --med;
    --nelem;
    --ncluv;
    --separ;
    --avsyl;
    --damer;
    --radus;
    --nrepr;
    --nsend; */

    /* nsend - holds index of the closest medoid */
    /* ncluv - holds information to which cluster given object belongd to */
    /* nsend[j] := i,  where x[i,] is the medoid to which x[j,] belongs */

    for (j = my_start; j < my_end; j++) {
      /* if you are not medoid */
      if (nrepr[j] == 0) {
        double dsmal = ss;
        /* find the closest medoid */ 
        for (k = 0; k < n_rows; k++) {
          if (nrepr[k] == 1) {
            if (dsmal > distance_matrix[DIST_INDEX(j,k,n_rows)]) {
              dsmal = distance_matrix[DIST_INDEX(j,k,n_rows)];
              ksmal = k;
            }
          }
        }
        nsend[j] = ksmal;
      } else {
        nsend[j] = j;
      }
    }
    
    mMPI_AllgatherINT(nsend, world_size, my_rank, my_start, my_end,
                      n_rows, MPI_COMM_WORLD);
  
    /* ncluv[j] := k , the cluster number  (k = 1..*kk) */
    jk = 1;
    nplac = nsend[0];
    for (j = my_start; j < my_end; j++) {
      ncluv[j] = 0;
      if (nsend[j] == nplac)
        ncluv[j] = 1;
    }

    mMPI_AllgatherINT(ncluv, world_size, my_rank, my_start, my_end,
                   n_rows, MPI_COMM_WORLD);  

    if(my_rank == MASTER_PROCESS)
      my_start = 1;
  
   /*==========================================*/
    
    for (ja = 1; ja < n_rows; ja++) {
      nplac = nsend[ja];
      if (ncluv[nplac] == 0) {
        ++jk;

        for (j = my_start; j < my_end; j++) {
          if (nsend[j] == nplac)
            ncluv[j] = jk;
        }

        if(my_rank == MASTER_PROCESS)
          my_start = 0;

        mMPI_AllgatherINT(ncluv, world_size, my_rank, my_start, my_end,
                          n_rows, MPI_COMM_WORLD);

        if(my_rank == MASTER_PROCESS)
          my_start = 1;

         if (jk == n_clusters)
          break;
      }
    }

    /*==========================================*/
      
    if(my_rank == MASTER_PROCESS)
      my_start = 0;
    
    if(all_stats) { /*	   analysis of the clustering. */
      int numl;

      int ntt, m;
      double ttt;
      
      for (k = my_kk_start; k < my_kk_end; k++) {
  
        ntt = -1, m = -1;
        ttt = -0.;
        radus[k] = -1.;
        R_CheckUserInterrupt();
 
        for (j = 0; j < n_rows; j++) {
          if (ncluv[j] == k+1) {
            double djm;
            ++ntt;
            m = nsend[j];
            djm = distance_matrix[DIST_INDEX(m,j,n_rows)];
            ttt += djm;

            if (radus[k] < djm)
              radus[k] = djm;
          }
        }

        avsyl[k] = ttt / (ntt+1);
        med[k] = m;
      }     

      mMPI_AllgatherDOUBLE(radus, world_size, my_rank, my_kk_start, my_kk_end,
                   n_clusters, MPI_COMM_WORLD);
      mMPI_AllgatherDOUBLE(avsyl, world_size, my_rank, my_kk_start, my_kk_end,
                   n_clusters, MPI_COMM_WORLD);
      mMPI_AllgatherINT(med, world_size, my_rank, my_kk_start, my_kk_end,
                   n_clusters, MPI_COMM_WORLD);
  
      if (n_clusters == 1) {
        damer[1] = *s;
        nrepr[1] = n_rows;
        return;
      }
     
      /*  ELSE	  kk > 1 : */
      
      /* numl = number of L-clusters. */
      numl = 0;

      /* Memory of previous process iteration is required */

      int nel = 0; 
          
      for (k = my_kk_start; k < my_kk_end; k++) {
        /*
          identification of cluster k:
          nelem= vector of object indices,
          nel  = number of objects
        */
        R_CheckUserInterrupt();

        nel = 0;      
        for (j = 0; j < n_rows; j++) {
          if (ncluv[j] == k+1) { 
            nelem[nel] = j;
            nel++;
          }
        }
              
        nrepr[k] = nel;
        
        if (nel == 1) {
          int nvn = nelem[1];
          damer[k] = 0.;
          separ[k] = ss;
          for (j = 0; j < n_rows; j++) {
            if (j != nvn) {
              if (separ[k] > distance_matrix[DIST_INDEX(nvn,j,n_rows)])
                separ[k] = distance_matrix[DIST_INDEX(nvn,j,n_rows)];
            }
          }
          
          /* Is cluster k
             1) an L-cluster	 or
             2) an L*-cluster ? */
          if (separ[k] == 0.)
            ++numl;
          
        } else { /*	       nel != 1 : */
        
          double dam = -1., sep = ss;
          Rboolean kand = TRUE;
          for (ja = 0; ja < nel; ja++) {
            int jb, nvna = nelem[ja];
            double aja = -1., ajb = ss;
            for (jb = 0; jb < n_rows; jb++) {			
              if (ncluv[jb] == k+1) {
                if (aja < distance_matrix[DIST_INDEX(nvna,jb,n_rows)])
                  aja = distance_matrix[DIST_INDEX(nvna,jb,n_rows)];
              } else {
                if (ajb > distance_matrix[DIST_INDEX(nvna,jb,n_rows)])
                  ajb = distance_matrix[DIST_INDEX(nvna,jb,n_rows)];
              }
            }
            if (kand && aja >= ajb)
              kand = FALSE;
            if (dam < aja)
              dam = aja;
            if (sep > ajb)
              sep = ajb;
          }
          separ[k] = sep;
          damer[k] = dam;
          if (kand) {
            ++numl;
            if (dam >= sep) /*	L-cluster */
              nisol[k] = 1;
            else/*		L*-cluster */
              nisol[k] = 2;
            continue /* k */;
          }
        }
        /* nel = 1 or (!kand) : */
        //nisol[k] = 0;
        
      } /* for(k) */

      mMPI_AllgatherDOUBLE(damer, world_size, my_rank, my_kk_start, my_kk_end,
                   n_clusters, MPI_COMM_WORLD);
      mMPI_AllgatherDOUBLE(separ, world_size, my_rank, my_kk_start, my_kk_end,
                   n_clusters, MPI_COMM_WORLD);
      mMPI_AllgatherINT(nisol, world_size, my_rank, my_kk_start, my_kk_end,
                   n_clusters, MPI_COMM_WORLD);
      mMPI_AllgatherINT(nrepr, world_size, my_rank, my_kk_start, my_kk_end,
                   n_clusters, MPI_COMM_WORLD);

    } /* all_stats */

} /* cstat */
