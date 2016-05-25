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
#include <R_ext/Print.h>/* for diagnostics */
#include <R_ext/Utils.h>/* for interrupting */

#include "../../common/mmap.h"
#include "../../../sprint.h"
#include "../../common/utils.h"
#include "pamedoids.h"
#include "bswap.h"
#include "cstat.h"
#include "silhouette.h"

int pamedoids(int n,...) {

  va_list ap;
  Rboolean all_stats, med_given, do_swap;

  int worldSize, worldRank, trace_lev, clusinf_dim1,
    n_rows, n_clusters, _map_file, my_start, my_end;

  int local_check = 0, global_check = 0;

  double max_distance, sky;

  int  *filesize, *nsend, *nrepr, *nelem, *med, *ncluv,
    *nisol;
  double *distance_matrix, *radus, *damer, *avsyl, *separ,
    *ttsyl,  *obj, *clusinf, *sylinf;
  char *filename;

  char *name;
  name = (char *) R_alloc(256, sizeof(char));

  /* Transient storage allocation. R will reclaim the memory
     at the end of the call to .C Note that this memory will
     be freed on error or user interrupt.
     See Section 6.1.1 "Writing R Extensions" 
  */
  filename = (char *) R_alloc(256, sizeof(char));
  filesize = (int *) R_alloc(1, sizeof(int));
        
  // Get size and rank from communicator
  MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
  MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);

  if(worldRank == MASTER_PROCESS) {
    
    if (n != 19) {
      DEBUG("rank 0 passed incorrect arguments into pamedoids!");
    }
    // Get input variables
    va_start(ap,n);

    _map_file = va_arg(ap,int);
    n_rows = va_arg(ap,int);
    n_clusters = va_arg(ap,int);
    filename = va_arg(ap, char*);
  }

  MPI_Bcast(&_map_file, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&n_rows, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&n_clusters, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(filename, 256, MPI_CHAR, 0, MPI_COMM_WORLD);

  /* Transient storage allocation. R will reclaim the memory
     at the end of the call to .C */
  if(_map_file) {
    distance_matrix = (double *) R_alloc(1, sizeof(double));
  } else {
    distance_matrix = (double *) R_alloc(n_rows*n_rows, sizeof(double));
  }
  
  if(worldRank == MASTER_PROCESS) {

    /* Master process reads all remaining arguments from the arguments list */
    distance_matrix = va_arg(ap,double*);
    nsend = va_arg(ap,int*);
    nrepr = va_arg(ap,int*);
    nelem = va_arg(ap,int*);
    radus = va_arg(ap,double*);
    damer = va_arg(ap,double*);
    avsyl = va_arg(ap,double*);
    separ = va_arg(ap,double*);
    ttsyl = va_arg(ap,double*);
    obj = va_arg(ap,double*);
    med = va_arg(ap,int*);
    ncluv = va_arg(ap,int*);
    clusinf = va_arg(ap,double*);
    sylinf = va_arg(ap,double*);
    nisol = va_arg(ap,int*);
  
    va_end(ap);   
  } else {

    nsend = (int *) R_alloc(n_rows, sizeof(double));
    nrepr = (int *) R_alloc(n_rows, sizeof(double));
    nelem = (int *) R_alloc(n_rows, sizeof(double));
    radus = (double *) R_alloc(n_rows, sizeof(double));
    damer = (double *) R_alloc(n_rows, sizeof(double));
    avsyl = (double *) R_alloc(n_rows, sizeof(double));
    separ = (double *) R_alloc(n_rows, sizeof(double));
    ttsyl  = (double *) R_alloc(1, sizeof(double));
    obj = (double *) R_alloc(2, sizeof(double));
    med = (int *) R_alloc(n_clusters, sizeof(int));
    ncluv = (int *) R_alloc(n_rows, sizeof(int));
    clusinf = (double *) R_alloc(n_clusters*5, sizeof(double));
    sylinf = (double *) R_alloc(n_rows*4, sizeof(double));
    nisol = (int*) R_alloc(n_clusters, sizeof(int));

    if ( (nsend == NULL) || (nrepr == NULL) || (nelem == NULL) || (radus == NULL) ||
         (damer == NULL) || (avsyl == NULL) || (separ == NULL) || (ttsyl == NULL) ||
         (obj == NULL) || (med == NULL) || (ncluv == NULL) || (clusinf == NULL) ||
         (sylinf == NULL) || (nisol == NULL) ) {
      local_check = 1;
      //ERR("**ERROR** : Memory allocation failed on slave process %d. Aborting.\n", worldRank);
    } else {
      local_check = 0;
    }
  }
  /* Check whether memory was successfully allocated on all worker nodes */
  MPI_Allreduce(&local_check, &global_check, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  /*  No need to free memory if allocation fails on one of the workers
      R_alloc will release it after .Call returns to R */
  if ( global_check != 0 ) return -1;
    
  if(_map_file) {

    if((distance_matrix = map_file(filename, filesize)) == NULL) {
      local_check = -1;
      //ERR("Rank: %d Error: Mapping of file %s failed!\n", worldRank, *filename);   
    }

     /* Check if all processes have successfully mapped the file to memory */
    MPI_Allreduce(&local_check, &global_check, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if ( global_check != 0 ) return -1;
    
  } else {
    /* Input data is kept in memory, distribute it among processes */
    MPI_Bcast(distance_matrix, n_rows*n_rows, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }

  MPI_Bcast(nsend, n_rows, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(nrepr, n_rows, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(nelem, n_rows, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(radus, n_rows, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(damer, n_rows, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(avsyl, n_rows, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(separ, n_rows, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(ttsyl, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(obj, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(med, n_clusters, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(ncluv, n_rows, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(clusinf, n_clusters*5, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(sylinf, n_rows*4, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(nisol, n_clusters, MPI_INT, 0, MPI_COMM_WORLD);

  /* Local variables */
  all_stats = (obj[0] == 0.);   /* if false, only return 'ncluv[]' */
  med_given = (med[0] != 0); /* if true, med[] contain initial medoids */
  do_swap = (nisol[0] != 0);

  /* We've figured out if we should do swapping, so initialise nisol
   * correctly */
  nisol[0] = 0;
  clusinf_dim1 = n_clusters;
  
  /* initialise medoids before broadcasting nrepr array */
  initMedoids(nrepr, med, med_given, n_rows, n_clusters);
  trace_lev = (int) obj[1];

  loopDistribute(worldRank, worldSize, n_rows, &my_start, &my_end);

  max_distance = getMaxDistance(distance_matrix, my_start, my_end, n_rows);
  
  bswap(worldRank, worldSize, my_start, my_end, n_clusters, n_rows, nrepr, med_given,
        do_swap, trace_lev, radus, damer, avsyl, distance_matrix, &sky,
        max_distance, obj);

  /* Compute Clustering & STATs if(all_stats): */
  cstat(worldRank, worldSize, my_start, my_end, n_clusters, n_rows, nsend, nrepr,
        all_stats, radus, damer, avsyl, separ, &max_distance, distance_matrix,
        ncluv, nelem, med, nisol);
  
  if(all_stats) {
    
    if(worldRank == MASTER_PROCESS) {
      for (int k = 0; k < n_clusters; k++) {
        
        clusinf[k]=		(double)       nrepr[k];
        clusinf[k + clusinf_dim1]	     = radus[k];
        clusinf[k + (clusinf_dim1 << 1)] = avsyl[k];
        clusinf[k + clusinf_dim1 * 3]    = damer[k];
        clusinf[k + (clusinf_dim1 << 2)] = separ[k];
      }
    }

    if (1 < n_clusters && n_clusters < n_rows) {
    
	    /* Compute Silhouette info : */
      silhouette(worldRank, worldSize, my_start, my_end,
                 n_clusters, n_rows, ncluv, nsend, nelem, nrepr,
                 radus, damer, avsyl, ttsyl, distance_matrix,
                 &max_distance, sylinf);
    }
  }
     
  if(trace_lev) Rprintf("end{cstat()}\n");
  
  return 0;
}

