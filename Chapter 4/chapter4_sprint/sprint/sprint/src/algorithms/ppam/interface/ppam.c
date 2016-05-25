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

#include "../../../sprint.h"
#include "../../../functions.h"
#include "ppam.h"

extern int pamedoids(int n,...);

/* ====   Partitioning Around Medoids   ====
   
   Accepts infromation from R. Broadcasts command code and data
   Response is returned back to R. 
*/

void ppam(
  int *_map_file, /* flag eq. 1 data is stored on disk, eq. 0 data is stored in memory */
  double *x, /* distance matrix */ 
  int *n_rows, /* nummber of rows in x */
  int *n_clusters, /* number of clusters */
  char **filename, /* path to a data file, in case data is stored on disk */
  /* below is the list of bswap and cstat attributes */
  int *nsend, /*logical*/ int *nrepr,
  int *nelem, double *radus, double *damer, double *avsyl, double *separ, double *ttsyl,
  double *obj, int *med, int *ncluv, double* clusinf, double *sylinf, int *nisol)
{

  int response,intCode;
  int worldSize, worldRank;

  enum commandCodes commandCode;

  // Check that MPI is initialized
  MPI_Initialized(&response);
  if (response) {
    DEBUG("\nMPI is init'ed in ppam\n");
  } else {
    DEBUG("\nMPI is NOT init'ed in ppam\n");
  }

   // Get size and rank from communicator
  MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
  MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);

  // Broadcast command to other processors
  commandCode = PPAM;
  intCode = (int)commandCode;
  MPI_Bcast(&intCode, 1, MPI_INT, 0, MPI_COMM_WORLD);

  // Call the partitioning around medoids function
  response = pamedoids(19, *_map_file, *n_rows, *n_clusters,
                       *filename, x, nsend, nrepr, nelem, radus, damer, avsyl,
                       separ, ttsyl, obj, med, ncluv, clusinf, sylinf, nisol);

  // TODO: return actual numerical response to be interpreted in ppam.R
  return;
  
}



