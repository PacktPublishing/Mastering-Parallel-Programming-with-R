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

void silhouette(int my_rank, int world_size, int my_start, int my_end,
                int n_clusters, int n_rows, int *ncluv,
                int *nsend, int *nelem, int *negbr,
                double *syl, double *srank, double *avsyl, double *ttsyl,
                double *distance_matrix, double *s, double *sylinf)
{
    int k, nsylr;
    
    /* pointers to sylinf[] columns -- sylinf[nn, 4] : */
    double *sylinf_2, *sylinf_3, *sylinf_4;
    sylinf_2 = sylinf	+ n_rows;
    sylinf_3 = sylinf_2 + n_rows;
    sylinf_4 = sylinf_3 + n_rows;

    /* Parameter adjustments */
    //--avsyl;
    //--ncluv;

    nsylr = 0;
    *ttsyl = 0.;
    for (k = 0; k < n_clusters; k++) {
	/* nelem[0:(ntt-1)] := indices (1-based) of obs. in cluster k : */
      int j,l, ntt = 0;
      for (j = 0; j < n_rows; j++) {
        if (ncluv[j] == k+1) {
          nelem[ntt] = j;
          ++ntt;
        }
      }

      int my_ntt_start, my_ntt_end;
      loopDistribute(my_rank, world_size, ntt, &my_ntt_start, &my_ntt_end);

	for (j = my_ntt_start; j < my_ntt_end; j++) {/* (j+1)-th obs. in cluster k */
	    int k_, nj = nelem[j];
	    double dysb = *s * 1.1 + 1.;
	    negbr[j] = -1;
	    /* for all clusters  k_ != k : */
	    for (k_ = 0; k_ < n_clusters; k_++) if (k_ != k) {
		double db = 0.;
		int nbb = 0;
		for (l = 0; l < n_rows; l++) if (ncluv[l] == k_+1) {
		    ++nbb;
		    if (l != nj)          
          db += distance_matrix[DIST_INDEX(nj,l,n_rows)];
		}
		db /= nbb; /* now  db(k_) := mean( d[j, l]; l in C_{k_} ) */
		if (dysb > db) {
		    dysb = db;
		    negbr[j] = k_;
		}
	    }/* negbr[j] := arg max_{k_} db(k_) */
	    if (ntt > 1) {
		double dysa = 0.;
		for (l = 0; l < ntt; ++l) {
		    int nl = nelem[l];
		    if (nj != nl)
			dysa += distance_matrix[DIST_INDEX(nj,nl,n_rows)];
		}
		dysa /= ntt - 1;
		if (dysa > 0.) {
		    if (dysb > 0.) {
			if (dysb > dysa)
			    syl[j] = 1. - dysa / dysb;
			else if (dysb < dysa)
			    syl[j] = dysb / dysa - 1.;
			else /* dysb == dysa: */
			    syl[j] = 0.;

			if (syl[j] < -1.)
			    syl[j] = -1.;
			else if (syl[j] > 1.)
			    syl[j] = 1.;

		    } else {
			syl[j] = -1.;
		    }
		}
		else /* dysa == 0 */ if (dysb > 0.)
		    syl[j] = 1.;
		else
		    syl[j] = 0.;
	    }
	    else { /*	  ntt == 1: */
		syl[j] = 0.;
	    }
	} /* for( j ) */

  mMPI_AllgatherDOUBLE(syl, world_size, my_rank, my_ntt_start, my_ntt_end,
                   ntt, MPI_COMM_WORLD);
  mMPI_AllgatherINT(negbr, world_size, my_rank, my_ntt_start, my_ntt_end,
                   ntt, MPI_COMM_WORLD);

  avsyl[k] = 0.;
	if (ntt == 0) /* this can happen when medoids are user-specified !*/
	    continue; /* next k */

	for (j = 0; j < ntt; ++j) {
	    int lang=-1 /*Wall*/;
	    double symax = -2.;
	    for (l = 0; l < ntt; ++l) {
		if (symax < syl[l]) {
		    symax = syl[l];
		    lang = l;
		}
	    }
	    nsend[j] = lang;
	    srank[j] = symax; /* = syl[lang] */
	    avsyl[k] += srank[j];
	    syl[lang] = -3.;
	}
	*ttsyl += avsyl[k];
	avsyl[k] /= ntt;

  
  
	if (ntt == 1) {
	    sylinf  [nsylr] = (double) k;
	    sylinf_2[nsylr] = (double) negbr[0];
	    sylinf_3[nsylr] = 0.;
	    sylinf_4[nsylr] = (double) nelem[0];
	    ++nsylr;
	} else {
	    for (j = 0; j < ntt; ++j) {
		int lplac = nsend[j];
		sylinf	[nsylr] = (double) k;
		sylinf_2[nsylr] = (double) negbr[lplac];
		sylinf_3[nsylr] = srank[j];
		sylinf_4[nsylr] = (double) nelem[lplac];
		++nsylr;
	    }
	}
    }

/* for (k) */
    *ttsyl /= n_rows;
}
