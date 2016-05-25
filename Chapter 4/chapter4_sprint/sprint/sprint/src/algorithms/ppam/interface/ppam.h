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
#ifndef _PPAM_H
#define _PPAM_H

void ppam(int *_do_swap, double *x, int *nrows, int *kk, char **filename, 
	 int *nsend, /*logical*/ int *nrepr, int *nelem,
	 double *radus, double *damer, double *avsyl, double *separ,
	 double *ttsyl, double *obj, int *med, int *ncluv,
          double *clusinf, double *sylinf, int *nisol);

void bswap(int kk, int nsam, int *nrepr,
	   Rboolean med_given, Rboolean do_swap, int trace_lev,
	   double *dysma, double *dysmb, double *beter,
	   double *dys, double *sky, double s, double *obj);

#endif
