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
 *                             2                                          *
 *  You should have received a copy of the GNU General Public License     *
 *  along with this program. If not, see <http://www.gnu.or/licenses/>.   *
 *                                                                        *
 **************************************************************************/
#ifndef _APPLY_H
#define _APPLY_H

SEXP matrixApply(
  SEXP result,
  SEXP data,
  SEXP margin,
  SEXP function,
  int worldRank,
  int worldSize);

SEXP listApply(
    SEXP result,
    SEXP data,
    SEXP function,
    int worldRank,
    int worldSize);

SEXP ffApply(SEXP result,
             SEXP data,
             SEXP margin,
             SEXP function,
             SEXP nrows,
             SEXP ncols,
             int worldRank,
             SEXP out_filename,
             int worldSize);

SEXP mappedObjectApply(
  SEXP result,
  double *data,
  SEXP margin,
  SEXP function,
  int nrows,
  int ncols,
  int worldRank,
  int worldSize);

void do_matrixApply(
  SEXP ans,
  SEXP data,
  SEXP margin,
  SEXP function,
  int my_start,
  int my_end,
  int* dimensions,
  int worldRank);

void doMappedObjectApply(
  SEXP ans,
  double *data,
  SEXP margin,
  SEXP function,
  int my_start,
  int my_end,
  int nrows,
  int ncols);

SEXP parseExpression(SEXP expressionSexp);
  
#endif
