/**************************************************************************
 *                                                                        *
 *  SPRINT: Simple Parallel R INTerface                                   *
 *  Copyright © 2012 The University of Edinburgh                          *
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

#ifndef _SERIALIZE_H
#define _SERIALIZE_H

/* When the serialize.c interface is made available to third party
 * packages this could be pulled entirely into the C layer.  As it is
 * we have to go through the interpreter. */
SEXP serialize_form(SEXP form);
SEXP unserialize_form(SEXP form);
SEXP getListElement(SEXP list, char *str);
void setListElement(SEXP list, char *str, SEXP value);
void reduce_combine(SEXP in, SEXP *out, SEXP (*combine_fn)(SEXP, SEXP),
                    int root, MPI_Comm comm);

#endif

