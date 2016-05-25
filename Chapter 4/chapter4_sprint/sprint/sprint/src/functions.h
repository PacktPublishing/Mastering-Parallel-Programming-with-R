/**************************************************************************
 *                                                                        *
 *  SPRINT: Simple Parallel R INTerface                                   *
 *  Copyright ? 2008,2009 The University of Edinburgh                     *
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

#ifndef _COMMANDS_H
#define _COMMANDS_H

/**
 * Lists all the functions available, ensure that TERMINATE is first and
 * LAST is last. If you add a command code you must add a command function
 * in sprint/functions.c
 **/

//enum commandCodes {TERMINATE = 0, PSVM, PCOR, PMAXT, PPAM, PAPPLY, PRANDOMFOREST, PBOOT, PSTRINGDIST, PTEST, INIT_RNG, RESET_RNG, PBOOTRP, PBOOTRPMULTI, LAST};
enum commandCodes {TERMINATE = 0, PCOR, PMAXT, PPAM, PAPPLY, PRANDOMFOREST, PBOOT, PSTRINGDIST, PTEST, INIT_RNG, RESET_RNG, PBOOTRP, PBOOTRPMULTI, PHELLO, LAST};

/**
 * Stereotype for interface functions. You almost certainly don't need to
 * mess with this.
 **/

typedef int (*commandFunction)(int n,...);

#endif
