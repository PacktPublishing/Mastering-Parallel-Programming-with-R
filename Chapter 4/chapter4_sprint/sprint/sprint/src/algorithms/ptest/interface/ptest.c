/**************************************************************************
 *                                                                        *
 *  SPRINT: Simple Parallel R INTerface                                   *
 *  Copyright © 2008,2009 The University of Edinburgh                     *
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

#include <Rdefines.h>
#include "../../../sprint.h"
#include "../../../functions.h"
extern int test(int n,...);

/* ******************************************************** *
 *  The stub for the R side of a very simple test command.  *
 *  Simply issues the command and returns.                  *
 * ******************************************************** */

SEXP ptest()
{
    SEXP result;
    char **func_results;
    int i, response, worldSize, intCode;
    enum commandCodes commandCode;

    MPI_Initialized(&response);
    if (response) {
        DEBUG("MPI is init'ed in ptest\n");
    } else {

        DEBUG("MPI is NOT init'ed in ptest\n");
        PROTECT(result = NEW_INTEGER(1));
        INTEGER(result)[0] = -1;
        UNPROTECT(1);

        return result;
    }

    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);

    func_results = (char **)malloc(sizeof(char*) * worldSize);

    // broadcast command to other processors
    commandCode = PTEST;
    intCode = (int)commandCode;
    MPI_Bcast(&intCode, 1, MPI_INT, 0, MPI_COMM_WORLD);

    response = test(1, func_results);

    PROTECT(result = allocVector(STRSXP, worldSize));

    for(i=0; i < worldSize; i++) {
        // add message to the response vector
        SET_STRING_ELT(result, i, mkChar(func_results[i]));
        free(func_results[i]);
    }
    free(func_results);

    UNPROTECT(1);
    return result;

}

