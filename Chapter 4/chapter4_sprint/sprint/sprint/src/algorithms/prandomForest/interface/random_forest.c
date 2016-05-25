/**************************************************************************
 *                                                                        *
 *  SPRINT: Simple Parallel R INTerface                                   *
 *  Copyright © 2008,2009,2010 The University of Edinburgh                *
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

#include "../../../sprint.h"
#include "../../../functions.h"
#include <stdarg.h>

extern int random_forest_driver(int n, ...);
SEXP prandomForest(SEXP args)
{
    SEXP result;

    enum commandCodes commandCode;
    int response,intCode;
    MPI_Comm comm;
    MPI_Initialized(&response);

    comm = MPI_COMM_WORLD;
    if ( response ) {
        DEBUG("MPI is initted in random_forest\n");
    } else {
        DEBUG("MPI not initted in random_forest\n");
        return ScalarInteger(-1);
    }

    /* We can run on a single process (the master does some work in
     * our implementation, so there's no need to check the size of
     * MPI_COMM_WORLD. */

    commandCode = PRANDOMFOREST;
    intCode = (int)commandCode;

    MPI_Bcast(&intCode, 1, MPI_INT, 0, comm);

    response = random_forest_driver(0, args, &result);
    return result;
}
