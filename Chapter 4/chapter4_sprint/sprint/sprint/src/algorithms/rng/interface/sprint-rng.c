/**************************************************************************
 *                                                                        *
 *  SPRINT: Simple Parallel R INTerface                                   *
 *  Copyright © 2008-2011 The University of Edinburgh                     *
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
 *  along with this program. If not, see <http://www.gnu.org/licenses/>.  *
 *                                                                        *
 **************************************************************************/

#include <stdint.h>
#include <stdarg.h>
#include <sys/time.h>
#include <Rdefines.h>
#include "../../../sprint.h"
#include "../../../functions.h"

extern int init_rng_worker(int n, ...);
extern int reset_rng_worker(int n, ...);

SEXP sprint_create_seed()
{
    SEXP result;
    struct timeval tv;

    PROTECT(result = allocVector(REALSXP, 1));

    gettimeofday (&tv, NULL);
    REAL(result)[0] = (((uint64_t) tv.tv_usec << 16) ^ tv.tv_sec);

    UNPROTECT(1);
    return result;
}

            
SEXP init_rng(SEXP seed)
{
    enum commandCodes command;
    int response,intCode;
    MPI_Comm comm;

    MPI_Initialized(&response);

    comm = MPI_COMM_WORLD;

    if ( response ) {
        DEBUG("MPI is initialized in init_rng\n");
    } else {
        DEBUG("MPI is not initialized in init_rng\n");
        return ScalarInteger(-1);
    }

    command = INIT_RNG;
    intCode = (int)command;

    MPI_Bcast(&intCode, 1, MPI_INT, 0, comm);

    init_rng_worker(1, seed);

    return R_NilValue;
}

SEXP reset_rng()
{
    enum commandCodes command;
    int response,intCode;
    MPI_Comm comm;

    MPI_Initialized(&response);

    comm = MPI_COMM_WORLD;

    if ( response ) {
        DEBUG("MPI is initialized in init_rng\n");
    } else {
        DEBUG("MPI is not initialized in init_rng\n");
        return ScalarInteger(-1);
    }

    command = RESET_RNG;
    intCode = (int)command;

    MPI_Bcast(&intCode, 1, MPI_INT, 0, comm);

    reset_rng_worker(1);

    return R_NilValue;
}
      
    
