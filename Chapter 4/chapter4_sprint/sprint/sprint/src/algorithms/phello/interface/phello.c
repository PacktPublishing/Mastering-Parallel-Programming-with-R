
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

#include <Rdefines.h>
#include "../../../sprint.h"
#include "../../../functions.h"
extern int hello(int n, ...);
/* ****************************************************************
 *  The stub for the R side of a very simple hello world command.	*
 *  Simply issues the command and returns 0 for successful 	     	*      
 *  completion of command or -1 for failure.                     	*
 * ****************************************************************/
// Note that all data from R is of type SEXP.
SEXP phello()
{
    SEXP result;
    int response, intCode;
    enum commandCodes commandCode;

    // Check MPI initialisaion
    MPI_Initialized(&response);
    if (response) {
        DEBUG("MPI is init'ed in phello\n");
    } else {
        DEBUG("MPI is NOT init'ed in phello\n");
        
        // return -1 if MPI is not initialised.
        PROTECT(result = NEW_INTEGER(1));
        INTEGER(result)[0] = -1;
        UNPROTECT(1);
        return result;
    }

    // broadcast command to other processors
    commandCode = PHELLO;
    intCode = (int)commandCode;
    DEBUG("commandCode in phello is %d \n", intCode);
    MPI_Bcast(&intCode, 1, MPI_INT, 0, MPI_COMM_WORLD);

    response = hello(0); // We are passing no arguments. 
    // If we wanted to pass 2 arguments, we'd write
    // response = hello(2, arg1, arg2);
  
    result =PROTECT(result = NEW_INTEGER(1));
    INTEGER(result)[0] = response;
    UNPROTECT(1);
    return result;

}

