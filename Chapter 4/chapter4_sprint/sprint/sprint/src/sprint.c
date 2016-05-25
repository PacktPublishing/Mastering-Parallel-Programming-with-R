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

/* ****************************************************** *
 *  The compute cluster framework.                        *
 *                                                        *
 *  This is compiled into a shared library which is then  *
 *  linked in to R.                                       *
 * ****************************************************** */

#include "sprint.h"
#include "functions.h"
#include <R.h>
#include <R_ext/Rdynload.h>
#include <Rinterface.h>
#include <Rembedded.h>

extern commandFunction commandLUT[];
static int mpi_init_flag = -1;

SEXP worker();

/* *********************************************************** *
 *  Initialise MPI environment. "Borrowed" from Rmpi (Hao Yu)  *
 * *********************************************************** */
void R_init_sprint(DllInfo *Dllinfo) {

    int worldSize, worldRank;
    int flag;
    char *dlhandle;
    static int fake_argc = 1;
    char *fake_argv[1];
    char *fake_argv0 = "R";
    int response;
    
  /* #work around for openMPI on some Linux versions. If this fails to work there will be dll errors.*/
#ifdef OPEN_MPI
  #if defined(__unix__) 
        dlhandle = dlopen("libmpi.so", RTLD_GLOBAL | RTLD_LAZY);
        if ( NULL == dlhandle ) {
            ERR("Failed to open libmpi.so library. %s\n", dlerror());
            return;
        }
  #endif
#endif

    MPI_Initialized(&flag);
    if (flag) {
        ERR("Flag: %d\n",flag);
        return;
    }
    else {
        fake_argv[0] = (char *)&fake_argv0;
        MPI_Init_thread(&fake_argc, (char ***)(void*)&fake_argv , MPI_THREAD_SERIALIZED, &response); 
        //MPI_Init(&fake_argc, (char ***)(void*)&fake_argv);
        DEBUG("%i: Starting up\n", worldRank);
        mpi_init_flag = 1;

        MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
        MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);
        if ( worldRank == 0 ) {
            MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
            MPI_Errhandler_set(MPI_COMM_SELF, MPI_ERRORS_RETURN);
        }
        return;
    }
}

/* *********************************************************************** *
 *                              Shutdown function                          *
 *                              -----------------                          *
 *  Called only by the master thread when pterminate() function is called  *
 * *********************************************************************** */
SEXP sprint_shutdown() {
    int worldRank,intCode;
    enum commandCodes commandCode;

    MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);
    if (mpi_init_flag == 1) {

        // Send terminate command
        commandCode = TERMINATE;
        intCode = (int)commandCode;
        DEBUG("%i: Shutting down %i\n", worldRank, commandCode);
        MPI_Bcast(&intCode, 1, MPI_INT, 0, MPI_COMM_WORLD);

        // Sprint should now do the finalize
        MPI_Finalize();

        // Set the flag
        mpi_init_flag = -1;
    }

    return ScalarInteger(1);
}

/* *************************************************************** *
 *                       Worker function                           *
 *                       ---------------                           *
 *  The node waits for a command file, reads its contents, acts    *
 *  upon that command and creates a result file. These files are   *
 *  concerned only with data for the compute cluster skeleton,     *
 *  they do not pass information needed by individual algorithms.  *
 *  They must sort that out themselves.                            *
 * *************************************************************** */

SEXP worker() {

    enum commandCodes commandCode;
    int intCode; // Need to pass int in MPI_Bcast instead of enum
    int response;
    int worldRank;
    int worldSize;

    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);

    if ( worldRank == 0 ) {
        return ScalarInteger(worldRank);
    }
    /* Start the command processing loop */

    do {

        // this matches the broadcast in interface/algorithm.c which is executed on rank 0
        // via R. All the rest of the processors should wait here

        intCode = (int)commandCode;
        MPI_Bcast(&intCode, 1, MPI_INT, 0, MPI_COMM_WORLD);
        commandCode = intCode;

        DEBUG("%i: Receiving command %i\n", worldRank, commandCode);
        /* Process the command */

        if (commandCode == TERMINATE) {
            /* TERMINATE is a special case. It is handled here rather than
             * via the look-up table.
             */
            DEBUG("%i: Worker shutting down.\n", worldRank);
        } else if ((commandCode > TERMINATE) && (commandCode < LAST)) {
            /* If the command code is in the range defined by the enum then
             * it is a real command so we use the look-up table to find out
             * which function to call.
             */

            DEBUG("%i: Performing command %i\n", worldRank, commandCode);
            response = commandLUT[commandCode](0);
            if (response != 0) {
                ERR("%i: Function returned ERR %i\n", worldRank, response);
            }
        } else {
            /* If we reach this point then it was a bum command and we should
             * not try to process it.
             */
            ERR("%i: Unrecognised command %i\n", worldRank, commandCode);
        }

    } while (commandCode != TERMINATE);

    /* Shut down MPI */
    MPI_Finalize();

    DEBUG("%i: End logging\n", worldRank);

    /* There used to be a call to R_CleanTempDir() here, but the temp dirs are
    cleaned automatically by R without this call (which fails the CRAN checks). */

    // Return the MPI rank of the process. This value is needed in R
    // to identify the 'worker' MPI processes and handle their shutdown.
    return ScalarInteger(worldRank);
}

