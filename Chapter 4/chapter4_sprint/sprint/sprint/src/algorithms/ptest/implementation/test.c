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

#include <stdio.h>
#include <stdarg.h>
#include <mpi.h>
#include <string.h>

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <Rdefines.h>

/**
 * The test command simple prints out rank and communicator size.
 **/

int test(int n,...)
{
  int worldSize;
  int worldRank;
  int i, message_length;
  char* message_buf;
  va_list ap;
  char **in_array;

  MPI_Status stat;

  SEXP response, R_fcall;
  SEXP hello_string, my_rank;

  MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
  MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);

  /* Creat R LANGSXP Vector, R function holder
     length of the vector is 1 + number of arguments */
  PROTECT(R_fcall = allocVector(LANGSXP, 3));
  
  /* Create R objects, argument holders passed to
     a function. */
  PROTECT(hello_string = allocVector(STRSXP, 1));
  PROTECT(my_rank = allocVector(INTSXP, 1));

  /* Initialise first argument */
  SET_STRING_ELT(hello_string, 0, mkChar("HELLO, FROM PROCESSOR:"));
  
  /* Initialise second argument */
  INTEGER(my_rank)[0] = worldRank;

  /* Build function object, first element is a function name
     followed by arguments */
  SETCAR(R_fcall, install("paste"));
  SETCADR(R_fcall, hello_string);
  SETCADDR(R_fcall, my_rank);

  /* Response container, Vector of strings, worldSize
     determines vector length */
  PROTECT(response = allocVector(STRSXP, 1));

  /* Pass the function to R evaluator */
  SET_STRING_ELT(response, 0, STRING_ELT(eval(R_fcall, R_GlobalEnv), 0));

  // master processor gathers results form slaves
  if (worldRank == 0) {

    // Get input variables
    va_start(ap, n);
    in_array = va_arg(ap,char**);
    va_end(ap);

    message_length = strlen(CHAR(STRING_ELT(response, 0))) + 1;
    in_array[0] = (char *)malloc(sizeof(char) * message_length);
    memcpy(in_array[0], CHAR(STRING_ELT(response, worldRank)), message_length);

    for (i=1; i<worldSize; i++) {
      // get the message length
      MPI_Recv(&message_length, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &stat);
      
      in_array[i] = (char *)malloc(sizeof(char) * message_length);

      // get the message
      MPI_Recv(in_array[i], message_length, MPI_CHAR, stat.MPI_SOURCE, 0, MPI_COMM_WORLD, &stat);

    }

    UNPROTECT(4);
    return 0;

    // slaves send their results
  } else {

    // determine and send the message length
    message_length = strlen(CHAR(STRING_ELT(response, 0))) + 1;
    MPI_Send(&message_length, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    
    message_buf = (char *)CHAR(STRING_ELT(response, 0));

    // send message
    MPI_Send(message_buf, message_length, MPI_CHAR, 0, 0, MPI_COMM_WORLD);

    // Slave processes expect an integer for response
    // They check the respopse to pick up errors
    UNPROTECT(4);

    return 0;

  }

}

