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

#include <R.h>

#include "../../../sprint.h"
#include "comms.h"

void gatherData(SEXP gathered_result, SEXP ans, int N, int my_start, int my_end, int worldRank)
{

  MPI_Status status;
  SEXP recv_buf;
  /* The length of all messages should be identical on all nodes */
  int message_length;
  
  if(worldRank == 0) {

    /* Depends if called by apply or lapply, we are sending back chunks of array or list elements */
    for(int i=my_end; i<N; i++) {

      MPI_Recv(&message_length, 1, MPI_INT, MPI_ANY_SOURCE, i, MPI_COMM_WORLD, &status);
      PROTECT(recv_buf = allocVector(REALSXP, message_length));
      MPI_Recv(REAL(recv_buf), message_length, MPI_DOUBLE, MPI_ANY_SOURCE, i, MPI_COMM_WORLD, &status);
      
      DEBUG("worldRank: %d i: %d, message_length: %d recv_buf: %f\n", worldRank, i, message_length, REAL(recv_buf)[0]);
      
      SET_VECTOR_ELT(ans, i, recv_buf);
      UNPROTECT(1);
    }
  } else if(worldRank != 0) {
    
    for(int i=0, tag=my_start; i<my_end-my_start; i++, tag++) {

      message_length = length(VECTOR_ELT(ans,i));
      
      MPI_Send(&message_length, 1, MPI_INT, 0, tag, MPI_COMM_WORLD);
      MPI_Send(REAL(VECTOR_ELT(ans,i)), message_length, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
    }
  }
  
  if(worldRank==0) {
    SET_VECTOR_ELT(gathered_result, 0, ans);
  }
  return;
}

/* Bcast function name or definition, cover case when definition is split into
   several lines and stored as a SEXP string vector */

void bcastRFunction(SEXP function, int function_nlines, int worldRank) {

  int i, line_length;
  char *msg_buffer;

  for(i=0; i<function_nlines; i++) {

    /* determine the length of the CHARSXP object in STRSXP vector
       and broadcast this value to all workers */
    line_length = 0;
    
    if (worldRank == MASTER_PROCESS) {
      /* add one to include string termination symbol */
      line_length = strlen((char *)CHAR((STRING_ELT(function,i)))) + 1;
    }
    
    MPI_Bcast(&line_length, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    /* allocate a message buffer to send the line of text, can't work on
       the STRSXP object directly as CHARSXP can be modified only by a setter
       function */
    msg_buffer = (char *) malloc(line_length);
    
    if (worldRank == MASTER_PROCESS) {
      strncpy(msg_buffer, (char *)CHAR((STRING_ELT(function,i))), line_length);
      msg_buffer[line_length-1] = '\0';
    }
    
    /* broadcast a line of function definition / name */
    MPI_Bcast(msg_buffer, line_length, MPI_CHAR, 0, MPI_COMM_WORLD);
    
    /* set the dunction STRSXP element on worker nodes */
    if (worldRank != MASTER_PROCESS) {
      SET_STRING_ELT(function, i, mkChar(msg_buffer));
    }
    
    /* free msg_buffer for the next iteration */
    free(msg_buffer);
  }
  return;
}
  
