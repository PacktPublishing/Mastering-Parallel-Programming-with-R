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

#ifndef _SPRINT_H
#define _SPRINT_H

#include <mpi.h>
#include <R.h>
#include <R_ext/Random.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <Rdefines.h>
#include <dlfcn.h>

/* Logging relies on the time library */
#include <sys/time.h>

/* A handy little macro to force logging */
#define LOG(...) \
        { \
          struct timeval tick; \
          gettimeofday(&tick, NULL); \
          Rprintf( "%i.%i - ", (int)tick.tv_sec, (int)tick.tv_usec); \
          Rprintf( __VA_ARGS__);  \
	  }

/* A handy little macto to force error logging */
#define ERR(...) \
        { \
          struct timeval tick; \
          gettimeofday(&tick, NULL); \
          error( "%i.%i - ", (int)tick.tv_sec, (int)tick.tv_usec);  \
          error( __VA_ARGS__); \
        }

/* A handy little macro to force debug logging */
#ifdef _DEBUG
#define DEBUG(...) \
        { \
          struct timeval tick; \
          gettimeofday(&tick, NULL); \
          Rprintf( "%i.%i - ", (int)tick.tv_sec, (int)tick.tv_usec); \
          Rprintf( __VA_ARGS__); \
        }
#else
#define DEBUG(...)
#endif

/* A handy macro to force profiling */
#ifdef _PROF
#define PROF(x, ...) if (!(x)) { Rprintf( __VA_ARGS__);  }
#else
#define PROF(...)
#endif

#endif
