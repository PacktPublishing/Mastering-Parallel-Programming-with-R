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

/**
 * Construct a look-up table of all the algorithms provided by
 * the cluster.
 **/

#include <stdio.h>
#include <stdarg.h>
#include "../../functions.h"
#include "../../sprint.h"

/*
 * Declare the various command functions as external
 */

extern int test(int n,...);
//extern int svm_call(int n,...);
extern int correlation(int n,...);
extern int permutation(int n,...);
extern int pamedoids(int n,...);
extern int apply(int n,...);
extern int random_forest_driver(int,...);
extern int boot(int,...);
extern int stringDist(int,...);
extern int init_rng_worker(int n, ...);
extern int reset_rng_worker(int n, ...);
extern int boot_rank_product(int n, ...);
extern int boot_rank_product_multi(int n, ...);
extern int hello(int n, ...);
/**
 * This is a dummy operation which can be used where a command code exists
 * but does not represent a useful function.
 **/

int voidCommand(int n,...)
{
  Rprintf("Void command called, I would not expect this to be called.\n");

  return 1;
}

/**
 * This array of function pointers ties up with the commandCode enumeration
 * found in src/functions.h
 **/

commandFunction commandLUT[] = {voidCommand,
//				svm_call,
                                correlation,
                                permutation,
                                pamedoids,
                                apply,
                                random_forest_driver,
                                boot,
                                stringDist,
                                test,
                                init_rng_worker,
                                reset_rng_worker,
                                boot_rank_product,
                                boot_rank_product_multi,
                                hello,
                                voidCommand};

