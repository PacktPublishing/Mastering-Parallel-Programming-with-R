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

#include "../../../sprint.h"
#include "../../../functions.h"

extern int permutation(int n,...);

void pmaxT(double *d, int *pnrow, int *pncol, int*L, double *pna, double *T, double *P,
           double *adjP, int *pB, int *index, char **options, int *extra, int *generator_flag) {

    enum commandCodes commandCode;
    int response,intCode;

    // Check that MPI is initialized
    MPI_Initialized(&response);
    if (response) {
        DEBUG("\nMPI is init'ed in pmaxT\n");
    } else {
        DEBUG("\nMPI is NOT init'ed in pmaxT\n");

        // It's a bit weird but I'll set the value of columns to a negative
        // value and the R code will check this
        *pncol = -1;
        return;
    }
    
    // Broadcast command to other processors
    commandCode = PMAXT;
    intCode = (int)commandCode;
    MPI_Bcast(&intCode, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Call the permutation function
    response = permutation(13, d, pnrow, pncol, L, pna, T, P, adjP, pB, index,
                           options, extra, generator_flag);

}

/* ************************************************************************************************************************** *
 *  double  *d                  -->  all input data                                                                           *
 *  int     *pnrow              -->  scalar integer value                                                                     *
 *  int     *pncol              -->  scalar integer value                                                                     *
 *  int     *L                  -->  vector of integer values of size: (1) pncol **OR** (2) pncol/2 --> if(test == "pairt")   *
 *  double  *pna                -->  scalar double value                                                                      *
 *  double  *T                  -->  allocated memory for *all* T values         |                                            *
 *  double  *P                  -->  allocated memory for *all* P values         | Probably I need to                         *
 *  double  *adjP               -->  allocated memory for *all* adjusted values  | remove these initializations               *
 *  int     *pB                 -->  scalar integer value                                                                     *
 *  int     *index              -->  allocated memory for the index                                                           *
 *  char    **options           -->  c(test, side, fixed.seed.sampling)                                                       *
 *                                   ----------------------------------                                                       *
 *                                   test == c("t", "f", "blockf", "pairt", "wilcoxon", "t.equalvar")                         *
 *                                   side == c("upper","abs","lower")                                                         *
 *                                   fixed.seed.sampling == c("y","n")                                                        *
 *  int     *extra              -->  scalar integer value                                                                     *
 *  int     *generator_flag     -->  scalar integer value                                                                     *
 * ************************************************************************************************************************** */

