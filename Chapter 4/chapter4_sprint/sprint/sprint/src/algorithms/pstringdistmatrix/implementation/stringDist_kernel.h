/**************************************************************************
 *                                                                        *
 *  SPRINT: Simple Parallel R INTerface                                   *
 *  Copyright © 2012 The University of Edinburgh                          *
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

#ifndef _STRINGDIST_KERNEL_H
#define _STRINGDIST_KERNEL_H

int computeStringDist(int worldRank, int worldSize, char* DNAStringSet, int *stringDist, char *out_filename,
                   int sample_width, int number_of_samples, int my_start, int my_end, int chunk_size);
int allocateMaxChunk(int worldRank, int my_start, int my_end, int *stringDist,
                     int number_of_samples);

#endif
