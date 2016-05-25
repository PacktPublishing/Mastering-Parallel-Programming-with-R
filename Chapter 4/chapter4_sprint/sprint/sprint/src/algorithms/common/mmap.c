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
 *                                                                        *
 *  You should have received a copy of the GNU General Public License     *
 *  along with this program. If not, see <http://www.gnu.or/licenses/>.   *
 *                                                                        *
 **************************************************************************/
#include <sys/mman.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include "mmap.h"
#include "../../sprint.h"

double* map_file(char* filename, int* filesize) {

  struct stat sb;
  double *p;
  int fd;

  /* Opens the file to be mapped, reads metadata, maps the file to
     memory and closes the file */
  
  fd = open (filename, O_RDONLY);
  if (fd == -1) {
	  error ("Could not open %s", filename);
    return NULL;
  }
  
  if (fstat (fd, &sb) == -1) {
    error("fstat");
    return NULL;
  }

  if (!S_ISREG (sb.st_mode)) {
    error ("%s is not a file \n", filename);
    return NULL;
  }

  p = mmap (0, sb.st_size, PROT_READ, MAP_SHARED, fd, 0);
  *filesize = sb.st_size;

  if (p == MAP_FAILED) {
    error("mmap");
    return NULL;
  }

  if(close(fd) == -1) {
    error("close");
    return NULL;
  }
  
  return p;
}

void unmap_file (double* p, char* filename) {

  struct stat sb;
  int fd;

  fd = open (filename, O_RDONLY);
  if (fd == -1) {
    error ("open");
    return;
  }

  if (fstat (fd, &sb) == -1) {
    error("fstat");
    return;
  }
                   
  if (munmap (p, sb.st_size) == -1) {
    error ("munmap");
    return;
  }

  if(close(fd) == -1) {
    error("close");
    return;
  }

  return;
}
