 Process this file with autoconf to produce a configure script.

dnl * ************************************************ *
dnl *  As short configure script just to find out      *
dnl *  which MPI compiler needs to be used to compile  *
dnl *  and link the MPI libraries to SPRINT.           *
dnl * ************************************************ *

AC_INIT(DESCRIPTION)

dnl * *************************************************** *
dnl *  Scans the system for all known MPI compilers and   *
dnl *  sets the values of $MPICC and $MPILIBS.            *
dnl * *************************************************** *

AC_EGREP_HEADER()

dnl * ********************************************************************* *
dnl * Taken from the Writing R Extensions Guide                             * 
dnl * (http://cran.r-project.org/doc/manuals/R-exts.html#Configure-example) *
dnl * ********************************************************************* *

: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi
CC=`"${R_HOME}/bin/R" CMD config CC`
R_CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
#CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`
SHLIB_LDFLAGS=`"${R_HOME}/bin/R" CMD config SHLIB_LDFLAGS`




AC_MSG_CHECKING([for C MPI support])

ACX_MPI([AC_MSG_RESULT([MPI found with success])],
        [AC_MSG_FAILURE("Unable to detect MPI compiler. Please use --with-wrapper-script option")])


dnl * *************************************************** *
dnl *  Optional parameter to the configure script         *
dnl *  in case the MPI compiler is not one of the         *
dnl *  supported compilers the ACX_MPI macro is scanning  *
dnl *  for. Rather unlikely though...                     *
dnl * *************************************************** *

AC_ARG_WITH(wrapper-script,
    [  --with-wrapper-script=<..>      Uses the cross compiler specified to
                                       compile the package ],
    [  if [[ "${withval}" != no ]];
       then
           AC_PROG_CC()
           MPICC=${withval}
           MPILIBS=""
       fi]
)

CC=$MPICC

AC_MSG_CHECKING([for C MPI-2 support])

AC_RUN_IFELSE(
            [AC_LANG_PROGRAM(
[[
#include <mpi.h>
]],
[[

if (MPI_VERSION > 1) exit(0);
else exit(1);

]])],
            [AC_MSG_RESULT(yes)],
            [AC_MSG_RESULT(no)
             AC_MSG_FAILURE([The SPRINT package requires full C MPI-2 or MPI-3 support. Try a more up to date MPI implementation (MPICH).])])



dnl AC_MSG_CHECKING([for correct MPI implementation])
dnl
dnl AC_RUN_IFELSE(
dnl             [AC_LANG_PROGRAM(
dnl [[
dnl #include <mpi.h>
dnl ]],
dnl [[
dnl 
dnl #ifdef OMPI_RELEASE_VERSION
dnl exit(1);
dnl #endif
dnl exit(0);
dnl ]])],
dnl             [AC_MSG_RESULT(yes)],
dnl             [AC_MSG_RESULT(no)
dnl             AC_MSG_FAILURE([The SPRINT package is not compatible with OpenMPI. Please install MPICH instead as described in the SPRINT User Guide ])])

dnl * ********************************************************* *
dnl *  Set the language to C++ and search for the MPI compiler  *
dnl * ********************************************************* *

AC_LANG_PUSH(C++)
AC_MSG_CHECKING([for C++ MPI support])
ACX_MPI([AC_MSG_RESULT([MPI found with success])],
        [AC_MSG_FAILURE("Unable to detect MPI compiler. Please use --with-wrapper-script option")])
CXX=$MPICXX
AC_MSG_CHECKING([for C++ MPI-2 support])

AC_RUN_IFELSE(
            [AC_LANG_PROGRAM(
[[
#include <mpi.h>
#include <stdlib.h>
]],
[[

if (MPI_VERSION > 1) exit(0);
else exit(1);

]])],
            [AC_MSG_RESULT(yes)],
            [AC_MSG_RESULT(no)
             AC_MSG_FAILURE([The SPRINT package requires full C++ MPI-2 support. Try a more up to date MPI implementation (MPICH2 or openMPI).])])

AC_LANG_POP()

dnl * *************************************************** *
dnl * Check for grep -E and egrep, in that order,         *
dnl * and set output variable EGREP to the first one      *
dnl * that is found.                                      *
dnl * *************************************************** *

dnl AC_PROG_EGREP()


dnl * *************************************************** *
dnl * Scan mpi.h header file for information about MPI-2  *
dnl * support. For some reason generic AC_EGREP_HEADER	  *
dnl * macro does not work for me. Also autoconf does not  *
dnl * support usuall grep escape strings.		  *
dnl * *************************************************** *

dnl AC_EGREP_HEADER(MPI_VERSION, mpi.h,
dnl		[AC_MSG_RESULT([yes])],
dnl 		[AC_MSG_FAILURE([no mpi-2 detected])])

dnl AC_MSG_CHECKING([for MPI_HOME])
dnl
dnl    if env | $EGREP -q 'MPI_HOME=.'
dnl    then 
dnl	AC_MSG_RESULT([yes])
dnl    else
dnl	AC_MSG_RESULT([no])
dnl	AC_MSG_ERROR([MPI_HOME env variable is not set. Please set the path of the MPI installation to use.])
dnl    fi

dnl AC_MSG_CHECKING([for MPI-2 support])

dnl    if $EGREP 'define MPI_VERSION' $MPI_HOME/include/mpi.h | grep -q '2'
dnl     then 
dnl	AC_MSG_RESULT([yes])
dnl    else
dnl	AC_MSG_RESULT([no])
dnl	AC_MSG_ERROR([This MPI implementation does not support the full MPI-2 functionality.])
dnl    fi

dnl * ************************************************ *
dnl *  If the MPI compiler is "cc" then that is the    *
dnl *  way to identify that we are in a Cray machine   *
dnl *  and we need to add the MPI library path to the  *
dnl *  MPILIBS variable.			       *
dnl * ************************************************ *

dnl * Commenting out as I don't believe this is necessary on a Cray anymore.
dnl if [[ "$MPICC" == cc ]]; then
dnl    MPILIBS="-L$MPICH_DIR/lib -lmpich"
dnl fi


dnl * ************************************************ *
dnl *  If it's the GCC CNU compiler used then add the  *
dnl *  "gnu99" C lang standar flag.                    *
dnl * ************************************************ *
AC_PROG_CC_C99
dnl if [[ "$GCC" == yes ]]; then
dnl    MPICC="$MPICC -std=gnu99 "
dnl fi
dnl determine if we are running in solaris
dnl AC_CONFIG_AUX_DIR([tools])
dnl AC_CANONICAL_HOST

dnl if [[ $host == *-*-solaris* ]]; then
dnl   MPICC="$MPICC -xc99 "
dnl fi

dnl * ************************************************ *
dnl *  Replacing these lines from makefile             *
dnl * ************************************************ *
dnl ifdef MPICXX
dnl 	$(MPICXX) $(SHLIB_LDFLAGS) -o $@ $^
dnl else
dnl 	$(CC) $(SHLIB_LDFLAGS) -o $@ $^ $(MPILIBS)
dnl endif

if [[ "$MPICXX" == "" ]]; then
    MPICXX=$(CC)
fi


dnl * ********************************************** *
dnl *  Substitute the flags to the Makevars.in file  *
dnl * ********************************************** *

AC_SUBST(MPILIBS)
AC_SUBST(CC)
AC_SUBST(MPICC)
AC_SUBST(MPICXX)
AC_SUBST(R_CFLAGS)
#AC_SUBST(CPPFLAGS)
AC_SUBST(SHLIB_LDFLAGS)

AC_OUTPUT(src/Makevars)

