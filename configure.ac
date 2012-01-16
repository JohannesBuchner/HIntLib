# Process this file with autoconf to produce a configure script.

###############################################################################
# Header
###############################################################################

AC_INIT(HIntLib, 0.0.3, [rudolf.schuerer@sbg.ac.at] )
AM_INIT_AUTOMAKE(HIntLib, 0.0.3)
AM_CONFIG_HEADER(config.h)

# Information on the package

AC_CONFIG_SRCDIR([src/HIntLib/integrator.h])
AC_CONFIG_FILES([Makefile src/Makefile src/HIntLib/Makefile share/Makefile tests/Makefile tests/toms659/Makefile tests/toms738/Makefile libtest/Makefile])
AC_LANG(C++)

###############################################################################
# Checks for programs
###############################################################################

# C++ Compiler

AC_PROG_CXX

# Fortran Compiler

fortran_check_progs=
fortran_check_dirs=

AC_PROG_F77

if test ! -z "${F77}"; then
   AC_MSG_CHECKING([whether ${F77} really works])
   AC_LANG_PUSH(Fortran 77)
   AC_TRY_COMPILE(,[
      INTEGER I
      I = 1
],
      [fortran_check_progs="compare_with_halton_in_toms659 compare_with_faure_in_toms659 compare_with_sobol_in_toms659 compare_with_nied_2_in_toms738 compare_with_nied_prime_in_toms738 compare_with_nied_pp_in_toms738"
      fortran_check_dirs="toms659 toms738"]
      AC_MSG_RESULT(yes), AC_MSG_RESULT(no))
   AC_LANG_POP(Fortran 77)
fi 

if test -z "$fortran_check_dirs"; then
   AC_MSG_WARN([Fortran Compiler does not work. Some tests in the test-suit will be skipped!])
fi

AC_SUBST(fortran_check_progs)
AC_SUBST(fortran_check_dirs)

# Other programs

AC_PROG_LN_S
AC_PROG_INSTALL

AM_PROG_LIBTOOL
if test -z "${GCC}"; then
AC_MSG_NOTICE([Patching libtool])
   sed -e 's/^with_gcc=$/with_gcc=no/' libtool > libtool.temp &&
   mv libtool.temp libtool
fi

AC_CHECK_PROG(HAVE_TAIL,tail,yes)
AC_CHECK_PROG(HAVE_SORT,sort,yes)

if test x"$HAVE_TAIL" != xyes -o x"$HAVE_SORT" != xyes; then
   AC_MSG_WARN([Either tail or sort is missing. Some programs in the test-suite will be skipped!])
   diff_check_progs=
else
   diff_check_progs=qrnsequence
fi
AC_SUBST(diff_check_progs)


###############################################################################
# Add compiler flags
###############################################################################

# Add MPI_HEADER_PATH to CXXFLAGS

AC_ARG_ENABLE([MPI], AC_HELP_STRING([--disable-MPI], [do not use MPI]),,)

if test "X${MPI_HEADER_PATH}" != "X"; then
  CXXFLAGS="${CXXFLAGS} -I${MPI_HEADER_PATH}"
  AC_MSG_NOTICE([Adding compiler flag -I${MPI_HEADER_PATH}])
fi

# Check for SGI CC and add compiler/linker options

AC_MSG_CHECKING([whether we are using the SGI C++ compiler CC])
AH_TEMPLATE([SGI],
            [Define if we are compiling with the SGI C++ compiler CC.])
AH_TEMPLATE([INSTANTIATE_STL],
       [Define if templates from the STL have to be instantiated manually.])

case $host/$CXX in
*-sgi-irix*/CC | "*-sgi-irix*/CC *" )
   AC_MSG_RESULT(yes)
   CXXFLAGS="${CXXFLAGS} -LANG:std"
   AC_DEFINE(SGI)
   AC_DEFINE(INSTANTIATE_STL)
   AC_MSG_NOTICE([Adding compiler flag -LANG:std])
   ;;
   *) AC_MSG_RESULT(no) ;;
esac

# Check for GNU g++ and ad compiler options

if test "$GXX" = yes; then
   AC_MSG_NOTICE([Adding compiler flags for the GNU C++ compiler g++])
   GNUFLAGS="-ffor-scope  -fstrict-aliasing -Wall -W -Wno-unused -Wctor-dtor-privacy -Wnon-virtual-dtor -Wreorder -Woverloaded-virtual -Wsign-promo"
fi
AC_SUBST(GNUFLAGS)

# A number of GCC-options we do NOT use for various reasons
# -malign-double    # breaks ifstream!!!
# -Wold-style-cast  # Produces many warnings in the STL
# -fast-math        # Problems with NaN handling in EstErr::getRelErr()

###############################################################################
# Checks for libraries.
###############################################################################

AC_CHECK_LIB(m,sqrt)
# AC_CHECK_LIB(mpi,MPI_Init)      # not used yet

###############################################################################
# Checks for header files.
###############################################################################

# Mandatory header files

AC_CHECK_HEADERS([stddef.h stdlib.h string.h \
                  string vector queue algorithm functional exception memory \
                  iomanip iosfwd iostream],,
   AC_MSG_WARN([Headerfile missig. Compilation will probably fail!!!]), [ ])

# Check for cmath of math.h

AC_CHECK_HEADERS([cmath],,
  AC_CHECK_HEADERS([math.h],,
     AC_MSG_WARN([Headerfile missing: We need either cmath or math.h!!!])
  , [ ])
, [ ])

# Check for limits. If it is missing, make sure we have limits.h and float.h

AC_CHECK_HEADERS([limits],,
  AC_CHECK_HEADERS([limits.h],,
    AC_MSG_WARN([Headerfile mising:  We need either limits or limits.h!!!]),[ ])
  AC_CHECK_HEADERS([float.h],,
    AC_MSG_WARN([Headerfile mising:  We need either limits or float.h!!!]),[ ])
, [ ])

# ostream is optional

AC_CHECK_HEADERS([ostream],,,[ ])

# Check for sstream. If its missing, we need streambuf to replace it
   
AC_CHECK_HEADERS([sstream],,
  AC_CHECK_HEADERS([streambuf],,
    AC_CHECK_HEADERS([streambuf.h],,
      AC_MSG_WARN([Headerfile mising:  We need either sstream or streambuf of streambu.h!!!])
    ,[ ])
  ,[ ])
,[ ])

# Do we have MPI?

LIB_HINTLIB_MPI=

if test x"$enable_MPI" != xno; then
AC_CHECK_HEADERS([mpi.h],
  LIB_HINTLIB_MPI=[libhintlibmpi.la],
  AC_MSG_WARN([Cannot find header file mpi.h. No parallel library will be built!
     
     Try     ./configure MPI_HEADER_PATH=path     to solve this problem!
     ]),
  [ ])
fi

AC_SUBST(LIB_HINTLIB_MPI)

###############################################################################
# Checks for types
###############################################################################

# Checks for typedefs, structures

AC_TYPE_SIZE_T
AC_CHECK_TYPES([ptrdiff_t])
AC_CHECK_TYPES([long long int, unsigned long long int])

# Size checks

AC_CHECK_SIZEOF(unsigned int)
AC_CHECK_SIZEOF(unsigned long int)
AC_CHECK_SIZEOF(unsigned long long int)

###############################################################################
# Do compile time  checks
###############################################################################

HL_STREAMS_SUPPORT_LOCAL
HL_FLT_ROUNDS_CONST

###############################################################################
# Do runtime checks
###############################################################################

HL_IEEE_MAGIC_WORKS

HL_ABS_AVAILABLE(std,std)
HL_ABS_AVAILABLE(,global)

AC_CHECK_FUNCS([floor pow sqrt])

AC_PROG_MAKE_SET

AC_DEFINE_UNQUOTED(DATADIR,"${prefix}/share/${PACKAGE}",[The directory which holds our data files])
AC_OUTPUT

