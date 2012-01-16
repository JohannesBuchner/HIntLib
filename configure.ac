# Process this file with autoconf to produce a configure script.

###############################################################################
# Header
###############################################################################

# autoconf setup

AC_INIT([HIntLib], 0.0.6, [rudolf.schuerer@sbg.ac.at], [HIntLib])
AC_CONFIG_SRCDIR([src/HIntLib/integrator.h])
AC_CONFIG_HEADERS(config.h)
AC_CONFIG_FILES([Makefile src/Makefile src/HIntLib/Makefile share/Makefile tests/Makefile tests/toms659/Makefile tests/toms738/Makefile libtest/Makefile])
dnl AC_CONFIG_FILES([Makefile src/Makefile src/HIntLib/Makefile share/Makefile tests/Makefile tests/toms659/Makefile tests/toms738/Makefile libtest/Makefile development/Makefile])
AC_LANG(C++)

# automake setup

AC_MSG_NOTICE([Initializing automake:])
AM_INIT_AUTOMAKE

# Information on the package

AC_CANONICAL_HOST

###############################################################################
# Checks for programs
###############################################################################

# C++ compiler

AC_MSG_NOTICE([Initializing C++ compiler:])
AC_PROG_CXX

# Fortran Compiler

AC_MSG_NOTICE([Initializing Fortran compiler:])
HL_PROG_F77_REALLY_WORKS

AM_CONDITIONAL(HAVE_FORTRAN,test "x${hl_cv_prog_f77_really_works}" = xyes)
if test "x${hl_cv_prog_f77_really_works}" = xno; then
  AC_MSG_WARN([Fortran Compiler does not work. Some tests in the test-suit will be skipped!])
fi

# C compiler
# Wo do not need a C compiler for building HIntLib.  However, LibTool requires
# it for some reason, so we can check for it right here anyway.

AC_MSG_NOTICE([Initializing C compiler:])
AC_PROG_CC

# LibTool

AC_MSG_NOTICE([Initializing LibTool:])
AC_DISABLE_STATIC
AC_PROG_LIBTOOL

# Other programs

AC_MSG_NOTICE([Checking for other programs])
AC_CHECK_PROG(HAVE_TAIL,tail,yes)
AC_CHECK_PROG(HAVE_SORT,sort,yes)

AM_CONDITIONAL(HAVE_TAIL_SORT,test "x$HAVE_TAIL" = xyes -a "x$HAVE_SORT" = xyes)
if test x"$HAVE_TAIL" != xyes -o x"$HAVE_SORT" != xyes; then
   AC_MSG_WARN([Either tail or sort is missing. Some programs in the test-suite will be skipped!])
fi


###############################################################################
# Add compiler flags
###############################################################################

# Add MPI_HEADER_PATH to CXXFLAGS

AC_ARG_ENABLE([MPI], AC_HELP_STRING([--disable-MPI], [do not use MPI]),,)

if test "X${MPI_HEADER_PATH}" != "X"; then
  CXXFLAGS="${CXXFLAGS} -I${MPI_HEADER_PATH}"
  CFLAGS="${CFLAGS} -I${MPI_HEADER_PATH}"
  AC_MSG_NOTICE([Adding compiler flag -I${MPI_HEADER_PATH}])
fi

# Check for SGI CC and add compiler/linker options

AC_MSG_CHECKING([whether we are using the SGI C++ compiler CC])
AH_TEMPLATE([SGI],
            [Define if we are compiling with the SGI C++ compiler CC.])

STL_OBJECT=

case $host/$CXX in
*-sgi-irix*/CC | *-sgi-irix*/CC\ * )
   AC_MSG_RESULT(yes)
   CXXFLAGS="${CXXFLAGS} -LANG:std -no_prelink -ptused -FE:eliminate_duplicate_inline_copies -FE:template_in_elf_section"
   AC_MSG_NOTICE([Adding compiler flags -LANG:std -no_prelink -ptused])
   if test "$LD" = "CC"; then
      AC_MSG_NOTICE([Adding linker flag -no_prelink])
      LDFLAGS="${LDFLAGS} -no_prelink -FE:eliminate_duplicate_inline_copies -FE:template_in_elf_section"
   fi
   AC_DEFINE(SGI)
# STL_OBJECT=stl.lo
   ;;
*) AC_MSG_RESULT(no) ;;
esac

AC_SUBST(STL_OBJECT)

# Cygwin stuff

NOUNDEFINED=

AC_MSG_CHECKING([whether we use Cygwin])

case $host in
*cygwin* )
   AC_MSG_RESULT(yes)
   NOUNDEFINED=-no-undefined
   ;;
*) AC_MSG_RESULT(no) ;;
esac

AC_SUBST(NOUNDEFINED)

# Turn on compiler warnings

AC_MSG_CHECKING([compiler flags to turn on warnings])

if test "$GXX" = yes; then
   COMPILERWARNINGS="-ffor-scope -Wall -W -Wno-unused -Wctor-dtor-privacy -Wnon-virtual-dtor -Wreorder -Woverloaded-virtual -Wsign-promo"
# A number of GCC-options we do NOT use for various reasons
# -malign-double    # breaks ifstream!!!
# -Wold-style-cast  # Produces many warnings in the STL
# -fast-math        # Problems with NaN handling in EstErr::getRelErr()

else
  case $host/$CXX in
*-sgi-irix*/CC | *-sgi-irix*/CC\ * )
   COMPILERWARNINGS="-fullwarn -woff 3439,3626,3649,3625,1424,1468,1209,1375,1210,1174,1506"
   ;;
*-linux-*/icpc | *-linux-*/icpc\ * )
   # COMPILERWARNINGS="-w1"
   # COMPILERWARNINGS="-w2 -Wall -wd981,383,810,530,279,488"
   COMPILERWARNINGS="-w2 -Wall -wd981,383,810,530,279,488,1418,1419"
   ;;
*)
   COMPILERWARNINGS=""
   ;;
esac
fi

AC_MSG_RESULT($COMPILERWARNINGS)
AC_SUBST(COMPILERWARNINGS)

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

# Check for cmath or math.h

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
      AC_MSG_WARN([Headerfile mising:  We need either sstream or streambuf or streambuf.h!!!])
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
# Checks for functions
###############################################################################

AC_CHECK_FUNCS([floor modf pow sqrt])

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

# --with-index=XXX

hl_index=auto

AC_ARG_WITH([index],
  AC_HELP_STRING([--with-index=XX],
     [determines the type of Index. Valid values are 32 (use u32), 64 (use u64), or auto (use u64 if unsigned long is at least 64 bits wide) [[auto]]]),
  AC_MSG_CHECKING([with-index argument])
  hl_index="$withval"
  AC_MSG_RESULT($hl_index)
,)

if test x"$hl_index" = x32; then
  AC_MSG_NOTICE([Index is set to u32.])
else
if test x"$hl_index" = x64; then
  AC_MSG_NOTICE([Index is set to u64.])
else
if test x"$hl_index" = xauto; then
  if test $ac_cv_sizeof_unsigned_long_int -ge 8 ; then
    hl_index=64
    AC_MSG_NOTICE([Index is set to u64 automatically.])
  else
    hl_index=32
    AC_MSG_NOTICE([Index is set to u32 automatically.])
  fi
else
   AC_MSG_ERROR([invalid argument --with-index=$hl_index])
fi
fi
fi
   
AC_DEFINE_UNQUOTED(INDEX,$hl_index,[Size of Index data type.])

# --with-real=XXX

hl_real=2

AC_ARG_WITH([real],
  AC_HELP_STRING([--with-real=XX],
     [determines the type of real. Valid values are float, double, or long (uses long double) [[double]]]),
  AC_MSG_CHECKING([with-real argument])

  if test x"$withval" = xfloat; then
    AC_MSG_RESULT(float)
    hl_real=1
  else
  if test x"$withval" = xdouble; then
    AC_MSG_RESULT(double)
    hl_real=2
  else
  if test x"$withval" = xlong; then
    AC_MSG_RESULT(long double)
    hl_real=3
  else
    AC_MSG_RESULT($withval)
    AC_MSG_ERROR([invalid argument --with-real=$withval])
  fi
  fi
  fi
,)

AC_DEFINE_UNQUOTED(REAL,$hl_real,[Type of real data type.])

###############################################################################
# Do compile time checks
###############################################################################

HL_STREAMS_SUPPORT_LOCAL
HL_FLT_ROUNDS_CONST
HL_NAMESPACE_REL_OPS
HL_UNREACHABLE_CALLS_REMOVED

###############################################################################
# Do runtime checks
###############################################################################

HL_IEEE_MAGIC_WORKS

HL_ABS_AVAILABLE(std,std)
HL_ABS_AVAILABLE(,global)

if test "x$prefix" != xNONE; then   
   AC_DEFINE_UNQUOTED(DATADIR,"${prefix}/share/${PACKAGE}",[The directory which holds our data files])
else
   AC_DEFINE_UNQUOTED(DATADIR,"${ac_default_prefix}/share/${PACKAGE}",[The directory which holds our data files])
fi

AC_DEFINE_UNQUOTED(SRCDATADIR,"`pwd`/share",[The directory which holds the uninstalled data files])

AC_OUTPUT

