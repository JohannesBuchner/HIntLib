# Process this file with autoconf to produce a configure script.

###############################################################################
# Header
###############################################################################

# autoconf setup

AC_INIT([HIntLib], 0.0.11, [rudolf.schuerer@sbg.ac.at], [HIntLib])

dnl SHARED_VERSION_NUMBER is the libtool library version used in connection
dnl with the -version-info option.
dnl
dnl 0.0. 9: 1
dnl 0.0.10: 2
dnl 0.0.11: 3
SHARED_VERSION_NUMBER=3
AC_SUBST(SHARED_VERSION_NUMBER)
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

case $host in
*cygwin* )
   if test $enable_static = yes ; then
      if test $enable_shared = yes ; then
        AC_MSG_WARN([

-------------------------------------------------------------------------------
You have configured HIntLib on a Cygwin system such that a static library as
well as a shared library (DLL) is built.
Whenever you compile an object file which is supposed to be linked against the
static HIntLib library you have to add '-DHINTLIB_STATIC_LIB_ONLY=1' to the
command line of the compiler.
If you want to avoid this, configure HIntLib to use either only a static
library (--disable-shared --endable-static) or only a DLL (default).
-------------------------------------------------------------------------------
])
      else
        AC_DEFINE(STATIC_LIB_ONLY,1,[Define to 1 if only a static library is built on Cygwin.])
      fi
   fi
   ;;
esac

# Turn on compiler warnings

AC_MSG_CHECKING([compiler flags to turn on warnings])

case $host/$CXX in
*/g++*)
   COMPILERWARNINGS="-ffor-scope -ftemplate-depth-25 -Wall -W -Wno-unused -Wctor-dtor-privacy -Wnon-virtual-dtor -Wreorder -Woverloaded-virtual -Wsign-promo"
# A number of GCC-options we do NOT use for various reasons
# -malign-double    # breaks ifstream!!!
# -Wold-style-cast  # Produces many warnings in the STL
# -fast-math        # Problems with NaN handling in EstErr::getRelErr()
   ;;
*-sgi-irix*/CC | *-sgi-irix*/CC\ * )
   COMPILERWARNINGS="-fullwarn -woff 3439,3626,3649,3625,1424,1468,1209,1375,1210,1174,1506"
   ;;
*-linux*/icpc | *-linux*/icpc\ * )
   # COMPILERWARNINGS="-w1"
   # COMPILERWARNINGS="-w2 -Wall -wd981,383,810,530,279,488"
   COMPILERWARNINGS="-w2 -Wall -wd981,383,810,530,279,488,1418,1419,444,1683"
   ;;
*)
   COMPILERWARNINGS=""
   ;;
esac

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

AC_MSG_NOTICE([Checking for header files:])

# Mandatory C++ and Unix header files

dnl unistd.h:  getopt() and relatives

AC_CHECK_HEADERS(
   [string vector queue algorithm utility functional exception memory complex \
    iomanip iosfwd iostream fstream unistd.h],,
   AC_MSG_ERROR([Header file missig!!!]), [ ])

# Mandatory C header files

dnl  cmath:   sin(), cos(), sqrt(), floor(), ceil(), exp(), log(), pow(), modf()
dnl  math.h:  lgamma()
dnl  cstring: strcpy()
dnl  cstddef: size_t, ptrdiff_t
dnl  cstdlib: rand(), srand(), RAND_MAX

HL_CHECK_C_OR_CPP_HEADER([math])
HL_CHECK_C_OR_CPP_HEADER([string])
HL_CHECK_C_OR_CPP_HEADER([stddef])
HL_CHECK_C_OR_CPP_HEADER([stdlib])

# limits

AC_CHECK_HEADERS([limits],,
  AC_CHECK_HEADERS([limits.h float.h],,
    AC_MSG_ERROR([Header file mising:  We need either limits or limits.h and float.h!!!]),[ ]), [ ])

# Check for sstream. If its missing, we need streambuf to replace it
   
hl_found=no
AC_CHECK_HEADERS([sstream streambuf streambuf.h],hl_found=yes;break,,[ ])

if test $hl_found = no; then
   AC_MSG_ERROR([Header file mising:  We need either sstream or streambuf or streambuf.h!!!])
fi

# Optional header files: ostream, istream, locale

AC_CHECK_HEADERS([ostream],,,[ ])
AC_CHECK_HEADERS([istream],,,[ ])
AC_CHECK_HEADERS([locale],,,[ ])

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

dnl AC_CHECK_FUNCS([floor modf pow sqrt])

HL_TEST_BUILTIN(parity)
HL_TEST_BUILTIN(parityl)
HL_TEST_BUILTIN(parityll)
HL_TEST_BUILTIN(popcount)
HL_TEST_BUILTIN(popcountl)
HL_TEST_BUILTIN(popcountll)
HL_TEST_BUILTIN(ffs)
HL_TEST_BUILTIN(ffsl)
HL_TEST_BUILTIN(ffsll)
HL_TEST_BUILTIN(clz)
HL_TEST_BUILTIN(clzl)
HL_TEST_BUILTIN(clzll)
HL_TEST_BUILTIN(ctz)
HL_TEST_BUILTIN(ctzl)
HL_TEST_BUILTIN(ctzll)

###############################################################################
# Checks for types
###############################################################################

# Checks for typedefs, structures

AC_TYPE_SIZE_T
AC_CHECK_TYPES([ptrdiff_t,long long int,unsigned long long int,long double,wchar_t])

# Size checks

AC_CHECK_SIZEOF(unsigned int)
AC_CHECK_SIZEOF(unsigned long int)
AC_CHECK_SIZEOF(unsigned long long int)
AC_CHECK_SIZEOF(wchar_t)

# --with-index=ARG

AC_ARG_WITH([index],
  AC_HELP_STRING([--with-index=ARG],
     [determines the type of Index. Valid values for ARG are `32' (use u32), `64' (use u64), and `auto' (use u64 if unsigned long is at least 64 bits wide) [[default=auto]]]),
  AC_MSG_CHECKING([with-index argument])
  hl_index="$withval"
  AC_MSG_RESULT($hl_index)
,hl_index=auto)

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
   
AC_DEFINE_UNQUOTED(INDEX,$hl_index,[Size of `Index' type.])

# --with-real=ARG

AC_ARG_WITH([real],
  AC_HELP_STRING([--with-real=ARG],
     [determines the type of real. Valid values for ARG are `float', `double', or `long' (uses long double) [[default=double]]]),
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
, hl_real=2
)

AC_DEFINE_UNQUOTED(REAL,$hl_real,[Type of `real' type.])

###############################################################################
# Do compile time and run time checks
###############################################################################

HL_STREAMS_SUPPORT_LOCALE
HL_FLT_ROUNDS_CONST
HL_UNREACHABLE_CALLS_REMOVED
HL_EQUAL_BUG
if test "x$ac_cv_type_long_double" = xyes; then
HL_COMPLEX_POW_BUG
fi
HL_IEEE_MAGIC_WORKS

# Determine wether we use wchar_t and UNICODE

hl_build_wchar=1

AC_ARG_ENABLE(wchar,
AS_HELP_STRING([--enable-wchar],[Turns on support for wchar_t with UNICODE, even if this capability is not detected])
AS_HELP_STRING([--disable-wchar],[Turns off support for wchar_t, even if this capability is detected]),[
  if test x"$enableval" != xyes -a x"$enableval" != xno; then
    AC_MSG_ERROR([invalid argument to --enable-wchar=$enableval])
  else
    AC_MSG_CHECKING([whether support for wchar_t and UNICODE should be built])
    hl_build_wchar=$enableval
    AC_MSG_RESULT([(--enable-wchar) $hl_build_wchar])
  fi],
  HL_UNICODE
  HL_OSTREAM_IS_BASIC_OSTREAM
  AC_MSG_CHECKING([whether support for wchar_t and UNICODE should be built])
  hl_build_wchar=yes
  test x"$hl_cv_unicode" == xyes || hl_build_wchar=no
  test x"$hl_cv_ostream_is_basic_ostream" == xyes || hl_build_wchar=no
  AC_MSG_RESULT([(autodetect) $hl_build_wchar])
)

if test x"$hl_build_wchar" = xyes ; then
  AC_DEFINE_UNQUOTED(BUILD_WCHAR,1,[Define to 1 if support for wchar_t should be built.])
fi

# --with-encoding-for-char=ARG

AC_MSG_CHECKING([which encoding should be used for char])
AC_ARG_WITH([encoding-for-char],
  AC_HELP_STRING([--with-encoding-for-char=ARG],
     [determines the character encoding used for char.  Valid values for ARG are `source' (only characters from the source character set), `latin1' (for Latin-1 or ISO-8859-1), `utf8' (for UTF-8 encoding), and `locale' (for switching between Latin-1 and UTF-8 depending on the locale of the stream at runtime).  The default value is `locale' if streams support locale, `latin1' for Cygwin, and `source' otherwise]),
  [case x"$withval" in
    xsource) AC_MSG_RESULT([only source character set])
             ;;
    xlatin1) hl_encoding_for_char=2
             AC_MSG_RESULT([Latin-1 or ISO-8859-1])
             ;;
    xutf8)   AC_DEFINE(ENCODING_UTF8,1,[Define to 1 if char uses UTF-8 encoding.])
             AC_MSG_RESULT([UTF-8])
             ;;
    xlocale) hl_encoding_for_char=4
             AC_MSG_RESULT([determined by the locale of the stream])
             if test x"$hl_cv_streams_support_locale" != xyes; then
                AC_MSG_ERROR([Invalid argument --with-encoding-for-char=locale because streams do not support locale!])
             fi
             ;;
    *) AC_MSG_ERROR([invalid argument --with-encoding-for-char=$withval]) ;;
  esac],[
  case $host in
    *cygwin* )
      hl_encoding_for_char=2
      AC_MSG_RESULT([Latin-1 or ISO-8859-1 (default for Cygwin)])
      ;;
    *)
      if test x"$hl_cv_streams_support_locale" = xyes; then
         hl_encoding_for_char=4
         AC_MSG_RESULT([determined by the locale of the stream (default)])
      else
         hl_encoding_for_char=1
         AC_MSG_RESULT([only source character set (default, because streams do not support locale)])
      fi
      ;;
  esac]
)

if test x"$hl_encoding_for_char" == x2 ; then
  AC_DEFINE(ENCODING_LATIN1,1,[Define to 1 if char uses Latin-1 encoding.])
fi
if test x"$hl_encoding_for_char" == x4 ; then
  AC_DEFINE(ENCODING_LOCALE,1,[Define to 1 if encoding of char should be determined by the locale.])
fi


# --with-character-set=ARG

AC_MSG_CHECKING([which set of characters should be used])
AC_ARG_WITH([character-set],
  AC_HELP_STRING([--with-character-set=ARG],
     [determines the set of characters used by HIntLib (independently of the encoding).  Valid values for ARG are `ascii' (for ASCII or ISO 646), `latin1' (for Latin-1 or ISO-8859-1), wgl4 (for all characters from the Windows Glyph List 4 (WGL4)), and `all' (for all UNICODE characters without restrictions) [[default=wgl4]]]),[
case x"$withval" in
  xascii)  hl_character_set=1
           AC_MSG_RESULT([ASCII only])
           ;;
  xlatin1) hl_character_set=2
           AC_MSG_RESULT([Latin-1 or ISO-8859-1])
           ;;
  xwgl4)   hl_character_set=3
           AC_MSG_RESULT([WGL4])
           ;;
  xall)    hl_character_set=4
           AC_MSG_RESULT([all])
           ;;
  *)       AC_MSG_ERROR([invalid argument --with-character_set=$withval]) ;;
esac],
  hl_character_set=3
  AC_MSG_RESULT([WGL4 (default)])
)

AC_DEFINE_UNQUOTED(CHARACTER_SET,$hl_character_set,[The set of valid characters (1=ASCII, 2=Latin-1, 3=WGL4, 4=all).])


# Which abs() do we use for floats?  Which one for ints?

HL_WHICH_ABS(::std::abs ::abs ::std::fabsf ::fabsf,float)
HL_WHICH_ABS(::std::abs ::abs ::std::fabs  ::fabs ,double)
if test "x$ac_cv_type_long_double" = xyes; then
HL_WHICH_ABS(::std::abs ::abs ::std::fabsl ::fabsl,long double)
else
hl_cv_abs_for_long_double_num=$hl_cv_abs_for_double_num
fi
HL_WHICH_ABS(::std::abs ::abs ::std::abs   ::abs,int)
HL_WHICH_ABS(::std::abs ::abs ::std::labs  ::labs,long int)
dnl if test "x$ac_cv_type_long_long_int" = xyes; then
dnl HL_WHICH_ABS(::std::abs ::abs ::std::llabs ::llabs,long long int)
dnl else
dnl hl_cv_abs_for_long_long_int_num=${hl_cv_abs_for_long_int_num}
dnl fi

AC_MSG_CHECKING([which strategy to use to get abs() right])
if test $hl_cv_abs_for_float_num -gt 2 \
     -o $hl_cv_abs_for_double_num -gt 2 \
     -o $hl_cv_abs_for_long_double_num -gt 2 \
     -o $hl_cv_abs_for_int_num -gt 2 \
     -o $hl_cv_abs_for_long_int_num -gt 2
then
   AC_MSG_RESULT([inline replacement])
   hl_cv_abs_strategy=1
else
 if test $hl_cv_abs_for_float_num -eq 1 \
      -a $hl_cv_abs_for_double_num -eq 1 \
      -a $hl_cv_abs_for_long_double_num -eq 1 \
      -a $hl_cv_abs_for_int_num -eq 1 \
      -a $hl_cv_abs_for_long_int_num -eq 1
 then
    AC_MSG_RESULT([using ::std::abs])
    hl_cv_abs_strategy=2
 else
  if test $hl_cv_abs_for_float_num -eq 2 \
       -a $hl_cv_abs_for_double_num -eq 2 \
       -a $hl_cv_abs_for_long_double_num -eq 2 \
       -a $hl_cv_abs_for_int_num -eq 2 \
       -a $hl_cv_abs_for_long_int_num -eq 2
  then
     AC_MSG_RESULT([using ::abs])
     hl_cv_abs_strategy=3
  else
     AC_MSG_RESULT([using ::std::abs and ::abs])
     hl_cv_abs_strategy=4
  fi
 fi
fi

AC_DEFINE_UNQUOTED(ABS_STRATEGY,$hl_cv_abs_strategy,
   [Strategy to get `abs()' right.])

# Which sin() do we use for floats?
# This result is also used for all other math functions covered by the C++
# standard.

HL_WHICH_MATH_FUN(sin,::std::sin ::sin ::std::sinf ::sinf,float)
HL_WHICH_MATH_FUN(sin,::std::sin ::sin, double)
if test "x$ac_cv_type_long_double" = xyes; then
HL_WHICH_MATH_FUN(sin,::std::sin ::sin ::std::sinl ::sinl,long double)
else
hl_cv_sin_for_long_double_num=$hl_cv_sin_for_double_num
fi

AC_MSG_CHECKING([which strategy to use to get math functions right])
if test $hl_cv_sin_for_float_num       -eq 1 \
     -a $hl_cv_sin_for_double_num      -eq 1 \
     -a $hl_cv_sin_for_long_double_num -eq 1
then
    AC_MSG_RESULT([using ::std::fname])
    hl_cv_math_fun_strategy=2
else
 if test $hl_cv_sin_for_float_num       -eq 2 \
      -a $hl_cv_sin_for_double_num      -eq 2 \
      -a $hl_cv_sin_for_long_double_num -eq 2
 then
    AC_MSG_RESULT([using ::fname])
    hl_cv_math_fun_strategy=3
 else
    AC_MSG_RESULT([inline replacement])
    hl_cv_math_fun_strategy=1
 fi
fi

AC_DEFINE_UNQUOTED(MATH_FUN_STRATEGY,$hl_cv_math_fun_strategy,
   [Strategy to get standard math functions right.])

# Which lgamma() shall we use?

# These tests are not executed, because lgamma() is not used

dnl HL_WHICH_MATH_FUN(lgamma,::std::lgamma ::lgamma ::std::lgammaf ::lgammaf,float)
dnl HL_WHICH_MATH_FUN(lgamma,::std::lgamma ::lgamma,double)
dnl if test "x$ac_cv_type_long_double" = xyes; then
dnl HL_WHICH_MATH_FUN(lgamma,::std::lgamma ::lgamma ::std::lgammal ::lgammal,long double)
dnl else
dnl hl_cv_lgamma_for_long_double_num=$hl_cv_lgamma_for_double_num
dnl fi

# Directories

if test "x$prefix" != xNONE; then   
   AC_DEFINE_UNQUOTED(DATADIR,"${prefix}/share/${PACKAGE}",[The directory which holds the installed data files.])
else
   AC_DEFINE_UNQUOTED(DATADIR,"${ac_default_prefix}/share/${PACKAGE}",[The directory which holds the installed data files.])
fi

AC_DEFINE_UNQUOTED(SRCDATADIR,"`pwd`/share",[The directory which holds the uninstalled data files.])

AC_OUTPUT

