AC_INIT([lalsimulation],[0.0.0.1],[lal-discuss@gravity.phys.uwm.edu])
AM_CONFIG_HEADER([include/config.h src/LALSimulationConfig.h])
AC_CONFIG_SRCDIR([src/LALSimulation.c])
AC_CONFIG_AUX_DIR([gnuscripts])
AC_CONFIG_MACRO_DIR([gnuscripts])
AC_CONFIG_FILES([ \
  lalsimulation.spec \
  lalsimulation.pc \
  Makefile \
  gnuscripts/Makefile \
  debian/Makefile \
  include/Makefile \
  include/lal/Makefile \
  src/Makefile \
  swig/Makefile \
  swig/octave/Makefile \
  swig/python/Makefile \
  test/Makefile
])
AM_INIT_AUTOMAKE([foreign])

AC_CANONICAL_HOST

LALSUITE_ARG_VAR

# set library api info
# 1. start with version information of "0:0:0" for each libtool library
# 2. update the version information only immediately before a public
#    release of your software. More frequent updates are unnecessary,
#    and only guarantee that the current interface number gets larger
#    faster.
# 3. if the library source code has changed at all since the last
#    update, then increment revision ("c:r:a" becomes "c:r+1:a").
# 4. if any interfaces have been added, removed, or changed since the
#    last update, increment current, and set revision to 0.
# 5. if any interfaces have been added since the last public release,
#    then increment age.
# 6. if any interfaces have been removed since the last public release,
#    then set age to 0.
AC_SUBST([LIBCURRENT],[0])
AC_SUBST([LIBREVISION],[0])
AC_SUBST([LIBAGE],[0])
AC_SUBST([LIBVERSION],[${LIBCURRENT}:${LIBREVISION}:${LIBAGE}])

# default install path
AC_PREFIX_DEFAULT([/opt/lscsoft/lalsimulation])

# nightly build
LALSUITE_ENABLE_NIGHTLY

# debug
LALSUITE_ENABLE_DEBUG

# put version/configure info in config header
lalsimulation_version_major=`echo "$VERSION" | cut -d. -f1`
lalsimulation_version_minor=`echo "$VERSION" | cut -d. -f2`
lalsimulation_version_micro=`echo "$VERSION" | cut -d. -f3`
lalsimulation_version_devel=`echo "$VERSION" | cut -d. -f4-`
test -z "$lalsimulation_version_micro" && lalsimulation_version_micro=0
test -z "$lalsimulation_version_devel" && lalsimulation_version_devel=0
lalsimulation_configure_date=`date +"%Y-%m-%dT%H:%M:%S%z"`
AC_DEFINE_UNQUOTED([LALSIMULATION_VERSION],["$VERSION"],[LALSimulation Version])
AC_DEFINE_UNQUOTED([LALSIMULATION_VERSION_MAJOR],[$lalsimulation_version_major],
                   [LALSimulation Version Major Number])
AC_DEFINE_UNQUOTED([LALSIMULATION_VERSION_MINOR],[$lalsimulation_version_minor],
                   [LALSimulation Version Minor Number])
AC_DEFINE_UNQUOTED([LALSIMULATION_VERSION_MICRO],[$lalsimulation_version_micro],
                   [LALSimulation Version Micro Number])
AC_DEFINE_UNQUOTED([LALSIMULATION_VERSION_DEVEL],[$lalsimulation_version_devel],
                   [LALSimulation Version Devel Number])
AC_DEFINE_UNQUOTED([LALSIMULATION_CONFIGURE_ARGS],["$ac_configure_args"],
                   [LALSimulation Configure Arguments])
AC_DEFINE_UNQUOTED([LALSIMULATION_CONFIGURE_DATE],["$lalsimulation_configure_date"],
                   [LALSimulation Configure Date])

# enable strict gcc flags
LALSUITE_ENABLE_GCC_FLAGS

# check for c compiler
m4_pattern_allow([AC_PROG_CC_C99])
m4_ifdef([AC_PROG_CC_C99],[AC_PROG_CC_C99],[LALSUITE_AC_PROG_CC_C99])

# use silent build rules if available
m4_ifdef([AM_SILENT_RULES],[AM_SILENT_RULES([yes])],)

# SWIG wrapping requires a c++ compiler
AC_PROG_CXX

# checks for programs
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_AWK
LT_INIT
AM_PATH_PYTHON([2.4],,[AC_MSG_ERROR([Python-2.4, or higher, is required])])
AC_CHECK_PROGS([M4],[gm4 m4],[m4])

# check for supported mac os x version
LALSUITE_OSX_VERSION_CHECK

# link tests using libtool
LALSUITE_USE_LIBTOOL

# check for header files
AC_HEADER_STDC

# check for gethostname in unistd.h
AC_MSG_CHECKING([for gethostname prototype in unistd.h])
AC_EGREP_HEADER([gethostname],[unistd.h],[AC_MSG_RESULT([yes])]
AC_DEFINE([HAVE_GETHOSTNAME_PROTOTYPE],[1],
          [Define if gethostname prototype is in unistd.h]),
          [AC_MSG_RESULT]([no]))

# check for required libraries
AC_CHECK_LIB([m],[main],,[AC_MSG_ERROR([could not find the math library])])

# check for gsl
PKG_CHECK_MODULES([GSL],[gsl],[true],[false])
CPPFLAGS="$CPPFLAGS $GSL_CFLAGS"
LIBS="$LIBS $GSL_LIBS"
AC_CHECK_LIB([gslcblas],[main],,[AC_MSG_ERROR([could not find the gsl library])])
AC_CHECK_LIB([gsl],[gsl_strerror],,[AC_MSG_ERROR([could not find the gsl library])])

# check for gsl headers
AC_CHECK_HEADERS([gsl/gsl_errno.h],,[AC_MSG_ERROR([could not find the gsl/gsl_errno.h header])])

# check for lal libraries and headers
LALSUITE_CHECK_LIB([LAL],[6.5.0],[LALVersion],[lal/LALStdio.h])
LALSUITE_CHECK_LIB([LALSupport],[6.5.0],[LALOpenDataFile],[lal/FileIO.h])

# check for SWIG wrapping
LALSUITE_ENABLE_SWIG

# add gcc specific flags
if test "$GCC" = yes; then
  CFLAGS="$CFLAGS $lal_gcc_flags"

  # add mac os x specific flags
  if test "x$MACOSX_VERSION" != "x"; then
    CFLAGS="$CFLAGS -mmacosx-version-min=10.4"
  fi
fi

# output
LAL_HEADER_LINK=${ac_aux_dir}/header-link.mk
AC_SUBST_FILE([LAL_HEADER_LINK])
AC_PROG_MAKE_SET
AC_OUTPUT
echo "
================================================================

        LALSimulation has now been successfully configured:

                SWIG wrappings for Octave are $SWIG_OCTAVE_ENABLE_VAL
                SWIG wrappings for Python are $SWIG_PYTHON_ENABLE_VAL

        and will be installed under the directory:

          ${prefix}

        Now run 'make' to build LALSimulation
        and run 'make install' to install LALSimulation

================================================================="