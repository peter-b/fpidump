#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])
AC_INIT([FPI Extraction Tool], [1.0], [p.brett@surrey.ac.uk], [fpidump])

AC_CONFIG_SRCDIR([fpidump.c])
AC_CONFIG_AUX_DIR([build-tools])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS([config.h])

AM_INIT_AUTOMAKE([foreign])

# We generate some files differently if a git repo is present
AX_GIT_REPO
AM_CONDITIONAL([HAVE_GIT_REPO], test "X$HAVE_GIT_REPO" = "Xyes")

# Checks for programs.
AC_PROG_CC
AC_PROG_CC_C99
AM_PROG_CC_C_O
PKG_PROG_PKG_CONFIG

# Checks for libraries.
PKG_CHECK_MODULES([GLIB], [glib-2.0 >= 2.6.0], ,
  AC_MSG_ERROR([GLIB 2.6.0 or later is required.]))
AC_CHECK_LIB([tiff], [TIFFOpen])

# Checks for header files.
AC_CHECK_HEADERS([stddef.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_HEADER_STDBOOL
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([strerror strchr])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
