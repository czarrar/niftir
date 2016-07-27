AC_INIT([niftir],[0.0.1])
AC_CONFIG_HEADER(src/config.h)
# AM_INIT_AUTOMAKE

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

# Checks for header files.
AC_HEADER_STDC

# Checks for libraries.

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE

# Checks for library functions.
AC_FUNC_MALLOC

dnl **********************************************
dnl * get infos about the installed pkg-config   *
dnl **********************************************
AC_PATH_PROG(PKG_CONFIG,pkg-config)

CHECK_ZLIB()
AC_SUBST(PACKAGE_CFLAGS)
AC_SUBST(PACKAGE_LIBS)
AC_SUBST(LIBS)
AC_SUBST(INCLUDE)
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
