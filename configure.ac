AC_INIT([rxbioclim], [0.1.0])
AC_CONFIG_AUX_DIR([tools])

# --- GDAL detection (optional) ---
AC_ARG_WITH([gdal-config],
  AS_HELP_STRING([--with-gdal-config=PATH],
    [full path to the gdal-config script]),
  [GDAL_CONFIG="${withval}"],
  [GDAL_CONFIG=""]
)

GDAL_CPPFLAGS=""
GDAL_LIBS=""
HAVE_GDAL=no

# Locate gdal-config
if test -z "${GDAL_CONFIG}"; then
  AC_PATH_PROG([GDAL_CONFIG], [gdal-config], [not-found])
else
  AC_MSG_CHECKING([for gdal-config at ${GDAL_CONFIG}])
  if test -x "${GDAL_CONFIG}"; then
    AC_MSG_RESULT([yes])
  else
    AC_MSG_RESULT([no])
    GDAL_CONFIG="not-found"
  fi
fi

if test "${GDAL_CONFIG}" != "not-found"; then

  # Query GDAL compiler and linker flags
  GDAL_CPPFLAGS=`${GDAL_CONFIG} --cflags 2>/dev/null`
  GDAL_LIBS=`${GDAL_CONFIG} --libs 2>/dev/null`

  # Check GDAL version >= 2.0.1
  GDAL_VERSION=`${GDAL_CONFIG} --version 2>/dev/null`
  AC_MSG_CHECKING([GDAL version >= 2.0.1])
  GDAL_VERSION_MAJOR=`echo "${GDAL_VERSION}" | cut -d. -f1`
  GDAL_VERSION_MINOR=`echo "${GDAL_VERSION}" | cut -d. -f2`
  GDAL_VERSION_PATCH=`echo "${GDAL_VERSION}" | cut -d. -f3`
  GDAL_VERSION_PATCH=${GDAL_VERSION_PATCH:-0}

  if test "${GDAL_VERSION_MAJOR}" -ge 3; then
    GDAL_OK=yes
  elif test "${GDAL_VERSION_MAJOR}" -eq 2 && test "${GDAL_VERSION_MINOR}" -ge 1; then
    GDAL_OK=yes
  elif test "${GDAL_VERSION_MAJOR}" -eq 2 && test "${GDAL_VERSION_MINOR}" -eq 0 && test "${GDAL_VERSION_PATCH}" -ge 1; then
    GDAL_OK=yes
  else
    GDAL_OK=no
  fi

  if test "${GDAL_OK}" = "yes"; then
    AC_MSG_RESULT([yes (${GDAL_VERSION})])

    # Test compilation and linking
    SAVED_CPPFLAGS="${CPPFLAGS}"
    SAVED_LIBS="${LIBS}"
    CPPFLAGS="${GDAL_CPPFLAGS} ${CPPFLAGS}"
    LIBS="${GDAL_LIBS} ${LIBS}"

    AC_LANG_PUSH([C++])
    AC_MSG_CHECKING([whether GDAL links correctly])
    AC_LINK_IFELSE(
      [AC_LANG_PROGRAM(
        [#include <gdal.h>],
        [GDALAllRegister();]
      )],
      [AC_MSG_RESULT([yes])
       HAVE_GDAL=yes],
      [AC_MSG_RESULT([no — GDAL link test failed; building without GDAL])
       GDAL_CPPFLAGS=""
       GDAL_LIBS=""]
    )
    AC_LANG_POP([C++])

    CPPFLAGS="${SAVED_CPPFLAGS}"
    LIBS="${SAVED_LIBS}"

    if test "${HAVE_GDAL}" = "yes"; then
      GDAL_CPPFLAGS="${GDAL_CPPFLAGS} -DHAVE_GDAL"
    fi
  else
    AC_MSG_RESULT([no (${GDAL_VERSION}) — need >= 2.0.1; building without GDAL])
    GDAL_CPPFLAGS=""
    GDAL_LIBS=""
  fi

else
  AC_MSG_WARN([gdal-config not found; building without GDAL (native engine disabled)])
fi

AC_SUBST([GDAL_CPPFLAGS])
AC_SUBST([GDAL_LIBS])

# --- CUDA detection (optional) ---
HAVE_CUDA=no
CUDA_CPPFLAGS=""
CUDA_LIBS=""
NVCC=""
CUDA_OBJECTS=""

AC_PATH_PROG([NVCC], [nvcc], [not-found])
if test "${NVCC}" != "not-found"; then
    CUDA_VERSION=`${NVCC} --version 2>/dev/null | grep "release" | sed 's/.*release //' | sed 's/,.*//'`
    AC_MSG_CHECKING([CUDA toolkit version])
    AC_MSG_RESULT([${CUDA_VERSION}])

    # Try to find CUDA include and lib paths
    CUDA_HOME="${CUDA_HOME:-/usr/local/cuda}"
    if test -d "${CUDA_HOME}/include" && test -d "${CUDA_HOME}/lib64"; then
        HAVE_CUDA=yes
        CUDA_CPPFLAGS="-I${CUDA_HOME}/include -DHAVE_CUDA"
        CUDA_LIBS="-L${CUDA_HOME}/lib64 -lcudart"
    elif test -d "${CUDA_HOME}/include" && test -d "${CUDA_HOME}/lib"; then
        HAVE_CUDA=yes
        CUDA_CPPFLAGS="-I${CUDA_HOME}/include -DHAVE_CUDA"
        CUDA_LIBS="-L${CUDA_HOME}/lib -lcudart"
    else
        AC_MSG_WARN([nvcc found but CUDA_HOME paths not detected; building without CUDA])
        NVCC=""
    fi
else
    NVCC=""
fi

if test "${HAVE_CUDA}" = "yes"; then
    AC_MSG_NOTICE([CUDA support enabled (version ${CUDA_VERSION})])
    CUDA_OBJECTS="bioclim_cuda.o"
else
    AC_MSG_NOTICE([Building without CUDA support (GPU acceleration disabled)])
fi

AC_SUBST([HAVE_CUDA])
AC_SUBST([CUDA_CPPFLAGS])
AC_SUBST([CUDA_LIBS])
AC_SUBST([NVCC])
AC_SUBST([CUDA_OBJECTS])

AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
