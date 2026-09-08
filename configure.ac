AC_INIT([xbioclim], [1.0.0])
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
  if test -z "${GDAL_VERSION_PATCH}"; then
    GDAL_VERSION_PATCH=0
  fi

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
# Accept --with-cuda / --without-cuda so users and CRAN reviewers can opt out
# of CUDA even when nvcc is present on the build machine.
AC_ARG_WITH([cuda],
  AS_HELP_STRING([--without-cuda],
    [disable CUDA GPU support even if nvcc is found (default: auto-detect)]),
  [with_cuda="${withval}"],
  [with_cuda=auto]
)

HAVE_CUDA=no
CUDA_CPPFLAGS=""
CUDA_LIBS=""
NVCC=""
CUDA_OBJECTS=""
CUDA_ARCH_FLAGS=""

if test "${with_cuda}" = "no"; then
  AC_MSG_NOTICE([CUDA disabled via --without-cuda])
else
  AC_PATH_PROG([NVCC], [nvcc], [not-found])
  if test "${NVCC}" != "not-found"; then
    CUDA_VERSION=`${NVCC} --version 2>/dev/null | grep "release" | sed 's/.*release //' | sed 's/,.*//'`
    AC_MSG_CHECKING([CUDA toolkit version])
    AC_MSG_RESULT([${CUDA_VERSION}])

    # Determine CUDA_HOME — derive from nvcc path if not set explicitly.
    if test -z "${CUDA_HOME}"; then
      # Try to infer from nvcc location: .../cuda/X.Y/bin/nvcc → .../cuda/X.Y
      CUDA_HOME=`dirname "${NVCC}"`   # .../bin
      CUDA_HOME=`dirname "${CUDA_HOME}"`  # .../cuda/X.Y or /usr/local/cuda
    fi

    # Locate the CUDA library directory
    CUDA_LIB_DIR=""
    if test -d "${CUDA_HOME}/lib64"; then
      CUDA_LIB_DIR="${CUDA_HOME}/lib64"
    elif test -d "${CUDA_HOME}/lib"; then
      CUDA_LIB_DIR="${CUDA_HOME}/lib"
    fi

    if test -d "${CUDA_HOME}/include" && test -n "${CUDA_LIB_DIR}"; then
      # Perform an actual compile test — verify nvcc produces object code
      AC_MSG_CHECKING([whether nvcc can compile a trivial CUDA program])
      if test -z "${TMPDIR}"; then
        CUDA_TMPDIR=/tmp
      else
        CUDA_TMPDIR="${TMPDIR}"
      fi
      CUDA_TEST_SRC="${CUDA_TMPDIR}/conftest_cuda_$$.cu"
      CUDA_TEST_OBJ="${CUDA_TMPDIR}/conftest_cuda_$$.o"
      echo '__global__ void _xbioclim_check(void) {}' > "${CUDA_TEST_SRC}"
      echo 'int main(void) { return 0; }'             >> "${CUDA_TEST_SRC}"
      if "${NVCC}" -x cu -I"${CUDA_HOME}/include" \
           -Xcompiler -fPIC \
           -c "${CUDA_TEST_SRC}" -o "${CUDA_TEST_OBJ}" \
           >/dev/null 2>&1; then
        AC_MSG_RESULT([yes])
        HAVE_CUDA=yes
        CUDA_CPPFLAGS="-I${CUDA_HOME}/include -DHAVE_CUDA"
        CUDA_LIBS="-L${CUDA_LIB_DIR} -lcudart"
        # Architecture flags — adapt to the CUDA toolkit version because
        # newer toolkits drop support for older GPU architectures:
        #   CUDA 11.x : sm_60 – sm_86
        #   CUDA 12.x : sm_60 – sm_90  (Pascal through Hopper)
        #   CUDA 13.0+: sm_75 – sm_90  (Maxwell/Pascal/Volta dropped)
        CUDA_MAJOR=`echo "${CUDA_VERSION}" | cut -d. -f1`
        AC_MSG_CHECKING([CUDA architecture flags for toolkit ${CUDA_MAJOR}.x])
        if test "${CUDA_MAJOR}" -ge 13 2>/dev/null; then
          # CUDA 13+: Turing (sm_75) is the oldest supported architecture.
          CUDA_ARCH_FLAGS="-gencode arch=compute_75,code=sm_75 -gencode arch=compute_80,code=sm_80 -gencode arch=compute_86,code=sm_86 -gencode arch=compute_89,code=sm_89 -gencode arch=compute_90,code=sm_90"
        elif test "${CUDA_MAJOR}" -ge 12 2>/dev/null; then
          # CUDA 12.x: Pascal (sm_60) through Hopper (sm_90).
          CUDA_ARCH_FLAGS="-gencode arch=compute_60,code=sm_60 -gencode arch=compute_70,code=sm_70 -gencode arch=compute_75,code=sm_75 -gencode arch=compute_80,code=sm_80 -gencode arch=compute_86,code=sm_86 -gencode arch=compute_89,code=sm_89 -gencode arch=compute_90,code=sm_90"
        else
          # CUDA 11.x and older: Pascal through Ampere.
          CUDA_ARCH_FLAGS="-gencode arch=compute_60,code=sm_60 -gencode arch=compute_70,code=sm_70 -gencode arch=compute_75,code=sm_75 -gencode arch=compute_80,code=sm_80 -gencode arch=compute_86,code=sm_86"
        fi
        AC_MSG_RESULT([${CUDA_ARCH_FLAGS}])
      else
        AC_MSG_RESULT([no — nvcc compile test failed; building without CUDA])
      fi
      rm -f "${CUDA_TEST_SRC}" "${CUDA_TEST_OBJ}"
    else
      AC_MSG_WARN([nvcc found but CUDA_HOME (${CUDA_HOME}) headers/libs not detected; building without CUDA])
    fi
  else
    NVCC=""
  fi

  if test "${with_cuda}" = "yes" && test "${HAVE_CUDA}" = "no"; then
    AC_MSG_ERROR([--with-cuda was requested but CUDA could not be configured])
  fi
fi

if test "${HAVE_CUDA}" = "yes"; then
  AC_MSG_NOTICE([CUDA support enabled (version ${CUDA_VERSION})])
  CUDA_OBJECTS="bioclim_cuda.o"
  # When CUDA is present, compile bioclim_cuda.cpp with nvcc instead of g++.
  CUDA_NVCC_RULE="bioclim_cuda.o: bioclim_cuda.cpp
	\$(NVCC) -x cu \$(CUDA_CPPFLAGS) \$(CUDA_ARCH_FLAGS) -Xcompiler -fPIC -c \$< -o \$@"
else
  AC_MSG_NOTICE([Building without CUDA support (GPU acceleration disabled)])
  CUDA_NVCC_RULE=""
fi

# Enumerate all .o files explicitly (no GNU Make extensions needed).
# bioclim_cuda.o is always included: when CUDA is absent, g++ compiles it
# to an empty object because all code is guarded by #ifdef HAVE_CUDA.
ALL_OBJECTS="bioclim.o bioclim_cuda.o bioclim_model.o bioclim_xt.o BioclimEngine.o cuda_utils.o gdal_io.o gdal_io_rcpp.o gdal_mask.o gdal_mask_rcpp.o rcpp_bioclim_mod.o RcppExports.o"

AC_SUBST([HAVE_CUDA])
AC_SUBST([CUDA_CPPFLAGS])
AC_SUBST([CUDA_LIBS])
AC_SUBST([NVCC])
AC_SUBST([CUDA_OBJECTS])
AC_SUBST([CUDA_ARCH_FLAGS])
AC_SUBST([CUDA_NVCC_RULE])
AC_SUBST([ALL_OBJECTS])

AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
