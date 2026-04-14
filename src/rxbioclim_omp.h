// rxbioclim_omp.h — safe OpenMP thread-count helper.
//
// CRAN policy requires packages to respect OMP_THREAD_LIMIT (typically set to
// 2 during R CMD check).  Use safe_omp_threads(n) instead of passing n
// directly to omp_set_num_threads() or num_threads() clauses.

#ifndef RXBIOCLIM_OMP_H
#define RXBIOCLIM_OMP_H

#ifdef _OPENMP
#include <omp.h>
#include <cerrno>
#include <climits>
#include <cstdlib>
#include <algorithm>

// Returns the number of threads to use, capped by OMP_THREAD_LIMIT (if set)
// and by omp_get_thread_limit() (the hard upper bound set at program start
// via OMP_THREAD_LIMIT or the runtime), but always at least 1.
//
// Note: omp_get_max_threads() reflects the current nthreads-var (which changes
// with omp_set_num_threads()), so it is not used here.  omp_get_thread_limit()
// is the immutable ceiling and is the correct upper bound to enforce.
inline int safe_omp_threads(int requested) {
    int limit = (requested < 1) ? 1 : requested;

    // Respect OMP_THREAD_LIMIT — set to 2 by CRAN during R CMD check.
    // Use strtol for safe parsing: distinguishes non-numeric input from "0",
    // and avoids undefined behavior on overflow.
    const char* env = std::getenv("OMP_THREAD_LIMIT");
    if (env != nullptr) {
        char* end = nullptr;
        errno = 0;
        long env_val = std::strtol(env, &end, 10);
        if (end != env && errno == 0 && env_val > 0 && env_val <= INT_MAX) {
            limit = std::min(limit, static_cast<int>(env_val));
        }
    }

    // Never exceed the hard thread-limit set by the OpenMP runtime.
    // omp_get_thread_limit() returns the value of OMP_THREAD_LIMIT as
    // interpreted by the runtime (INT_MAX when unset); it is immutable during
    // a program run, unlike omp_get_max_threads() which tracks nthreads-var.
    int omp_limit = omp_get_thread_limit();
    if (omp_limit > 0) limit = std::min(limit, omp_limit);

    return (limit < 1) ? 1 : limit;
}

#else
// Without OpenMP, always return 1.
inline int safe_omp_threads(int /*requested*/) { return 1; }
#endif

#endif // RXBIOCLIM_OMP_H
