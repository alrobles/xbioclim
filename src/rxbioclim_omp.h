// rxbioclim_omp.h — safe OpenMP thread-count helper.
//
// CRAN policy requires packages to respect OMP_THREAD_LIMIT (typically set to
// 2 during R CMD check).  Use safe_omp_threads(n) instead of passing n
// directly to omp_set_num_threads() or num_threads() clauses.

#ifndef RXBIOCLIM_OMP_H
#define RXBIOCLIM_OMP_H

#ifdef _OPENMP
#include <omp.h>
#include <cstdlib>
#include <algorithm>

// Returns the number of threads to use, capped by OMP_THREAD_LIMIT (if set)
// and by the OpenMP runtime maximum, but always at least 1.
inline int safe_omp_threads(int requested) {
    int limit = (requested < 1) ? 1 : requested;

    // Respect OMP_THREAD_LIMIT — set to 2 by CRAN during R CMD check.
    const char* env = std::getenv("OMP_THREAD_LIMIT");
    if (env != nullptr) {
        int env_limit = std::atoi(env);
        if (env_limit > 0) limit = std::min(limit, env_limit);
    }

    // Never exceed what the OpenMP runtime currently allows.
    limit = std::min(limit, omp_get_max_threads());

    return (limit < 1) ? 1 : limit;
}

#else
// Without OpenMP, always return 1.
inline int safe_omp_threads(int /*requested*/) { return 1; }
#endif

#endif // RXBIOCLIM_OMP_H
