// rxbioclim_omp.h — OpenMP thread-count helper for CRAN compliance.
//
// CRAN requires that packages respect OMP_THREAD_LIMIT and default to
// at most 2 threads during R CMD check.  This helper caps the requested
// thread count accordingly.

#ifndef RXBIOCLIM_OMP_H
#define RXBIOCLIM_OMP_H

#ifdef _OPENMP
#include <omp.h>
#include <cstdlib>

inline int rxbioclim_safe_threads(int requested) {
    const char* limit = std::getenv("OMP_THREAD_LIMIT");
    int max_allowed = limit ? std::atoi(limit) : omp_get_max_threads();
    if (max_allowed < 1) max_allowed = 1;
    if (requested < 1) requested = 1;
    return (requested < max_allowed) ? requested : max_allowed;
}
#else
inline int rxbioclim_safe_threads(int requested) {
    return (requested < 1) ? 1 : requested;
}
#endif

#endif  // RXBIOCLIM_OMP_H
