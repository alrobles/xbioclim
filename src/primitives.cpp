#include "xbioclim/primitives.hpp"

#include <cmath>
#include <cstddef>
#include <limits>

#ifdef XBIOCLIM_USE_OPENMP
#include <omp.h>
#endif

namespace xbioclim {

Array1D row_mean(const Array2D& A) {
    return xt::mean(A, /*axis=*/1);
}

Array1D row_sum(const Array2D& A) {
    return xt::sum(A, /*axis=*/1);
}

Array1D row_max(const Array2D& A) {
    return xt::amax(A, /*axis=*/1);
}

Array1D row_min(const Array2D& A) {
    return xt::amin(A, /*axis=*/1);
}

Array1D row_std(const Array2D& A) {
    // Population standard deviation: ddof=0, denominator N=12 (xtensor default)
    // Use initializer-list form {1} instead of bare integer 1 for axis argument.
    return xt::stddev(A, {1});
}

Array2D elementwise_diff(const Array2D& A, const Array2D& B) {
    return A - B;
}

IndexArray rolling_quarter_argmax(const Array2D& A) {
    const std::size_t N = A.shape(0);
    IndexArray result = xt::zeros<std::size_t>({N});

#ifdef XBIOCLIM_USE_OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (std::size_t p = 0; p < N; ++p) {
        float best = std::numeric_limits<float>::lowest();
        std::size_t bi = 0;
        bool has_nan = false;
        for (std::size_t i = 0; i < 12; ++i) {
            float a = A(p, i);
            float b = A(p, (i + 1) % 12);
            float c = A(p, (i + 2) % 12);
            if (std::isnan(a) || std::isnan(b) || std::isnan(c)) {
                has_nan = true;
                break;
            }
            float s = a + b + c;
            if (s > best) {
                best = s;
                bi = i;
            }
        }
        result(p) = has_nan ? std::numeric_limits<std::size_t>::max() : bi;
    }
    return result;
}

IndexArray rolling_quarter_argmin(const Array2D& A) {
    const std::size_t N = A.shape(0);
    IndexArray result = xt::zeros<std::size_t>({N});

#ifdef XBIOCLIM_USE_OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (std::size_t p = 0; p < N; ++p) {
        float best = std::numeric_limits<float>::max();
        std::size_t bi = 0;
        bool has_nan = false;
        for (std::size_t i = 0; i < 12; ++i) {
            float a = A(p, i);
            float b = A(p, (i + 1) % 12);
            float c = A(p, (i + 2) % 12);
            if (std::isnan(a) || std::isnan(b) || std::isnan(c)) {
                has_nan = true;
                break;
            }
            float s = a + b + c;
            if (s < best) {
                best = s;
                bi = i;
            }
        }
        result(p) = has_nan ? std::numeric_limits<std::size_t>::max() : bi;
    }
    return result;
}

Array1D quarter_mean(const Array2D& A, const IndexArray& starts) {
    const std::size_t N = A.shape(0);
    Array1D result = xt::zeros<float>({N});
    static constexpr std::size_t SENTINEL = std::numeric_limits<std::size_t>::max();

#ifdef XBIOCLIM_USE_OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (std::size_t p = 0; p < N; ++p) {
        const std::size_t i = starts(p);
        if (i == SENTINEL) {
            result(p) = std::numeric_limits<float>::quiet_NaN();
        } else {
            result(p) = (A(p, i % 12)
                       + A(p, (i + 1) % 12)
                       + A(p, (i + 2) % 12)) / 3.0f;
        }
    }
    return result;
}

} // namespace xbioclim
