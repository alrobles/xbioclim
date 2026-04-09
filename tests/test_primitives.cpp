#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "xbioclim/primitives.hpp"

using Catch::Matchers::WithinAbs;
using namespace xbioclim;

// Helper: build a [N, 12] array where row p has values [1..12] scaled by p+1
static Array2D make_ramp(std::size_t N = 1) {
    Array2D A = xt::zeros<float>({N, std::size_t(12)});
    for (std::size_t p = 0; p < N; ++p)
        for (std::size_t m = 0; m < 12; ++m)
            A(p, m) = static_cast<float>(m + 1);
    return A;
}

TEST_CASE("row_mean returns correct mean for ramp data", "[primitives]") {
    Array2D A = make_ramp(3);
    Array1D result = row_mean(A);
    REQUIRE(result.size() == 3);
    // mean(1..12) = 6.5
    for (std::size_t p = 0; p < 3; ++p)
        CHECK_THAT(result(p), WithinAbs(6.5f, 1e-4f));
}

TEST_CASE("row_sum returns correct sum for ramp data", "[primitives]") {
    Array2D A = make_ramp(2);
    Array1D result = row_sum(A);
    REQUIRE(result.size() == 2);
    // sum(1..12) = 78
    for (std::size_t p = 0; p < 2; ++p)
        CHECK_THAT(result(p), WithinAbs(78.0f, 1e-4f));
}

TEST_CASE("row_max returns correct maximum", "[primitives]") {
    Array2D A = make_ramp(2);
    Array1D result = row_max(A);
    REQUIRE(result.size() == 2);
    for (std::size_t p = 0; p < 2; ++p)
        CHECK_THAT(result(p), WithinAbs(12.0f, 1e-4f));
}

TEST_CASE("row_min returns correct minimum", "[primitives]") {
    Array2D A = make_ramp(2);
    Array1D result = row_min(A);
    REQUIRE(result.size() == 2);
    for (std::size_t p = 0; p < 2; ++p)
        CHECK_THAT(result(p), WithinAbs(1.0f, 1e-4f));
}

TEST_CASE("row_std returns population SD for ramp data", "[primitives]") {
    Array2D A = make_ramp(1);
    Array1D result = row_std(A);
    // Population SD of [1..12]: sqrt(143/12) ≈ 3.4521
    const float expected = std::sqrt(143.0f / 12.0f);
    CHECK_THAT(result(0), WithinAbs(expected, 1e-3f));
}

TEST_CASE("elementwise_diff produces A-B", "[primitives]") {
    Array2D A = make_ramp(1);
    // B has all entries = 1
    Array2D B = xt::ones<float>({std::size_t(1), std::size_t(12)});
    Array2D diff = elementwise_diff(A, B);
    for (std::size_t m = 0; m < 12; ++m)
        CHECK_THAT(diff(0, m), WithinAbs(static_cast<float>(m), 1e-4f));
}

TEST_CASE("rolling_quarter_argmax returns correct index for ramp pr",
          "[primitives]") {
    // pr = [1..12] — wettest quarter starts at index 9 (months 10,11,12)
    Array2D A = make_ramp(1);
    IndexArray idx = rolling_quarter_argmax(A);
    REQUIRE(idx(0) == 9u);
}

TEST_CASE("rolling_quarter_argmin returns correct index for ramp pr",
          "[primitives]") {
    // pr = [1..12] — driest quarter starts at index 0 (months 1,2,3)
    Array2D A = make_ramp(1);
    IndexArray idx = rolling_quarter_argmin(A);
    REQUIRE(idx(0) == 0u);
}

TEST_CASE("quarter_mean gives correct mean for wettest quarter",
          "[primitives]") {
    // With ramp data, wettest start index = 9 → months 10,11,12 → mean = 11
    Array2D A = make_ramp(1);
    IndexArray starts = xt::zeros<std::size_t>({std::size_t(1)});
    starts(0) = 9u;
    Array1D result = quarter_mean(A, starts);
    CHECK_THAT(result(0), WithinAbs(11.0f, 1e-3f));
}

TEST_CASE("quarter_mean handles circular wrap-around", "[primitives]") {
    // Start index 11 → months 12, 1, 2 → values 12, 1, 2 → mean = 5
    Array2D A = make_ramp(1);
    IndexArray starts = xt::zeros<std::size_t>({std::size_t(1)});
    starts(0) = 11u;
    Array1D result = quarter_mean(A, starts);
    CHECK_THAT(result(0), WithinAbs(5.0f, 1e-3f));
}
