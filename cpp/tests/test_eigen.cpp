#include <cmath>
#include <vector>

#include "check.hpp"
#include "modegrid/eigen.hpp"

void test_eigen(Suite& suite) {
    using namespace mg;
    // second difference with fixed ends: eigenvalues -4 sin^2(j pi / (2 (n + 1))), vectors sin(i j pi / (n + 1))
    const std::size_t n = 40;
    const double pi = std::acos(-1.0);
    Banded a(n, 1);
    for (std::size_t i = 0; i < n; ++i) {
        a(i, i) = -2.0;
        if (i + 1 < n) a(i, i + 1) = 1.0, a(i + 1, i) = 1.0;
    }
    const EigenPair pair = largest_eigenpair(a, 0.01);
    const double exact = -4.0 * std::pow(std::sin(pi / (2.0 * (n + 1))), 2);
    CHECK(pair.converged);
    CHECK_NEAR(pair.value, exact, 1e-12);
    double norm = 0.0, worst = 0.0;
    for (double v : pair.vector) norm += v * v;
    CHECK_NEAR(norm, 1.0, 1e-12);
    const double scale = std::sqrt(2.0 / (n + 1));
    for (std::size_t i = 0; i < n; ++i) worst = std::fmax(worst, std::fabs(pair.vector[i] - scale * std::sin((i + 1) * pi / (n + 1))));
    CHECK(worst < 1e-8);
    CHECK(pair.vector[n / 2] > 0.0);

    // the residual A x - value x vanishes
    const std::vector<double> ax = a.multiply(pair.vector);
    double residual = 0.0;
    for (std::size_t i = 0; i < n; ++i) residual = std::fmax(residual, std::fabs(ax[i] - pair.value * pair.vector[i]));
    CHECK(residual < 1e-10);

    // a non-symmetric matrix with real eigenvalues: diagonal scaling of the same operator
    Banded skew(n, 1);
    for (std::size_t i = 0; i < n; ++i) {
        skew(i, i) = -2.0;
        if (i + 1 < n) skew(i, i + 1) = 2.0, skew(i + 1, i) = 0.5;
    }
    CHECK_NEAR(largest_eigenpair(skew, 0.01).value, exact, 1e-12);

    // too few iterations: reported as not converged
    CHECK(!largest_eigenpair(a, 5.0, 1e-15, 2).converged);
}
