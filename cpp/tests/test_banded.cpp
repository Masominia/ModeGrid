#include <vector>

#include "check.hpp"
#include "modegrid/banded.hpp"

void test_banded(Suite& suite) {
    using namespace mg;
    // a diagonally dominant pentadiagonal matrix
    const std::size_t n = 12;
    Banded a(n, 2);
    for (std::size_t i = 0; i < n; ++i) {
        a(i, i) = 6.0 + 0.1 * static_cast<double>(i);
        if (i + 1 < n) a(i, i + 1) = -1.0, a(i + 1, i) = -1.5;
        if (i + 2 < n) a(i, i + 2) = 0.5, a(i + 2, i) = 0.25;
    }
    CHECK(a.in_band(3, 5));
    CHECK(!a.in_band(3, 6));
    std::vector<double> x(n);
    for (std::size_t i = 0; i < n; ++i) x[i] = 1.0 + 0.3 * static_cast<double>(i * i % 7);
    std::vector<double> b = a.multiply(x);
    // the first row by hand
    CHECK_NEAR(b[0], 6.0 * x[0] - x[1] + 0.5 * x[2], 1e-12);

    Banded lu = a;
    CHECK_THROWS(lu.solve(b));  // not factorised yet
    lu.factorise();
    lu.solve(b);
    double worst = 0.0;
    for (std::size_t i = 0; i < n; ++i) worst = std::fmax(worst, std::fabs(b[i] - x[i]));
    CHECK(worst < 1e-12);

    // shifting the diagonal
    Banded shifted = a;
    shifted.add_to_diagonal(2.0);
    CHECK_NEAR(shifted(4, 4), a(4, 4) + 2.0, 1e-12);
    CHECK_NEAR(shifted(4, 5), a(4, 5), 1e-12);

    // a band wider than the matrix is clipped
    Banded small(3, 10);
    CHECK(small.half_width() == 2);
    // singular matrix
    Banded zero(3, 1);
    CHECK_THROWS(zero.factorise());
    std::vector<double> wrong(5, 1.0);
    CHECK_THROWS(a.multiply(wrong));
}
