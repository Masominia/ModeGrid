#include <cmath>

#include "check.hpp"
#include "modegrid/eigen.hpp"
#include "modegrid/wave_operator.hpp"

void test_wave_operator(Suite& suite) {
    using namespace mg;
    const double pi = std::acos(-1.0), lambda = 1.55, k = 2.0 * pi / lambda;
    // numbering follows the shorter side
    const Grid tall = make_grid(4, 9, 0.1, 0.1, 0, 0, 1.0), wide = make_grid(9, 4, 0.1, 0.1, 0, 0, 1.0);
    CHECK(Ordering(tall).half_width() == 4 && Ordering(tall).index(1, 2) == 9);
    CHECK(Ordering(wide).half_width() == 4 && Ordering(wide).index(2, 1) == 9);

    // uniform medium in a metal box: the lowest mode is a product of half sines
    const int nx = 24, ny = 16;
    const double dx = 0.05, dy = 0.04, n = 2.0;
    const Grid box = make_grid(nx, ny, dx, dy, 0, 0, n * n);
    const double lx = std::pow(2.0 / dx * std::sin(pi / (2.0 * nx)), 2), ly = std::pow(2.0 / dy * std::sin(pi / (2.0 * ny)), 2);
    for (Field f : {Field::Te, Field::Tm, Field::Scalar}) {
        const EigenPair p = largest_eigenpair(assemble(box, lambda, f, Symmetry::None), 1.001 * k * k * n * n);
        CHECK_NEAR(p.value, k * k * n * n - lx - ly, 1e-8);
    }
    // an even mirror doubles the box: half a cosine across twice the width
    const double lx_even = std::pow(2.0 / dx * std::sin(pi / (4.0 * nx)), 2);
    CHECK_NEAR(largest_eigenpair(assemble(box, lambda, Field::Te, Symmetry::Even), 1.001 * k * k * n * n).value, k * k * n * n - lx_even - ly, 1e-8);
    // an odd mirror is a wall
    CHECK_NEAR(largest_eigenpair(assemble(box, lambda, Field::Te, Symmetry::Odd), 1.001 * k * k * n * n).value, k * k * n * n - lx - ly, 1e-8);

    // interface coefficients: across a vertical interface the te operator weights the neighbour by 2 eps_q / (eps_p + eps_q)
    Grid step = make_grid(4, 6, 0.1, 0.1, 0, 0, 1.0);
    for (int iy = 0; iy < 6; ++iy) step.at(2, iy) = step.at(3, iy) = 4.0;
    const Ordering order(step);
    const Banded te = assemble(step, lambda, Field::Te, Symmetry::None), tm = assemble(step, lambda, Field::Tm, Symmetry::None);
    CHECK_NEAR(te(order.index(1, 3), order.index(2, 3)), 2.0 * 4.0 / 5.0 / 0.01, 1e-9);
    CHECK_NEAR(te(order.index(2, 3), order.index(1, 3)), 2.0 * 1.0 / 5.0 / 0.01, 1e-9);
    CHECK_NEAR(tm(order.index(1, 3), order.index(2, 3)), 100.0, 1e-9);
    CHECK_NEAR(te(order.index(1, 3), order.index(1, 4)), 100.0, 1e-9);
    // away from the interface and the walls the diagonal is k^2 eps - 2/dx^2 - 2/dy^2
    Grid plain = make_grid(6, 6, 0.1, 0.1, 0, 0, 1.0);
    const Banded sc = assemble(plain, lambda, Field::Scalar, Symmetry::None);
    CHECK_NEAR(sc(Ordering(plain).index(2, 2), Ordering(plain).index(2, 2)), k * k - 400.0, 1e-9);
}
