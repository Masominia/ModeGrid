// A structure uniform in x must reproduce the analytic slab index.
#include <cmath>

#include "check.hpp"
#include "modegrid/mode.hpp"

// Exact index of the fundamental mode of a symmetric slab by bisection on the even-mode equation.
static double slab_exact(double lambda, double t, double nc, double ncl, bool tm) {
    const double k = 2.0 * std::acos(-1.0) / lambda;
    const double weight = tm ? (nc * nc) / (ncl * ncl) : 1.0;
    double lo = ncl + 1e-12, hi = nc - 1e-12;
    for (int i = 0; i < 200; ++i) {
        const double n = 0.5 * (lo + hi);
        const double kappa = k * std::sqrt(nc * nc - n * n), gamma = k * std::sqrt(n * n - ncl * ncl);
        // kappa t / 2 - atan(weight gamma / kappa) falls with n
        if (0.5 * kappa * t - std::atan(weight * gamma / kappa) > 0.0) lo = n;
        else hi = n;
    }
    return 0.5 * (lo + hi);
}

void test_slab(Suite& suite) {
    using namespace mg;
    const double lambda = 1.55, k = 2.0 * std::acos(-1.0) / lambda, width = 8.0;
    struct Case {
        double t, nc;
        bool tm;
        double tol;
    };
    for (const Case& c : {Case{0.22, 3.476, false, 1.5e-3}, Case{0.22, 3.476, true, 3e-3}, Case{0.40, 1.996, false, 1e-3}, Case{0.40, 1.996, true, 1e-3}}) {
        Grid g = make_grid(40, static_cast<int>(std::lround((2.4 + c.t) / 0.01)), 0.1, 0.01, 0.0, -1.2, 1.444 * 1.444);
        paint(g, -100.0, 100.0, 0.0, c.t, c.nc * c.nc);
        const Field field = c.tm ? Field::Tm : Field::Te;
        const Mode m = solve_mode(g, lambda, field, Symmetry::Even);
        const double exact = slab_exact(lambda, c.t, c.nc, 1.444, c.tm);
        // between the mirror and the wall the field is half a cosine
        const double expected = std::sqrt(exact * exact - std::pow(std::acos(-1.0) / width / k, 2));
        CHECK_NEAR(m.n_eff, expected, c.tol);
        // not confined sideways: below the slab index found on the same grid
        CHECK(!m.guided);
        CHECK_NEAR(cladding_index(g, lambda, field, true), exact, c.tol);
    }
}
