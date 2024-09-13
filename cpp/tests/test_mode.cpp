#include <cmath>

#include "check.hpp"
#include "modegrid/mode.hpp"

static mg::Grid strip(double dx, double dy, bool half) {
    mg::Grid g = mg::make_grid(half ? 48 : 96, 91, dx, dy, half ? 0.0 : -1.2, -0.8, 1.444 * 1.444);
    mg::paint(g, -0.225, 0.225, 0.0, 0.22, 3.476 * 3.476);
    return g;
}

void test_mode(Suite& suite) {
    using namespace mg;
    const Grid half = strip(0.025, 0.02, true), full = strip(0.025, 0.02, false);
    const Mode te = solve_mode(half, 1.55, Field::Te, Symmetry::Even);
    CHECK(te.converged && te.guided);
    CHECK(te.n_eff > 2.3 && te.n_eff < 2.5);
    const Mode tm = solve_mode(half, 1.55, Field::Tm, Symmetry::Even), scalar = solve_mode(half, 1.55, Field::Scalar, Symmetry::Even);
    CHECK(tm.n_eff < te.n_eff && te.n_eff < scalar.n_eff);
    CHECK(tm.guided);

    // the whole structure gives the same mode as half of it with a mirror
    const Mode whole = solve_mode(full, 1.55, Field::Te, Symmetry::None);
    CHECK_NEAR(whole.n_eff, te.n_eff, 1e-9);
    const ModeShape a = shape_of(half, te, true), b = shape_of(full, whole, false);
    CHECK_NEAR(a.confinement, b.confinement, 1e-7);
    CHECK_NEAR(a.area, b.area, 1e-7);
    CHECK_NEAR(a.rms_x, b.rms_x, 1e-7);
    CHECK_NEAR(a.rms_y, b.rms_y, 1e-7);
    CHECK(a.confinement > 0.5 && a.confinement < 0.9);
    CHECK(a.area > 0.08 && a.area < 0.3);
    CHECK(a.centre_y > 0.05 && a.centre_y < 0.17);
    CHECK(std::fabs(b.peak_x) < 0.03);

    // the antisymmetric mode of a 450 nm strip is not confined; the solver finds a box mode instead
    const Mode odd = solve_mode(half, 1.55, Field::Te, Symmetry::Odd);
    CHECK(!odd.guided && odd.n_eff < 1.444);
    // a wider strip guides it
    Grid wide = make_grid(60, 91, 0.025, 0.02, 0.0, -0.8, 1.444 * 1.444);
    paint(wide, -0.4, 0.4, 0.0, 0.22, 3.476 * 3.476);
    const Mode odd_wide = solve_mode(wide, 1.55, Field::Te, Symmetry::Odd);
    CHECK(odd_wide.guided && odd_wide.n_eff < solve_mode(wide, 1.55, Field::Te, Symmetry::Even).n_eff);

    // longer wavelength: lower index
    CHECK(solve_mode(half, 1.60, Field::Te, Symmetry::Even).n_eff < te.n_eff);

    // cladding index: material index for uniform walls, slab index where a slab reaches the wall
    CHECK_NEAR(cladding_index(half, 1.55, Field::Te, true), 1.444, 1e-12);
    Grid rib = make_grid(60, 222, 0.025, 0.01, 0.0, -1.0, 1.444 * 1.444);
    paint(rib, -10.0, 10.0, 0.0, 0.09, 3.476 * 3.476);
    paint(rib, -0.25, 0.25, 0.0, 0.22, 3.476 * 3.476);
    const double slab = cladding_index(rib, 1.55, Field::Te, true);
    CHECK(slab > 2.0 && slab < 2.2);
    const Mode rib_mode = solve_mode(rib, 1.55, Field::Te, Symmetry::Even);
    CHECK(rib_mode.guided && rib_mode.n_eff > slab);
}
