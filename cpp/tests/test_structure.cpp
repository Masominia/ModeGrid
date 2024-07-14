#include <cmath>

#include "check.hpp"
#include "modegrid/structure.hpp"

void test_structure(Suite& suite) {
    using namespace mg;
    Structure s;
    s.dx = 0.025, s.dy = 0.02;
    s.rects.push_back({material_by_name("si"), -0.225, 0.225, 0.0, 0.22});
    const Grid full = grid_of(s, 1.55, false), half = grid_of(s, 1.55, true);
    CHECK(full.nx == 120 && full.ny == 111 && half.nx == 60);
    CHECK_NEAR(full.x0, -1.5, 1e-12);
    CHECK_NEAR(half.x0, 0.0, 1e-12);
    CHECK_NEAR(full.y0, -1.0, 1e-12);
    // the core occupies 18 by 11 cells
    int core = 0;
    for (double e : full.eps) core += e > 12.0 ? 1 : 0;
    CHECK(core == 18 * 11);
    CHECK_NEAR(full.max_eps(), 3.476 * 3.476, 1e-12);
    CHECK_NEAR(half.at(0, 55), 3.476 * 3.476, 1e-12);
    CHECK_NEAR(half.at(9, 55), 1.444 * 1.444, 1e-12);

    // dispersion enters through the wavelength
    CHECK(grid_of(s, 1.50, true).max_eps() > half.max_eps());

    // a window that is not a whole number of cells is refused
    s.width = 3.01;
    CHECK_THROWS(grid_of(s, 1.55, false));
    s.width = 3.0, s.ymax = 1.23;
    CHECK_THROWS(grid_of(s, 1.55, false));
}
