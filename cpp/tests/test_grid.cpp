#include "check.hpp"
#include "modegrid/grid.hpp"

void test_grid(Suite& suite) {
    using namespace mg;
    Grid g = make_grid(10, 8, 0.1, 0.05, -0.5, 0.0, 2.0);
    CHECK(g.cells() == 80);
    CHECK_NEAR(g.x(0), -0.45, 1e-12);
    CHECK_NEAR(g.y(7), 0.375, 1e-12);
    CHECK_THROWS(make_grid(1, 8, 0.1, 0.1, 0, 0, 1));
    CHECK_THROWS(make_grid(4, 8, 0.0, 0.1, 0, 0, 1));

    // an aligned rectangle replaces whole cells
    paint(g, -0.2, 0.2, 0.1, 0.2, 9.0);
    CHECK_NEAR(g.at(3, 2), 9.0, 1e-12);
    CHECK_NEAR(g.at(6, 3), 9.0, 1e-12);
    CHECK_NEAR(g.at(2, 2), 2.0, 1e-12);
    CHECK_NEAR(g.at(3, 4), 2.0, 1e-12);
    CHECK_NEAR(g.max_eps(), 9.0, 1e-12);

    // a rectangle covering a quarter of a cell mixes in a quarter of its permittivity
    Grid h = make_grid(4, 4, 1.0, 1.0, 0.0, 0.0, 1.0);
    paint(h, 0.5, 1.0, 0.5, 1.0, 5.0);
    CHECK_NEAR(h.at(0, 0), 2.0, 1e-12);
    // painting is ordered: the later rectangle wins where it covers
    paint(h, 0.0, 4.0, 0.0, 1.0, 3.0);
    CHECK_NEAR(h.at(0, 0), 3.0, 1e-12);
    // a rectangle outside the window changes nothing
    paint(h, 10.0, 11.0, 0.0, 1.0, 7.0);
    CHECK_NEAR(h.max_eps(), 3.0, 1e-12);

    // walls: the mirror side is left out when asked
    Grid w = make_grid(4, 4, 1.0, 1.0, 0.0, 0.0, 1.0);
    w.at(0, 2) = 6.0;
    CHECK_NEAR(w.max_wall_eps(false), 6.0, 1e-12);
    CHECK_NEAR(w.max_wall_eps(true), 1.0, 1e-12);
}
