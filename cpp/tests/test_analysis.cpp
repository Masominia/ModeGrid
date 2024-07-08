#include <cmath>

#include "check.hpp"
#include "modegrid/analysis.hpp"
#include "modegrid/parser.hpp"

void test_analysis(Suite& suite) {
    using namespace mg;
    const std::string head = "grid dx=0.025 dy=0.02\nwindow width=2.4 ymin=-0.8 ymax=1.02\n";
    const Analysis a = analyse(parse_structure_text(head + "rect si -0.225 0.225 0 0.22\nsymmetry even\n"));
    CHECK(a.mode.guided && a.mode.converged);
    // group index well above the phase index, as for any high-contrast strip
    CHECK(a.group_index > a.mode.n_eff + 1.2 && a.group_index < 4.6);
    CHECK_NEAR(a.n_cladding, 1.444, 1e-12);
    CHECK(a.grid.nx == 48 && a.grid.ny == 91);

    // without dispersion in the materials the group index is lower
    const Analysis flat = analyse(parse_structure_text(head + "background 1.444\nrect 3.476 -0.225 0.225 0 0.22\nsymmetry even\n"));
    CHECK_NEAR(flat.mode.n_eff, a.mode.n_eff, 1e-9);
    CHECK(flat.group_index < a.group_index - 0.05);

    // coupler: even above odd, both around the single-guide index, splitting falls with the gap
    auto coupler = [&](double gap) {
        char text[300];
        std::snprintf(text, sizeof text, "grid dx=0.025 dy=0.02\nwindow width=3.6 ymin=-0.8 ymax=1.02\nrect si %g %g 0 0.22\nrect si %g %g 0 0.22\ncoupler\n", gap / 2, gap / 2 + 0.45, -gap / 2 - 0.45, -gap / 2);
        return analyse(parse_structure_text(text));
    };
    const Analysis near = coupler(0.2), far = coupler(0.3);
    CHECK(near.n_even > near.n_odd && near.odd_guided);
    CHECK(near.n_even > a.mode.n_eff && near.n_odd < a.mode.n_eff);
    CHECK_NEAR(near.coupling_length, 1.55 / (2.0 * (near.n_even - near.n_odd)), 1e-9);
    CHECK(far.coupling_length > 1.8 * near.coupling_length);
    CHECK_NEAR(0.5 * (far.n_even + far.n_odd), a.mode.n_eff, 2e-3);

    CHECK_THROWS(analyse(parse_structure_text("window width=3.01 ymin=-1 ymax=1\n")));
}
