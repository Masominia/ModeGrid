#include "check.hpp"
#include "modegrid/parser.hpp"

void test_parser(Suite& suite) {
    using namespace mg;
    const Structure d = parse_structure_text("");
    CHECK_NEAR(d.wavelength, 1.55, 1e-12);
    CHECK(d.symmetry == Symmetry::None && d.field == Field::Te && !d.coupler && d.rects.empty());

    const Structure s = parse_structure_text(
        "# note\ntitle My rib\nwavelength 1.31\ngrid dx=0.025 dy=0.01\nwindow width=4 ymin=-1 ymax=1.2\nbackground air\n"
        "layer sio2 -1 0\nlayer si 0 0.09\nrect 3.5 -0.25 0.25 0 0.22  # core\nsymmetry odd\nfield tm\nexpect neff 2.5 0.1\n");
    CHECK(s.title == "My rib");
    CHECK_NEAR(s.wavelength, 1.31, 1e-12);
    CHECK_NEAR(s.dx, 0.025, 1e-12);
    CHECK_NEAR(s.dy, 0.01, 1e-12);
    CHECK_NEAR(s.width, 4.0, 1e-12);
    CHECK_NEAR(s.ymax, 1.2, 1e-12);
    CHECK_NEAR(s.background.n0, 1.0, 1e-12);
    CHECK(s.rects.size() == 3);
    CHECK(s.rects[0].xmin < -1e8 && s.rects[0].xmax > 1e8);
    CHECK_NEAR(s.rects[2].material.n0, 3.5, 1e-12);
    CHECK_NEAR(s.rects[2].xmax, 0.25, 1e-12);
    CHECK(s.symmetry == Symmetry::Odd && s.field == Field::Tm);
    CHECK(s.expectations.size() == 1 && s.expectations[0].kind == "neff" && s.expectations[0].line == 12);
    CHECK(parse_structure_text("coupler\n").coupler);

    const char* bad[][2] = {
        {"wavelength 1.55\nfrobnicate", "line 2: unknown statement"},
        {"wavelength -1", "must be positive"},
        {"wavelength", "usage: wavelength"},
        {"grid dx=0", "cell sizes"},
        {"grid dx", "key=value"},
        {"window width=3 ymin=1 ymax=0", "positive size"},
        {"background unobtainium", "neither a known material"},
        {"layer si 0.2 0.1", "ymin < ymax"},
        {"rect si 0 1 0", "usage: rect"},
        {"rect si 1 0 0 1", "xmin < xmax"},
        {"rect si a 1 0 1", "not a number"},
        {"symmetry sideways", "none, even or odd"},
        {"field circular", "te, tm or scalar"},
        {"expect", "needs a kind"},
        {"coupler\nsymmetry even", "remove the symmetry"},
    };
    for (const auto& b : bad) CHECK(throws_with([&] { parse_structure_text(b[0]); }, b[1]));
    CHECK_THROWS(load_structure("/no/such/file.wg"));
}
