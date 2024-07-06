#include "check.hpp"
#include "modegrid/expect.hpp"
#include "modegrid/parser.hpp"

static std::vector<mg::Check> run(const std::string& text) {
    const mg::Structure s = mg::parse_structure_text("grid dx=0.05 dy=0.04\nwindow width=2.4 ymin=-0.8 ymax=1.04\n" + text);
    return mg::evaluate(s, mg::analyse(s));
}

void test_expect(Suite& suite) {
    const std::string strip = "rect si -0.25 0.25 0 0.24\nsymmetry even\n";
    const auto c = run(strip + "expect neff 2.5 0.2\nexpect neff 1.0 0.01\nexpect guided yes\nexpect guided no\nexpect cells 1104\n");
    CHECK(c.size() == 5);
    CHECK(c[0].pass && !c[1].pass && c[2].pass && !c[3].pass && c[4].pass);
    CHECK(c[0].name == "neff 2.5 0.2");
    CHECK(c[4].detail.find("24 by 46") != std::string::npos);

    const auto all = run(strip + "expect neff_between 1.444 3.0\nexpect group_index 4 0.6\nexpect confinement 0.7 0.2\nexpect area 0.15 0.5\n"
                                 "expect rms_x 0.15 0.1\nexpect rms_y 0.1 0.1\nexpect centre_y 0.12 0.05\n");
    for (const auto& check : all) CHECK(check.pass);

    const auto pair = run("rect si 0.1 0.55 0 0.24\nrect si -0.55 -0.1 0 0.24\ncoupler\nexpect n_even 2.5 0.2\nexpect n_odd 2.5 0.2\nexpect splitting 0.03 0.9\nexpect coupling_length 25 0.9\n");
    for (const auto& check : pair) CHECK(check.pass);

    const auto bad = run(strip + "expect neff 2.5\nexpect wobble 1 2\nexpect neff abc 1\nexpect guided maybe\nexpect n_odd 2 1\n");
    CHECK(bad.size() == 5);
    for (const auto& check : bad) CHECK(!check.pass);
    CHECK(bad[0].detail.find("line 5") != std::string::npos && bad[0].detail.find("takes 2 arguments") != std::string::npos);
    CHECK(bad[1].detail.find("unknown expectation") != std::string::npos);
    CHECK(bad[2].detail.find("not a number") != std::string::npos);
    CHECK(bad[3].detail.find("yes or no") != std::string::npos);
    CHECK(bad[4].detail.find("needs a coupler") != std::string::npos);
}
