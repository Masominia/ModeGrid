#include "check.hpp"
#include "modegrid/strutil.hpp"

void test_strutil(Suite& suite) {
    using namespace mg;
    CHECK(split_words("  rect si  0 1 ").size() == 4);
    CHECK(strip_comment("  grid dx=0.02 # fine ") == "grid dx=0.02");
    CHECK(strip_comment("# only a comment").empty());
    CHECK_NEAR(to_number("-1.5e-1"), -0.15, 1e-15);
    CHECK_THROWS(to_number("12abc"));
    CHECK_THROWS(to_number(""));
    const Params p = parse_params({"grid", "dx=0.02", "dy=0.01"}, 1);
    CHECK_NEAR(get(p, "dx", 9.0), 0.02, 1e-15);
    CHECK_NEAR(get(p, "missing", 9.0), 9.0, 1e-15);
    CHECK_NEAR(require(p, "dy", "grid"), 0.01, 1e-15);
    CHECK_THROWS(require(p, "dz", "grid"));
    CHECK_THROWS(parse_params({"grid", "dx"}, 1));
    CHECK_THROWS(parse_params({"grid", "=3"}, 1));
}
