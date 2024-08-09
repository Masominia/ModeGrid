#include "check.hpp"
#include "modegrid/material.hpp"

void test_material(Suite& suite) {
    using namespace mg;
    CHECK_NEAR(material_by_name("si").index(1.55), 3.476, 1e-12);
    CHECK_NEAR(material_by_name("oxide").index(1.55), 1.444, 1e-12);
    CHECK_NEAR(material_by_name("nitride").index(1.55), 1.996, 1e-12);
    CHECK_NEAR(material_by_name("air").index(1.3), 1.0, 1e-12);
    CHECK(material_by_name("si").index(1.5) > material_by_name("si").index(1.6));
    CHECK_NEAR(material_by_name("si").group_index(1.55), 3.476 + 1.55 * 0.0865, 1e-12);
    const Material constant = material_by_name("2.1");
    CHECK_NEAR(constant.index(1.2), 2.1, 1e-12);
    CHECK_NEAR(constant.group_index(1.2), 2.1, 1e-12);
    CHECK_THROWS(material_by_name("mithril"));
    CHECK_THROWS(material_by_name("0.8"));
}
