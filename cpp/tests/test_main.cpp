#include <cstdio>

#include "check.hpp"

int main() {
    Suite suite;
    struct Entry {
        const char* name;
        void (*run)(Suite&);
    };
    const Entry entries[] = {{"strutil", test_strutil}, {"material", test_material}, {"grid", test_grid}, {"banded", test_banded}, {"eigen", test_eigen}, {"wave_operator", test_wave_operator},
                             {"mode", test_mode}, {"slab", test_slab}, {"structure", test_structure}, {"parser", test_parser}, {"analysis", test_analysis}, {"expect", test_expect}};
    for (const Entry& e : entries) {
        const int before = suite.failures, checks = suite.checks;
        e.run(suite);
        std::printf("%-14s %3d checks  %s\n", e.name, suite.checks - checks, suite.failures == before ? "ok" : "FAILED");
    }
    std::printf("%d checks, %d failures\n", suite.checks, suite.failures);
    return suite.failures == 0 ? 0 : 1;
}
