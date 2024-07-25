// Command line front end.
//
//   modegrid run <structure>     metrics, checks, the permittivity and the field on the grid
//   modegrid check <structure>   metrics and checks; exit status 1 if a check fails
//
// Output lines: `key value`, `check <PASS|FAIL> <name> :: <detail>`, and for `run` one `eps` and one
// `field` line per grid row, bottom row first.
#include <cmath>
#include <cstdio>
#include <exception>
#include <string>

#include "modegrid/analysis.hpp"
#include "modegrid/expect.hpp"
#include "modegrid/parser.hpp"

int main(int argc, char** argv) {
    const std::string command = argc > 1 ? argv[1] : "";
    if (argc != 3 || (command != "run" && command != "check")) {
        std::fprintf(stderr, "usage: modegrid <run|check> <structure file>\n");
        return 2;
    }
    try {
        const mg::Structure s = mg::load_structure(argv[2]);
        const mg::Analysis a = mg::analyse(s);
        std::printf("title %s\n", s.title.c_str());
        std::printf("wavelength %.6f\nnx %d\nny %d\ndx %.6f\ndy %.6f\nx0 %.6f\ny0 %.6f\n", s.wavelength, a.grid.nx, a.grid.ny, a.grid.dx, a.grid.dy, a.grid.x0, a.grid.y0);
        std::printf("neff %.10f\ngroup_index %.8f\nn_cladding %.8f\nguided %d\nconverged %d\niterations %d\n", a.mode.n_eff, a.group_index, a.n_cladding, a.mode.guided ? 1 : 0, a.mode.converged ? 1 : 0, a.mode.iterations);
        std::printf("confinement %.8f\narea %.8f\nrms_x %.8f\nrms_y %.8f\ncentre_y %.8f\n", a.shape.confinement, a.shape.area, a.shape.rms_x, a.shape.rms_y, a.shape.centre_y);
        if (s.coupler) std::printf("n_even %.10f\nn_odd %.10f\nodd_guided %d\nsplitting %.10f\ncoupling_length %.6f\n", a.n_even, a.n_odd, a.odd_guided ? 1 : 0, a.n_even - a.n_odd, a.coupling_length);
        int failed = 0;
        for (const mg::Check& c : mg::evaluate(s, a)) {
            std::printf("check %s %s :: %s\n", c.pass ? "PASS" : "FAIL", c.name.c_str(), c.detail.c_str());
            failed += c.pass ? 0 : 1;
        }
        if (command == "check") return failed == 0 ? 0 : 1;
        for (int iy = 0; iy < a.grid.ny; ++iy) {
            std::printf("eps");
            for (int ix = 0; ix < a.grid.nx; ++ix) std::printf(" %.6g", a.grid.at(ix, iy));
            std::printf("\nfield");
            for (int ix = 0; ix < a.grid.nx; ++ix) std::printf(" %.5e", a.mode.field[static_cast<std::size_t>(iy) * a.grid.nx + ix]);
            std::printf("\n");
        }
        return 0;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "modegrid: %s\n", e.what());
        return 2;
    }
}
