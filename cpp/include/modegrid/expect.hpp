// Expectations of a structure file, checked against its analysis.
//
//   neff <value> <tol>            effective index
//   neff_between <low> <high>
//   group_index <value> <tol>
//   confinement <fraction> <tol>  share of the field power in the highest-index material
//   area <um^2> <reltol>          effective mode area
//   rms_x | rms_y <um> <tol>      root mean square half-widths
//   centre_y <um> <tol>
//   guided yes|no                 index above the cladding index
//   cells <count>                 number of unknowns
//   n_even | n_odd <value> <tol>  supermodes of a coupler
//   splitting <value> <reltol>    n_even - n_odd
//   coupling_length <um> <reltol>
#pragma once

#include <string>
#include <vector>

#include "modegrid/analysis.hpp"

namespace mg {

struct Check {
    std::string name;
    bool pass = false;
    std::string detail;
};

// A malformed expectation becomes a failed check with its line number; it does not stop the others.
std::vector<Check> evaluate(const Structure& s, const Analysis& a);

}  // namespace mg
