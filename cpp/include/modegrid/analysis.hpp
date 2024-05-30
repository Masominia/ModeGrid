// Everything computed for one structure file.
#pragma once

#include "modegrid/mode.hpp"
#include "modegrid/structure.hpp"

namespace mg {

struct Analysis {
    Grid grid;                 // at the design wavelength
    Mode mode;                 // the mode asked for (the even supermode of a coupler)
    ModeShape shape;
    double group_index = 0.0;  // n - lambda dn/dlambda by a central difference over wavelength
    double n_cladding = 0.0;   // the index a confined mode must exceed
    // couplers only
    double n_even = 0.0, n_odd = 0.0;
    bool odd_guided = false;
    double coupling_length = 0.0;  // lambda / (2 (n_even - n_odd)), um
};

// Half-step of the wavelength difference used for the group index, um.
constexpr double kGroupStep = 0.005;

Analysis analyse(const Structure& s);

}  // namespace mg
