#include "modegrid/analysis.hpp"

#include <cmath>

namespace mg {

Analysis analyse(const Structure& s) {
    const Symmetry symmetry = s.coupler ? Symmetry::Even : s.symmetry;
    const bool half = symmetry != Symmetry::None;
    Analysis a;
    a.grid = grid_of(s, s.wavelength, half);
    a.mode = solve_mode(a.grid, s.wavelength, s.field, symmetry);
    a.shape = shape_of(a.grid, a.mode, half);
    a.n_cladding = cladding_index(a.grid, s.wavelength, s.field, half);
    const double lo = s.wavelength - kGroupStep, hi = s.wavelength + kGroupStep;
    const double n_lo = solve_mode(grid_of(s, lo, half), lo, s.field, symmetry).n_eff;
    const double n_hi = solve_mode(grid_of(s, hi, half), hi, s.field, symmetry).n_eff;
    a.group_index = a.mode.n_eff - s.wavelength * (n_hi - n_lo) / (2.0 * kGroupStep);
    if (s.coupler) {
        const Mode odd = solve_mode(a.grid, s.wavelength, s.field, Symmetry::Odd);
        a.n_even = a.mode.n_eff, a.n_odd = odd.n_eff, a.odd_guided = odd.guided;
        a.coupling_length = s.wavelength / (2.0 * std::fabs(a.n_even - a.n_odd));
    }
    return a;
}

}  // namespace mg
