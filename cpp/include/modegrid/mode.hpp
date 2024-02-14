// One guided mode of a grid and the numbers that describe it.
#pragma once

#include <vector>

#include "modegrid/grid.hpp"
#include "modegrid/wave_operator.hpp"

namespace mg {

struct Mode {
    double lambda = 0.0;
    double n_eff = 0.0;
    bool guided = false;         // index above the cladding index, see cladding_index()
    bool converged = false;
    int iterations = 0;
    std::vector<double> field;   // row by row like Grid::eps, unit Euclidean norm, peak positive
};

// The index a mode must exceed to be confined by the structure and not by the walls of the window.
// Top and bottom walls contribute their highest material index. A side wall that cuts through layers (the slab
// of a rib waveguide) contributes the index of the slab mode of that column, found on the same grid; a uniform
// side wall contributes its material index. `symmetric` leaves out the mirror at x = x0.
double cladding_index(const Grid& g, double lambda, Field field, bool symmetric);

// The mode of highest effective index for the given polarisation and symmetry.
Mode solve_mode(const Grid& g, double lambda, Field field, Symmetry symmetry);

struct ModeShape {
    double confinement = 0.0;  // share of the field power in the cells of highest index
    double area = 0.0;         // (integral E^2)^2 / integral E^4, whole structure, um^2
    double rms_x = 0.0;        // root mean square half-width about the centroid, um
    double rms_y = 0.0;
    double centre_y = 0.0;     // centroid height, um
    double peak_x = 0.0, peak_y = 0.0;  // position of the largest field
};

// `symmetric` says the grid holds the right half of a structure mirrored about x = x0.
ModeShape shape_of(const Grid& g, const Mode& mode, bool symmetric);

}  // namespace mg
