// A waveguide cross-section as a list of painted rectangles, and its grid at a wavelength.
#pragma once

#include <string>
#include <vector>

#include "modegrid/grid.hpp"
#include "modegrid/material.hpp"
#include "modegrid/wave_operator.hpp"

namespace mg {

struct Rect {
    Material material;
    double xmin, xmax, ymin, ymax;
};

struct Expectation {
    std::string kind;
    std::vector<std::string> args;
    int line = 0;
};

struct Structure {
    std::string title;
    double wavelength = 1.55;
    double dx = 0.02, dy = 0.02;
    double width = 3.0;             // full width of the window, centred on x = 0
    double ymin = -1.0, ymax = 1.22;
    Material background{1.444, -0.0120};
    std::vector<Rect> rects;        // painted in order, later ones on top
    Symmetry symmetry = Symmetry::None;
    Field field = Field::Te;
    bool coupler = false;           // solve the even and the odd supermode
    std::vector<Expectation> expectations;
};

// The grid at a wavelength. With a symmetry (or `half` set) only x >= 0 is kept.
// Throws std::invalid_argument if the window is not a whole number of cells.
Grid grid_of(const Structure& s, double lambda, bool half);

}  // namespace mg
