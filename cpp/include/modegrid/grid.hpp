// A uniform grid of rectangular cells, each with one permittivity. Fields live at the cell centres.
#pragma once

#include <cstddef>
#include <vector>

namespace mg {

enum class Symmetry { None, Even, Odd };

struct Grid {
    int nx = 0, ny = 0;
    double dx = 0.0, dy = 0.0;
    double x0 = 0.0, y0 = 0.0;  // lower left corner of the window
    std::vector<double> eps;    // row by row, ny rows of nx cells

    std::size_t cells() const { return static_cast<std::size_t>(nx) * ny; }
    double& at(int ix, int iy) { return eps[static_cast<std::size_t>(iy) * nx + ix]; }
    double at(int ix, int iy) const { return eps[static_cast<std::size_t>(iy) * nx + ix]; }
    double x(int ix) const { return x0 + (ix + 0.5) * dx; }
    double y(int iy) const { return y0 + (iy + 0.5) * dy; }

    double max_eps() const;
    // Largest permittivity among the cells on the outer walls. `symmetric` leaves out the wall at x = x0.
    double max_wall_eps(bool symmetric) const;
};

// A window filled with one permittivity. Throws std::invalid_argument on a non-positive size.
Grid make_grid(int nx, int ny, double dx, double dy, double x0, double y0, double eps);

// Paints a rectangle: each cell takes the new permittivity in proportion to the area the rectangle covers.
void paint(Grid& g, double xmin, double xmax, double ymin, double ymax, double eps);

}  // namespace mg
