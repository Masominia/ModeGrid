// The discretised wave operator whose eigenvalues are the squared propagation constants.
#pragma once

#include <cstddef>

#include "modegrid/banded.hpp"
#include "modegrid/grid.hpp"

namespace mg {

// Te: the dominant field is E_x, discontinuous across vertical interfaces.
// Tm: the dominant field is E_y, discontinuous across horizontal interfaces.
// Scalar: no polarisation, continuous field.
enum class Field { Te, Tm, Scalar };

// Cells are numbered along the shorter side of the grid so that the band is as narrow as possible.
struct Ordering {
    int nx, ny;
    bool by_rows;  // true: index = iy * nx + ix
    explicit Ordering(const Grid& g) : nx(g.nx), ny(g.ny), by_rows(g.nx <= g.ny) {}
    std::size_t index(int ix, int iy) const { return by_rows ? static_cast<std::size_t>(iy) * nx + ix : static_cast<std::size_t>(ix) * ny + iy; }
    std::size_t half_width() const { return static_cast<std::size_t>(by_rows ? nx : ny); }
};

// Assembles the operator at the given wavelength. The outer walls are metallic (zero field half a cell
// beyond the last centre); with a symmetry, the wall at x = x0 is a mirror: even or odd field.
Banded assemble(const Grid& g, double lambda, Field field, Symmetry symmetry);

}  // namespace mg
