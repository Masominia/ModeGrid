#include "modegrid/grid.hpp"

#include <algorithm>
#include <stdexcept>

namespace mg {

Grid make_grid(int nx, int ny, double dx, double dy, double x0, double y0, double eps) {
    if (nx < 2 || ny < 2) throw std::invalid_argument("the grid needs at least two cells each way");
    if (!(dx > 0.0) || !(dy > 0.0)) throw std::invalid_argument("cell sizes must be positive");
    Grid g;
    g.nx = nx, g.ny = ny, g.dx = dx, g.dy = dy, g.x0 = x0, g.y0 = y0;
    g.eps.assign(g.cells(), eps);
    return g;
}

double Grid::max_eps() const { return *std::max_element(eps.begin(), eps.end()); }

double Grid::max_wall_eps(bool symmetric) const {
    double m = 0.0;
    for (int ix = 0; ix < nx; ++ix) m = std::max({m, at(ix, 0), at(ix, ny - 1)});
    for (int iy = 0; iy < ny; ++iy) {
        m = std::max(m, at(nx - 1, iy));
        if (!symmetric) m = std::max(m, at(0, iy));
    }
    return m;
}

static double overlap(double lo, double hi, double a, double b) { return std::max(0.0, std::min(hi, b) - std::max(lo, a)); }

void paint(Grid& g, double xmin, double xmax, double ymin, double ymax, double eps) {
    for (int iy = 0; iy < g.ny; ++iy) {
        const double fy = overlap(g.y0 + iy * g.dy, g.y0 + (iy + 1) * g.dy, ymin, ymax) / g.dy;
        if (fy <= 0.0) continue;
        for (int ix = 0; ix < g.nx; ++ix) {
            const double fx = overlap(g.x0 + ix * g.dx, g.x0 + (ix + 1) * g.dx, xmin, xmax) / g.dx;
            if (fx <= 0.0) continue;
            double f = fx * fy;
            if (f > 1.0 - 1e-9) f = 1.0;  // a cell covered up to rounding is covered
            if (f < 1e-9) continue;
            g.at(ix, iy) = g.at(ix, iy) * (1.0 - f) + f * eps;
        }
    }
}

}  // namespace mg
