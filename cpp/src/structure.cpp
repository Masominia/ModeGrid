#include "modegrid/structure.hpp"

#include <cmath>
#include <stdexcept>

namespace mg {

static int whole_cells(double length, double step, const char* what) {
    const double count = length / step;
    const int n = static_cast<int>(std::lround(count));
    if (std::fabs(count - n) > 1e-6 || n < 2) throw std::invalid_argument(std::string("the window is not a whole number of cells in ") + what);
    return n;
}

Grid grid_of(const Structure& s, double lambda, bool half) {
    const double span = half ? 0.5 * s.width : s.width;
    const int nx = whole_cells(span, s.dx, "x"), ny = whole_cells(s.ymax - s.ymin, s.dy, "y");
    const double nb = s.background.index(lambda);
    Grid g = make_grid(nx, ny, s.dx, s.dy, half ? 0.0 : -0.5 * s.width, s.ymin, nb * nb);
    for (const Rect& r : s.rects) {
        const double n = r.material.index(lambda);
        paint(g, r.xmin, r.xmax, r.ymin, r.ymax, n * n);
    }
    return g;
}

}  // namespace mg
