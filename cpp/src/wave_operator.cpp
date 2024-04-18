#include "modegrid/wave_operator.hpp"

#include <cmath>

namespace mg {

Banded assemble(const Grid& g, double lambda, Field field, Symmetry symmetry) {
    const double pi = std::acos(-1.0);
    const double k2 = std::pow(2.0 * pi / lambda, 2);
    const Ordering order(g);
    Banded a(g.cells(), order.half_width());
    const int step[4][2] = {{-1, 0}, {1, 0}, {0, -1}, {0, 1}};
    for (int iy = 0; iy < g.ny; ++iy) {
        for (int ix = 0; ix < g.nx; ++ix) {
            const double p = g.at(ix, iy);
            const std::size_t me = order.index(ix, iy);
            double diagonal = k2 * p;
            for (int s = 0; s < 4; ++s) {
                const int jx = ix + step[s][0], jy = iy + step[s][1];
                const bool along_x = step[s][1] == 0;
                const double h = along_x ? g.dx : g.dy;
                const bool inside = jx >= 0 && jx < g.nx && jy >= 0 && jy < g.ny;
                const double q = inside ? g.at(jx, jy) : p;
                // across the direction of the dominant field the normal flux is continuous, not the field
                const bool jump = (field == Field::Te && along_x) || (field == Field::Tm && !along_x);
                const double off = (jump ? 2.0 * q / (p + q) : 1.0) / (h * h);
                diagonal -= (jump ? 2.0 * p / (p + q) : 1.0) / (h * h);
                if (inside) {
                    a(me, order.index(jx, jy)) = off;
                } else if (jx < 0 && symmetry == Symmetry::Even) {
                    diagonal += off;  // mirror image with the same sign
                } else {
                    diagonal -= off;  // wall or odd mirror: image with the opposite sign
                }
            }
            a(me, me) = diagonal;
        }
    }
    return a;
}

}  // namespace mg
