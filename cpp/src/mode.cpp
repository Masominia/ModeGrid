#include "modegrid/mode.hpp"

#include <cmath>

#include "modegrid/eigen.hpp"

namespace mg {

static double column_index(const Grid& g, int ix, double lambda, Field field) {
    double lo = g.at(ix, 0), hi = lo;
    for (int iy = 1; iy < g.ny; ++iy) lo = std::fmin(lo, g.at(ix, iy)), hi = std::fmax(hi, g.at(ix, iy));
    if (hi - lo < 1e-12) return std::sqrt(hi);
    // one-dimensional version of the operator along y, walls at both ends
    const double pi = std::acos(-1.0);
    const double k = 2.0 * pi / lambda, h2 = g.dy * g.dy;
    const bool jump = field == Field::Tm;
    Banded a(static_cast<std::size_t>(g.ny), 1);
    for (int iy = 0; iy < g.ny; ++iy) {
        const double p = g.at(ix, iy);
        double diagonal = k * k * p;
        for (int jy : {iy - 1, iy + 1}) {
            const bool inside = jy >= 0 && jy < g.ny;
            const double q = inside ? g.at(ix, jy) : p;
            const double off = (jump ? 2.0 * q / (p + q) : 1.0) / h2;
            diagonal -= (jump ? 2.0 * p / (p + q) : 1.0) / h2;
            if (inside) a(iy, jy) = off;
            else diagonal -= off;
        }
        a(iy, iy) = diagonal;
    }
    const EigenPair pair = largest_eigenpair(a, 1.0005 * k * k * hi);
    return pair.value > 0.0 ? std::sqrt(pair.value) / k : 0.0;
}

double cladding_index(const Grid& g, double lambda, Field field, bool symmetric) {
    double eps = 0.0;
    for (int ix = 0; ix < g.nx; ++ix) eps = std::fmax(eps, std::fmax(g.at(ix, 0), g.at(ix, g.ny - 1)));
    double n = std::fmax(std::sqrt(eps), column_index(g, g.nx - 1, lambda, field));
    if (!symmetric) n = std::fmax(n, column_index(g, 0, lambda, field));
    return n;
}

Mode solve_mode(const Grid& g, double lambda, Field field, Symmetry symmetry) {
    const double pi = std::acos(-1.0);
    const double k = 2.0 * pi / lambda;
    const Banded a = assemble(g, lambda, field, symmetry);
    // every eigenvalue lies below k^2 eps_max; a small margin keeps the shifted matrix definite
    const EigenPair pair = largest_eigenpair(a, 1.0005 * k * k * g.max_eps());
    Mode m;
    m.lambda = lambda;
    m.n_eff = pair.value > 0.0 ? std::sqrt(pair.value) / k : 0.0;
    m.converged = pair.converged;
    m.iterations = pair.iterations;
    m.guided = m.n_eff > cladding_index(g, lambda, field, symmetry != Symmetry::None) + 1e-9;
    const Ordering order(g);
    m.field.resize(g.cells());
    for (int iy = 0; iy < g.ny; ++iy)
        for (int ix = 0; ix < g.nx; ++ix) m.field[static_cast<std::size_t>(iy) * g.nx + ix] = pair.vector[order.index(ix, iy)];
    return m;
}

ModeShape shape_of(const Grid& g, const Mode& mode, bool symmetric) {
    const double top = g.max_eps() * (1.0 - 1e-9);
    double s2 = 0.0, s4 = 0.0, core = 0.0, sx = 0.0, sy = 0.0, sxx = 0.0, syy = 0.0, peak = -1.0;
    ModeShape out;
    for (int iy = 0; iy < g.ny; ++iy) {
        for (int ix = 0; ix < g.nx; ++ix) {
            const double e = mode.field[static_cast<std::size_t>(iy) * g.nx + ix], p = e * e;
            s2 += p, s4 += p * p;
            if (g.at(ix, iy) >= top) core += p;
            sx += p * g.x(ix), sy += p * g.y(iy);
            sxx += p * g.x(ix) * g.x(ix), syy += p * g.y(iy) * g.y(iy);
            if (std::fabs(e) > peak) peak = std::fabs(e), out.peak_x = g.x(ix), out.peak_y = g.y(iy);
        }
    }
    out.confinement = core / s2;
    out.area = (symmetric ? 2.0 : 1.0) * s2 * s2 / s4 * g.dx * g.dy;
    // a mirrored structure has its centroid on the mirror
    const double cx = symmetric ? g.x0 : sx / s2;
    out.centre_y = sy / s2;
    out.rms_x = std::sqrt(std::fmax(0.0, sxx / s2 - 2.0 * cx * sx / s2 + cx * cx));
    out.rms_y = std::sqrt(std::fmax(0.0, syy / s2 - out.centre_y * out.centre_y));
    return out;
}

}  // namespace mg
