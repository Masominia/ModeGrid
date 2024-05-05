#include "modegrid/eigen.hpp"

#include <cmath>

namespace mg {

static double dot(const std::vector<double>& a, const std::vector<double>& b) {
    double s = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) s += a[i] * b[i];
    return s;
}

static void normalise(std::vector<double>& x) {
    const double norm = std::sqrt(dot(x, x));
    for (double& v : x) v /= norm;
}

EigenPair largest_eigenpair(const Banded& a, double shift, double tolerance, int max_iterations) {
    const std::size_t n = a.size();
    EigenPair out;
    // a start with no sign changes has a component along the fundamental mode of any polarisation
    out.vector.assign(n, 1.0);
    for (std::size_t i = 0; i < n; ++i) out.vector[i] += 0.3 * std::sin(0.37 * static_cast<double>(i));
    normalise(out.vector);

    double sigma = shift, previous = 0.0;
    bool moved = false;
    Banded lu = a;
    lu.add_to_diagonal(-sigma);
    lu.factorise();
    for (int it = 1; it <= max_iterations; ++it) {
        std::vector<double> y = out.vector;
        lu.solve(y);
        const double mu = dot(out.vector, y);  // Rayleigh quotient of the inverse
        const double value = sigma + 1.0 / mu;
        normalise(y);
        out.vector.swap(y);
        out.value = value;
        out.iterations = it;
        if (it > 1 && std::fabs(value - previous) <= tolerance * std::fabs(value)) {
            out.converged = true;
            break;
        }
        if (!moved && it >= 6 && std::fabs(value - previous) <= 1e-4 * std::fabs(shift - value)) {
            // close enough to trust the estimate: bring the shift to a twentieth of its distance
            sigma = value + 0.05 * (shift - value);
            lu = a;
            lu.add_to_diagonal(-sigma);
            lu.factorise();
            moved = true;
        }
        previous = value;
    }
    double peak = 0.0;
    for (double v : out.vector)
        if (std::fabs(v) > std::fabs(peak)) peak = v;
    if (peak < 0.0)
        for (double& v : out.vector) v = -v;
    return out;
}

}  // namespace mg
