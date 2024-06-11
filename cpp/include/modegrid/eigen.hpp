// The largest eigenvalue of the wave operator by shift-and-invert iteration.
#pragma once

#include <vector>

#include "modegrid/banded.hpp"

namespace mg {

struct EigenPair {
    double value = 0.0;
    std::vector<double> vector;  // unit Euclidean norm, largest entry positive
    int iterations = 0;
    bool converged = false;
};

// Iterates x <- (A - shift)^-1 x from a smooth start. `shift` must lie above every eigenvalue, so that the
// largest one is the closest; the shift is moved once towards the estimate to speed up convergence.
// Stops when the eigenvalue changes by less than `tolerance` relative to its size.
EigenPair largest_eigenpair(const Banded& a, double shift, double tolerance = 1e-13, int max_iterations = 300);

}  // namespace mg
