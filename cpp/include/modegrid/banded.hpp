// A banded square matrix with LU factorisation in place (no pivoting).
#pragma once

#include <cstddef>
#include <vector>

namespace mg {

class Banded {
public:
    Banded(std::size_t n, std::size_t half_width);

    std::size_t size() const { return n_; }
    std::size_t half_width() const { return m_; }

    // Element (i, j); j must lie within the band.
    double& operator()(std::size_t i, std::size_t j) { return data_[i * (2 * m_ + 1) + (j + m_ - i)]; }
    double operator()(std::size_t i, std::size_t j) const { return data_[i * (2 * m_ + 1) + (j + m_ - i)]; }
    bool in_band(std::size_t i, std::size_t j) const { return (i > j ? i - j : j - i) <= m_; }

    void add_to_diagonal(double value);
    std::vector<double> multiply(const std::vector<double>& x) const;

    // Replaces the matrix by its LU factors. Throws std::runtime_error on a zero pivot.
    // Without pivoting this is safe for the definite matrices the solver produces, not in general.
    void factorise();
    // Solves with the factors, overwriting `b`. Call factorise() first.
    void solve(std::vector<double>& b) const;

private:
    std::size_t n_, m_;
    std::vector<double> data_;
    bool factorised_ = false;
};

}  // namespace mg
