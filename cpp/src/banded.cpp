#include "modegrid/banded.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace mg {

Banded::Banded(std::size_t n, std::size_t half_width) : n_(n), m_(std::min(half_width, n ? n - 1 : 0)), data_(n * (2 * m_ + 1), 0.0) {}

void Banded::add_to_diagonal(double value) {
    for (std::size_t i = 0; i < n_; ++i) (*this)(i, i) += value;
}

std::vector<double> Banded::multiply(const std::vector<double>& x) const {
    if (x.size() != n_) throw std::invalid_argument("vector length does not match the matrix");
    std::vector<double> y(n_, 0.0);
    for (std::size_t i = 0; i < n_; ++i) {
        const std::size_t lo = i > m_ ? i - m_ : 0, hi = std::min(n_ - 1, i + m_);
        double s = 0.0;
        for (std::size_t j = lo; j <= hi; ++j) s += (*this)(i, j) * x[j];
        y[i] = s;
    }
    return y;
}

void Banded::factorise() {
    for (std::size_t k = 0; k < n_; ++k) {
        const double pivot = (*this)(k, k);
        if (std::fabs(pivot) < 1e-300) throw std::runtime_error("zero pivot in the banded factorisation");
        const std::size_t last = std::min(n_ - 1, k + m_);
        for (std::size_t i = k + 1; i <= last; ++i) {
            const double factor = (*this)(i, k) / pivot;
            if (factor == 0.0) continue;
            (*this)(i, k) = factor;
            double* row_i = &data_[i * (2 * m_ + 1) + (m_ - i)];  // row_i[j] is element (i, j)
            const double* row_k = &data_[k * (2 * m_ + 1) + (m_ - k)];
            for (std::size_t j = k + 1; j <= last; ++j) row_i[j] -= factor * row_k[j];
        }
    }
    factorised_ = true;
}

void Banded::solve(std::vector<double>& b) const {
    if (!factorised_) throw std::logic_error("solve called before factorise");
    if (b.size() != n_) throw std::invalid_argument("vector length does not match the matrix");
    for (std::size_t i = 1; i < n_; ++i) {
        const std::size_t lo = i > m_ ? i - m_ : 0;
        double s = b[i];
        for (std::size_t j = lo; j < i; ++j) s -= (*this)(i, j) * b[j];
        b[i] = s;
    }
    for (std::size_t i = n_; i-- > 0;) {
        const std::size_t hi = std::min(n_ - 1, i + m_);
        double s = b[i];
        for (std::size_t j = i + 1; j <= hi; ++j) s -= (*this)(i, j) * b[j];
        b[i] = s / (*this)(i, i);
    }
}

}  // namespace mg
