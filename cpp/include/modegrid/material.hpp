// Refractive indices with first-order dispersion around 1.55 um.
#pragma once

#include <string>

namespace mg {

struct Material {
    double n0 = 1.0;      // index at the reference wavelength
    double slope = 0.0;   // dn / dlambda per micrometre
    double index(double lambda) const { return n0 + slope * (lambda - 1.55); }
    double group_index(double lambda) const { return index(lambda) - lambda * slope; }
};

// si, sio2, sin, air (and their long names), or a number for a constant index of at least 1.
// Throws std::invalid_argument otherwise.
Material material_by_name(const std::string& text);

}  // namespace mg
