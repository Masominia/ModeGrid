#include "modegrid/material.hpp"

#include <stdexcept>

#include "modegrid/strutil.hpp"

namespace mg {

Material material_by_name(const std::string& text) {
    if (text == "si" || text == "silicon") return {3.476, -0.0865};
    if (text == "sio2" || text == "silica" || text == "oxide") return {1.444, -0.0120};
    if (text == "sin" || text == "nitride") return {1.996, -0.0240};
    if (text == "air") return {1.0, 0.0};
    double n = 0.0;
    try {
        n = to_number(text);
    } catch (const std::invalid_argument&) {
        throw std::invalid_argument("'" + text + "' is neither a known material nor an index");
    }
    if (n < 1.0) throw std::invalid_argument("an index must be at least 1, got " + text);
    return {n, 0.0};
}

}  // namespace mg
