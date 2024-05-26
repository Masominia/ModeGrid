// Reader for structure files.
//
//   title <text>
//   wavelength <um>
//   grid dx=<um> dy=<um>
//   window width=<um> ymin=<um> ymax=<um>
//   background <material>
//   layer <material> <ymin> <ymax>                 a slab across the whole window
//   rect <material> <xmin> <xmax> <ymin> <ymax>
//   symmetry none|even|odd                         mirror at x = 0; only x >= 0 is solved
//   field te|tm|scalar
//   coupler                                        solve the even and the odd supermode
//   expect <kind> <args...>
//
// Materials: si, sio2, sin, air or a number. `#` starts a comment.
#pragma once

#include <istream>
#include <string>

#include "modegrid/structure.hpp"

namespace mg {

// Throws std::invalid_argument with the line number on any error.
Structure parse_structure(std::istream& in);
Structure parse_structure_text(const std::string& text);
Structure load_structure(const std::string& path);

}  // namespace mg
