#include "modegrid/parser.hpp"

#include <fstream>
#include <sstream>
#include <stdexcept>

#include "modegrid/strutil.hpp"

namespace mg {

static void need(const std::vector<std::string>& w, std::size_t count, const char* usage) {
    if (w.size() != count) throw std::invalid_argument(std::string("usage: ") + usage);
}

static void statement(Structure& s, const std::string& line, const std::vector<std::string>& w, int line_number) {
    const std::string& key = w[0];
    if (key == "title") {
        s.title = line.size() > 6 ? line.substr(6) : "";
    } else if (key == "wavelength") {
        need(w, 2, "wavelength <um>");
        s.wavelength = to_number(w[1]);
        if (!(s.wavelength > 0.0)) throw std::invalid_argument("wavelength must be positive");
    } else if (key == "grid") {
        const Params p = parse_params(w, 1);
        s.dx = get(p, "dx", s.dx), s.dy = get(p, "dy", s.dy);
        if (!(s.dx > 0.0) || !(s.dy > 0.0)) throw std::invalid_argument("cell sizes must be positive");
    } else if (key == "window") {
        const Params p = parse_params(w, 1);
        s.width = get(p, "width", s.width), s.ymin = get(p, "ymin", s.ymin), s.ymax = get(p, "ymax", s.ymax);
        if (!(s.width > 0.0) || !(s.ymax > s.ymin)) throw std::invalid_argument("the window must have a positive size");
    } else if (key == "background") {
        need(w, 2, "background <material>");
        s.background = material_by_name(w[1]);
    } else if (key == "layer") {
        need(w, 4, "layer <material> <ymin> <ymax>");
        const Rect r{material_by_name(w[1]), -1e9, 1e9, to_number(w[2]), to_number(w[3])};
        if (!(r.ymax > r.ymin)) throw std::invalid_argument("a layer needs ymin < ymax");
        s.rects.push_back(r);
    } else if (key == "rect") {
        need(w, 6, "rect <material> <xmin> <xmax> <ymin> <ymax>");
        const Rect r{material_by_name(w[1]), to_number(w[2]), to_number(w[3]), to_number(w[4]), to_number(w[5])};
        if (!(r.xmax > r.xmin) || !(r.ymax > r.ymin)) throw std::invalid_argument("a rectangle needs xmin < xmax and ymin < ymax");
        s.rects.push_back(r);
    } else if (key == "symmetry") {
        need(w, 2, "symmetry none|even|odd");
        if (w[1] == "none") s.symmetry = Symmetry::None;
        else if (w[1] == "even") s.symmetry = Symmetry::Even;
        else if (w[1] == "odd") s.symmetry = Symmetry::Odd;
        else throw std::invalid_argument("symmetry must be none, even or odd");
    } else if (key == "field") {
        need(w, 2, "field te|tm|scalar");
        if (w[1] == "te") s.field = Field::Te;
        else if (w[1] == "tm") s.field = Field::Tm;
        else if (w[1] == "scalar") s.field = Field::Scalar;
        else throw std::invalid_argument("field must be te, tm or scalar");
    } else if (key == "coupler") {
        need(w, 1, "coupler");
        s.coupler = true;
    } else if (key == "expect") {
        if (w.size() < 2) throw std::invalid_argument("expect needs a kind");
        s.expectations.push_back({w[1], std::vector<std::string>(w.begin() + 2, w.end()), line_number});
    } else {
        throw std::invalid_argument("unknown statement '" + key + "'");
    }
}

Structure parse_structure(std::istream& in) {
    Structure s;
    std::string raw;
    int line_number = 0;
    while (std::getline(in, raw)) {
        ++line_number;
        const std::string line = strip_comment(raw);
        if (line.empty()) continue;
        try {
            statement(s, line, split_words(line), line_number);
        } catch (const std::invalid_argument& e) {
            throw std::invalid_argument("line " + std::to_string(line_number) + ": " + e.what());
        }
    }
    if (s.coupler && s.symmetry != Symmetry::None) throw std::invalid_argument("a coupler solves both symmetries itself: remove the symmetry statement");
    return s;
}

Structure parse_structure_text(const std::string& text) {
    std::istringstream in(text);
    return parse_structure(in);
}

Structure load_structure(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::invalid_argument("cannot read " + path);
    return parse_structure(in);
}

}  // namespace mg
