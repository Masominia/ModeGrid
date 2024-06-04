#include "modegrid/expect.hpp"

#include <cmath>
#include <cstdio>
#include <stdexcept>

#include "modegrid/strutil.hpp"

namespace mg {

static std::string fmt(const char* pattern, double a, double b, double c) {
    char buffer[200];
    std::snprintf(buffer, sizeof buffer, pattern, a, b, c);
    return buffer;
}

static std::vector<double> numbers(const Expectation& e, std::size_t count) {
    if (e.args.size() != count) throw std::invalid_argument("'" + e.kind + "' takes " + std::to_string(count) + " argument" + (count == 1 ? "" : "s"));
    std::vector<double> v;
    for (const std::string& a : e.args) v.push_back(to_number(a));
    return v;
}

static Check absolute(double got, const Expectation& e) {
    const std::vector<double> v = numbers(e, 2);
    return {"", std::fabs(got - v[0]) <= v[1], fmt("%.8f, expected %.8f, off by %.2e", got, v[0], got - v[0])};
}

static Check relative(double got, const Expectation& e) {
    const std::vector<double> v = numbers(e, 2);
    return {"", std::fabs(got / v[0] - 1.0) <= v[1], fmt("%.6g, expected %.6g, relative error %.2e", got, v[0], got / v[0] - 1.0)};
}

static void need_coupler(const Structure& s, const Expectation& e) {
    if (!s.coupler) throw std::invalid_argument("'" + e.kind + "' needs a coupler statement");
}

static Check one(const Structure& s, const Analysis& a, const Expectation& e) {
    const std::string& k = e.kind;
    if (k == "neff") return absolute(a.mode.n_eff, e);
    if (k == "group_index") return absolute(a.group_index, e);
    if (k == "confinement") return absolute(a.shape.confinement, e);
    if (k == "rms_x") return absolute(a.shape.rms_x, e);
    if (k == "rms_y") return absolute(a.shape.rms_y, e);
    if (k == "centre_y") return absolute(a.shape.centre_y, e);
    if (k == "area") return relative(a.shape.area, e);
    if (k == "neff_between") {
        const std::vector<double> v = numbers(e, 2);
        return {"", a.mode.n_eff >= v[0] && a.mode.n_eff <= v[1], fmt("%.6f, expected between %.6f and %.6f", a.mode.n_eff, v[0], v[1])};
    }
    if (k == "guided") {
        if (e.args.size() != 1 || (e.args[0] != "yes" && e.args[0] != "no")) throw std::invalid_argument("'guided' takes yes or no");
        return {"", a.mode.guided == (e.args[0] == "yes"), fmt("index %.6f, cladding %.6f, difference %.2e", a.mode.n_eff, a.n_cladding, a.mode.n_eff - a.n_cladding)};
    }
    if (k == "cells") {
        const std::vector<double> v = numbers(e, 1);
        return {"", static_cast<double>(a.grid.cells()) == v[0], fmt("%.0f cells (%.0f by %.0f)", static_cast<double>(a.grid.cells()), a.grid.nx, a.grid.ny)};
    }
    if (k == "n_even") return need_coupler(s, e), absolute(a.n_even, e);
    if (k == "n_odd") return need_coupler(s, e), absolute(a.n_odd, e);
    if (k == "splitting") return need_coupler(s, e), relative(a.n_even - a.n_odd, e);
    if (k == "coupling_length") return need_coupler(s, e), relative(a.coupling_length, e);
    throw std::invalid_argument("unknown expectation '" + k + "'");
}

std::vector<Check> evaluate(const Structure& s, const Analysis& a) {
    std::vector<Check> out;
    for (const Expectation& e : s.expectations) {
        std::string name = e.kind;
        for (const std::string& arg : e.args) name += " " + arg;
        Check c;
        try {
            c = one(s, a, e);
        } catch (const std::invalid_argument& err) {
            c = {"", false, "line " + std::to_string(e.line) + ": " + err.what()};
        }
        c.name = name;
        out.push_back(c);
    }
    return out;
}

}  // namespace mg
