#include "modegrid/strutil.hpp"

#include <sstream>
#include <stdexcept>

namespace mg {

std::vector<std::string> split_words(const std::string& line) {
    std::vector<std::string> words;
    std::istringstream ss(line);
    std::string w;
    while (ss >> w) words.push_back(w);
    return words;
}

std::string strip_comment(const std::string& line) {
    std::string s = line.substr(0, line.find('#'));
    const auto first = s.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return "";
    return s.substr(first, s.find_last_not_of(" \t\r\n") - first + 1);
}

double to_number(const std::string& text) {
    std::size_t used = 0;
    double v = 0.0;
    try {
        v = std::stod(text, &used);
    } catch (const std::exception&) {
        throw std::invalid_argument("'" + text + "' is not a number");
    }
    if (used != text.size()) throw std::invalid_argument("'" + text + "' is not a number");
    return v;
}

Params parse_params(const std::vector<std::string>& words, std::size_t first) {
    Params p;
    for (std::size_t i = first; i < words.size(); ++i) {
        const auto eq = words[i].find('=');
        if (eq == std::string::npos || eq == 0) throw std::invalid_argument("expected key=value, got '" + words[i] + "'");
        p[words[i].substr(0, eq)] = to_number(words[i].substr(eq + 1));
    }
    return p;
}

double get(const Params& p, const std::string& key, double fallback) {
    const auto it = p.find(key);
    return it == p.end() ? fallback : it->second;
}

double require(const Params& p, const std::string& key, const std::string& owner) {
    const auto it = p.find(key);
    if (it == p.end()) throw std::invalid_argument(owner + " needs " + key + "=...");
    return it->second;
}

}  // namespace mg
