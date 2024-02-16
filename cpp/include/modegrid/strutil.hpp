// Small text helpers for the structure parser and the command line.
#pragma once

#include <map>
#include <string>
#include <vector>

namespace mg {

using Params = std::map<std::string, double>;

std::vector<std::string> split_words(const std::string& line);

// Removes everything from the first '#' and surrounding blanks.
std::string strip_comment(const std::string& line);

// Parses `key=value` words into numbers. Throws std::invalid_argument on a malformed word.
Params parse_params(const std::vector<std::string>& words, std::size_t first);

double to_number(const std::string& text);

// Value of `key`, or `fallback` if absent.
double get(const Params& p, const std::string& key, double fallback);

// Value of `key`; throws std::invalid_argument naming `owner` if absent.
double require(const Params& p, const std::string& key, const std::string& owner);

}  // namespace mg
