// A very small test harness: CHECK records a failure and carries on, each test file registers one function.
#pragma once

#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>

struct Suite {
    int checks = 0;
    int failures = 0;
    void record(bool ok, const char* expr, const char* file, int line) {
        ++checks;
        if (!ok) {
            ++failures;
            std::printf("FAILED %s:%d  %s\n", file, line, expr);
        }
    }
};

#define CHECK(cond) suite.record((cond), #cond, __FILE__, __LINE__)
#define CHECK_NEAR(x, y, tol) suite.record(std::fabs((x) - (y)) <= (tol), #x " ~ " #y, __FILE__, __LINE__)
#define CHECK_THROWS(expr)                                 \
    do {                                                   \
        bool threw_ = false;                               \
        try {                                              \
            (void)(expr);                                  \
        } catch (const std::exception&) {                  \
            threw_ = true;                                 \
        }                                                  \
        suite.record(threw_, #expr " throws", __FILE__, __LINE__); \
    } while (0)

// True if the expression throws an exception whose message contains `needle`.
template <typename F>
bool throws_with(F&& f, const std::string& needle) {
    try {
        f();
    } catch (const std::exception& e) {
        return std::string(e.what()).find(needle) != std::string::npos;
    }
    return false;
}

void test_strutil(Suite& suite);
void test_material(Suite& suite);
void test_grid(Suite& suite);
void test_banded(Suite& suite);
void test_eigen(Suite& suite);
void test_wave_operator(Suite& suite);
void test_mode(Suite& suite);
void test_slab(Suite& suite);
void test_structure(Suite& suite);
void test_parser(Suite& suite);
void test_analysis(Suite& suite);
void test_expect(Suite& suite);
