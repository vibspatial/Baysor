#pragma once

#include <vector>
#include <cstdint>
#include <random>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <string>
#include <unordered_map>

#include "baysor/utils/xoshiro.h"

namespace baysor {

inline constexpr double kPi = 3.141592653589793238462643383279502884;

// --- Fast math inlines (match Julia @fastmath fmax/fmin/fsort) ---
inline double fmax(double a, double b) { return a > b ? a : b; }
inline double fmin(double a, double b) { return a < b ? a : b; }
inline float fmax(float a, float b) { return a > b ? a : b; }
inline float fmin(float a, float b) { return a < b ? a : b; }
inline int imax(int a, int b) { return a > b ? a : b; }
inline int imin(int a, int b) { return a < b ? a : b; }

// --- Counting / histogramming ---

/// Count occurrences of each value in [1..max_value]. Index 0 of result = count of value 1.
/// If max_value <= 0, it's inferred from the data.
/// If drop_zero is true, zero-valued entries in `values` are silently skipped.
std::vector<int> count_array(const std::vector<int>& values, int max_value = -1, bool drop_zero = false);

/// Count with float weights
std::vector<double> count_array_weighted(const std::vector<int>& values,
                                         const std::vector<double>& weights,
                                         int max_value = -1, bool drop_zero = false);

// --- Group-by / split ---

/// Group indices 0..N-1 by factor[i], returning vector-of-vectors.
/// Result[k] = indices where factor[i] == k+1 (1-based factor values).
/// If drop_zero, factor values of 0 are skipped.
std::vector<std::vector<int>> split_ids(const std::vector<int>& factor,
                                         int max_factor = -1,
                                         bool drop_zero = false);

// --- Weighted sampling ---

/// Weighted categorical sample from weights[0..n-1], returns index
int fsample(const double* weights, int n, std::mt19937& rng);
int fsample(const double* weights, int n, Xoshiro256pp& rng);
Xoshiro256pp& global_xoshiro_rng();
void reset_global_xoshiro_rng(std::uint64_t seed = 1);

/// Sample from arr[weights], returns arr[sampled_index]
inline int fsample(const int* arr, const double* weights, int n, std::mt19937& rng) {
    return arr[fsample(weights, n, rng)];
}
inline int fsample(const int* arr, const double* weights, int n, Xoshiro256pp& rng) {
    return arr[fsample(weights, n, rng)];
}

// --- Weighted statistics ---
double wmean(const double* values, const double* weights, int n);
std::pair<double, double> wmean_std(const double* values, const double* weights, int n);

// --- Convergence ---
struct ConvergenceResult {
    double max_diff;
    double change_frac;
};

ConvergenceResult estimate_difference_l0(
    const double* m1, const double* m2, int rows, int cols,
    const double* col_weights = nullptr, double change_threshold = 1e-7
);

// --- String utilities ---
std::vector<std::string> split_string_list(const std::string& s, char sep = ',');
std::string strip(const std::string& s);

} // namespace baysor
