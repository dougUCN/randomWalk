#ifndef MC_H_
#define MC_H_

#include <random>

typedef std::mt19937_64 TMCGenerator; ///< typedef to default random-number generator

/**
 * Create a piecewise linear distribution from a function with single parameter
 *
 * @param f Function returning a double using a single double parameter
 * @param range_min Lower range limit of distribution
 * @param range_max Upper range limit of distribution
 *
 * @return Return piecewise linear distribution
 */
template<typename UnaryFunction>
std::piecewise_linear_distribution<double> parse_distribution(UnaryFunction f, const double range_min, const double range_max);

// v^2 distribution
double v2Spectrum(const double v);

// v^3 distribution
double v3Spectrum(const double v);

extern std::piecewise_linear_distribution<double> v2_distribution;

extern std::piecewise_linear_distribution<double> v3_distribution;

extern std::uniform_int_distribution<int> zeroOrOne;

extern std::uniform_real_distribution<double> zeroToOne;

#endif /*MC_H_*/
