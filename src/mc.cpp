#include "mc.hpp"

#include <iostream>

const int PIECEWISE_LINEAR_DIST_INTERVALS = 1000;
const double UCN_E_MIN = 1.;
const double UCN_E_MAX = 215.;

template<typename UnaryFunction>
std::piecewise_linear_distribution<double> parse_distribution(UnaryFunction f, const double range_min, const double range_max){
	if (range_min > range_max)
		throw std::runtime_error("Invalid range used to parse distribution");
	double rmin = range_min;
	double rmax = range_max;
	int nw = PIECEWISE_LINEAR_DIST_INTERVALS;
	if (range_min == range_max){
		rmax = std::nextafter(range_max, std::numeric_limits<double>::max()); // use slightly larger double value for range_max, else distribution will default to range 0..1
		nw = 1;
	}

	return std::piecewise_linear_distribution<double>(nw, rmin, rmax, f);
}

// v^2 distribution
double v2Spectrum(const double v){
	return v*v;
}

// v^3 distribution
double v3Spectrum(const double v){
	return v*v*v;
}

std::piecewise_linear_distribution<double> v2_distribution = parse_distribution(v2Spectrum, UCN_E_MIN, UCN_E_MAX);

std::piecewise_linear_distribution<double> v3_distribution = parse_distribution(v3Spectrum, UCN_E_MIN, UCN_E_MAX);

std::uniform_int_distribution<int> zeroOrOne(0,1);

std::uniform_real_distribution<double> zeroToOne(0, std::nextafter(1, std::numeric_limits<double>::max()));
