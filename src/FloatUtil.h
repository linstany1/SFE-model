#ifndef FLOAT_UTIL_H
#define FLOAT_UTIL_H

#include <cmath>
#include "Config.h"

namespace FloatUtil {
    // Floating point equality comparison
    inline bool equal(double a, double b) {
        return std::abs(a - b) < Config::EPSILON;
    }
    
    // Less than comparison with epsilon tolerance
    inline bool less(double a, double b) {
        return a < b - Config::EPSILON;
    }
    
    // Less than or equal comparison
    inline bool lessOrEqual(double a, double b) {
        return a < b + Config::EPSILON;
    }
    
    // Greater than comparison with epsilon tolerance
    inline bool greater(double a, double b) {
        return a > b + Config::EPSILON;
    }
    
    // Greater than or equal comparison
    inline bool greaterOrEqual(double a, double b) {
        return a > b - Config::EPSILON;
    }
    
    // Check if value is approximately zero
    inline bool isZero(double a) {
        return std::abs(a) < Config::EPSILON;
    }
    
    // Safe division avoiding division by zero
    inline double safeDivide(double numerator, double denominator, double defaultValue = 0.0) {
        if (std::abs(denominator) < Config::EPSILON) {
            return defaultValue;
        }
        return numerator / denominator;
    }
    
    // Clamp value to [min, max] range
    template<typename T>
    inline T clamp(T value, T minVal, T maxVal) {
        if (value < minVal) return minVal;
        if (value > maxVal) return maxVal;
        return value;
    }
    
    // Safe logarithm avoiding log(0) or log(negative)
    inline double safeLog(double x) {
        if (x < Config::MIN_LIGHT_FOR_LOG) {
            return std::log(Config::MIN_LIGHT_FOR_LOG);
        }
        return std::log(x);
    }
    
    // Safe square root
    inline double safeSqrt(double x) {
        return std::sqrt(std::max(0.0, x));
    }
}

#endif // FLOAT_UTIL_H
