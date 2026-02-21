#ifndef CONFIG_H
#define CONFIG_H

// For MSVC compatibility: must define _USE_MATH_DEFINES before cmath
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include <cmath>

// Define M_PI if not available (MSVC does not define it by default)
#ifndef M_PI
#define M_PI 3.1415926
#endif

#include <cmath>

namespace Config {


    // ===========================================
    // Spatial Scales
    // ===========================================
    constexpr double PLOT_SIZE = 30.0;           // Plot edge length (m)
    constexpr double CELL_SIZE = 1.0;            // Micro-grid cell size (m)
    constexpr int GRID_DIM = 30;                 // Grid dimension (30x30 cells)
    constexpr double HASH_BUCKET_SIZE = 5.0;     // Spatial hash bucket size (m)
    constexpr int HASH_DIM = 6;                  // 30m / 5m = 6 buckets per dimension
    
    // ===========================================
    // Physical Constants
    // ===========================================
    constexpr double BREAST_HEIGHT = 1.37;       // DBH measurement height (m)
    constexpr double HERB_HEIGHT = 1.0;          // Fixed herb layer height (m)
    constexpr double INITIAL_SAPLING_HEIGHT = 0.05; // New seedling initial height (m)
    constexpr double MIN_SAPLING_DIST = 0.1;     // Minimum sapling spacing (m)
    constexpr double MIN_HERB_BIOMASS = 0.01;    // Minimum herb biomass (kg/mÂ²)
    
    // ===========================================
    // Seed Parameters
    // ===========================================
    constexpr double N_UNLIMITED = 100.0;        // Seed saturation threshold (seeds/mÂ²)
    constexpr int SEED_LIFESPAN = 3;             // Seed bank lifespan (years)
    constexpr double LAI_SAT = 3.0;              // Source strength saturation LAI
    constexpr double SEED_BACKGROUND_PROB = 0.001; // Long-distance background probability
    constexpr double SEED_NEAR_PROB = 0.95;      // Near-field saturation probability
    
    // ===========================================
    // Water Balance Parameters
    // ===========================================
    constexpr double WHC_TOP_DEFAULT = 180.0;    // MODIFIED: Top layer WHC (mm), ~30cm depth
    constexpr double C_ET = 130.0;               // Maximum water supply rate (mm/month)
    constexpr double INTERCEPTION_MIN = 0.06;    // Minimum interception coefficient
    constexpr double INTERCEPTION_MAX = 0.23;    // Maximum interception coefficient
    
    // NEW: Herb root distribution (fixed, not species-specific)
    constexpr double K_TOP_HERB = 1.0;           // Herb roots 100% in top layer
    constexpr double K_SUB_HERB = 0.0;           // Herb has no deep roots
    
    // ===========================================
    // Growing Degree Day Parameters
    // ===========================================
    constexpr double K_GDD_TREE = 0.0013;        // Tree temperature response coefficient
    constexpr double K_GDD_HERB = 0.0013;        // Herb temperature response coefficient
    constexpr double T_BASE = 5;               // GDD base temperature (Â°C)
    
    // ===========================================
    // Mortality Parameters
    // ===========================================
    constexpr int STRESS_THRESHOLD_YEARS = 3;    // Consecutive stress years for sapling death
    constexpr double DEFAULT_STRESS_THRESHOLD = 0.05; // Default f_env stress threshold (lowered from 0.1)
    constexpr double P_BASE_MORTALITY = 0.001;   // Base background mortality rate
    
    // NEW: Natural mortality parameters (Sigmoid curve)
    // P_nat = P_base + (1 - P_base) / (1 + exp(k * (inflection - Age_rel)))
    constexpr double MORTALITY_K = 20.0;         // Steepness of mortality curve
    constexpr double MORTALITY_INFLECTION = 0.95; // Age_rel at which P_nat reaches ~50%
    
    // NEW: Stress mortality parameters
    // P_stress = P_max_stress * (Intensity)^2 * DurationFactor
    // Intensity = max(0, 1 - f_env/stress_threshold)
    // DurationFactor = min(1.0, (stress_years/resilience_limit)^3)
    constexpr double P_MAX_STRESS = 0.33;        // Maximum stress mortality probability
    constexpr int RESILIENCE_LIMIT = 5;          // Years of stress tolerance (was 3)
    
    // NEW: Recovery mechanism parameters
    // When f_env >= stress_threshold, tree recovers
    // If f_env >= GOOD_RECOVERY_THRESHOLD: stress_years -= 2
    // Else: stress_years -= 1
    constexpr double GOOD_RECOVERY_THRESHOLD = 0.60; // Threshold for fast recovery
    
    // ===========================================
    // Simulation Control
    // ===========================================
    constexpr int SPIN_UP_YEARS = 500;           // Spin-up period duration
    constexpr int SATURATION_YEARS = 50;         // Saturated seed rain years
    constexpr int CLIMATE_BASE_WINDOW = 30;      // Base climate sampling window
    constexpr int DEFAULT_OUTPUT_INTERVAL = 5;   // Default output interval (years)
    
    // ===========================================
    // Light Calculation Parameters
    // ===========================================
    constexpr double DEFAULT_EXTINCTION_KE = 0.5; // Default extinction coefficient
    constexpr double MIN_LIGHT_FOR_LOG = 1e-6;   // Minimum light for log calculation
    constexpr double MAX_CROWN_RATIO = 0.9;      // Maximum crown base height ratio
    
    // ===========================================
    // Numerical Constants
    // ===========================================
    constexpr double EPSILON = 1e-9;             // Floating point comparison epsilon
    
    // ===========================================
    // [Phase 4] Toroidal Geometry Utilities
    // ===========================================
    // The plot uses periodic (toroidal) boundary conditions to eliminate
    // edge effects. Trees near one boundary interact with trees near the
    // opposite boundary as if the plot were infinitely tiled.
    
    // Toroidal 1D distance: shortest path on a ring of circumference PLOT_SIZE
    inline double toroidalDelta(double a, double b) {
        double d = a - b;
        double half = PLOT_SIZE * 0.5;
        if (d > half) d -= PLOT_SIZE;
        else if (d < -half) d += PLOT_SIZE;
        return d;
    }
    
    // Toroidal 2D Euclidean distance
    inline double toroidalDist(double x1, double y1, double x2, double y2) {
        double dx = toroidalDelta(x1, x2);
        double dy = toroidalDelta(y1, y2);
        return std::sqrt(dx * dx + dy * dy);
    }
    
    // ===========================================
    // Leap Year Support Functions
    // ===========================================
    inline constexpr bool isLeapYear(int year) {
        return (year % 4 == 0 && year % 100 != 0) || (year % 400 == 0);
    }
    
    // Non-leap year days by month (0-indexed placeholder, 1-12 actual months)
    constexpr int DAYS_NORMAL[13] = {
        0,   // Placeholder for 1-based indexing
        31,  // January
        28,  // February (non-leap)
        31,  // March
        30,  // April
        31,  // May
        30,  // June
        31,  // July
        31,  // August
        30,  // September
        31,  // October
        30,  // November
        31   // December
    };
    
    // Get days in month with leap year support
    inline constexpr int getDaysInMonth(int month, int year) {
        if (month == 2 && isLeapYear(year)) {
            return 29;
        }
        return DAYS_NORMAL[month];
    }

    // Legacy array for backward compatibility (non-leap year)
    constexpr int DAYS_IN_MONTH[13] = {
        0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31
    };
}

#endif // CONFIG_H