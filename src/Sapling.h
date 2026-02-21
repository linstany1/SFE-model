#ifndef SAPLING_H
#define SAPLING_H

#include <cmath>
#include "Config.h"

// Forward declaration
struct SpeciesProfile;

class Sapling {
public:
    // ===========================================
    // Identity & Position
    // ===========================================
    int id;                      // Global unique ID
    int species_id;              // Species ID
    double x, y;                 // Local coordinates [0, 30)
    
    // ===========================================
    // State Variables
    // ===========================================
    double d0;                   // NEW: Basal diameter (cm)
    double height;               // Current height (m)
    int age;                     // Sapling age (years)
    int stress_years;            // Consecutive stress years counter
    
    // ===========================================
    // Environmental Response (Updated each year)
    // ===========================================
    double f_env;                // Combined environmental response factor
    double available_light;      // Available light at sapling height
    
    // Constructor
    Sapling();
    
    // ===========================================
    // NEW: Initialization from Establishment
    // ===========================================
    
    // Initialize sapling with establishment conditions
    // d0 = deltaBD_pot × f_env × random_variation
    // height = ya × d0 × random_variation
    void initializeFromEstablishment(double f_env_val, const SpeciesProfile& sp,
                                      double d0_variation, double height_variation);
    
    // ===========================================
    // Leaf Area Calculations
    // ===========================================
    
    // Calculate leaf area using equivalent DBH from height-diameter ratio
    double calcLeafArea(const SpeciesProfile& sp) const;
    
    // Calculate foliage weight
    double calcFoliageWeight(const SpeciesProfile& sp) const;
    
    // ===========================================
    // Growth (NEW: Linear d0-height relationship)
    // ===========================================
    
    // Apply annual basal diameter and height growth
    // Δd0 = deltaBD_pot × f_env
    // height = ya × d0
    void grow(const SpeciesProfile& sp, double f_env_current);
    
    // ===========================================
    // Recruitment Check
    // ===========================================
    
    // Check if sapling should be recruited to adult
    bool checkRecruitment() const;
    
    // Calculate initial DBH upon recruitment (DBH = height / sapHD)
    double calcInitialDBH(const SpeciesProfile& sp) const;
    
    // ===========================================
    // Grid Location
    // ===========================================
    
    // Get grid cell indices
    int getGridU() const;
    int getGridV() const;
    
    // ===========================================
    // State Management
    // ===========================================
    
    // Update stress counter based on f_env and species stress_threshold
    void updateStressCounter(double stress_threshold);
    
    // Check if sapling should die from stress (3 consecutive years)
    bool checkStressDeath() const;
    
    // Increment age
    void incrementAge();
    
    // ===========================================
    // Light Type Classification
    // ===========================================
    
    // Check if this is a "tall" sapling (above herb layer)
    bool isTallSapling() const;
    
    // Check if this is a "short" sapling (within herb layer)
    bool isShortSapling() const;
};

#endif // SAPLING_H
