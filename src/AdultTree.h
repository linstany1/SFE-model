#ifndef ADULT_TREE_H
#define ADULT_TREE_H

#include "Config.h"
#include <cmath>


// Forward declaration
struct SpeciesProfile;

class AdultTree {
public:
    // ===========================================
    // Identity & Position
    // ===========================================
    int id;                      // Global unique ID
    int species_id;              // Species ID
    double x, y;                 // Local coordinates [0, 30)
    
    // ===========================================
    // Structure & Morphology
    // ===========================================
    double dbh;                  // Diameter at breast height (cm)
    double height;               // Tree height (m)
    double crown_base_height;    // Crown base height (m)
    double crown_radius;         // Maximum crown radius (m)
    double leaf_area;            // Leaf area (m²)
    double fol_weight;           // Actual foliage biomass (kg) - persistent state variable
    double potential_fol_weight;  // Potential (allometric max) foliage biomass (kg)
    double leaf_area_density;    // Leaf area density in crown (m²/m³)
    
    // ===========================================
    // State Variables
    // ===========================================
    int age;                     // Tree age (years)
    int stress_years;            // Consecutive stress years counter
    
    // ===========================================
    // Environmental Response (Updated each year)
    // ===========================================
    double f_env;                // Combined environmental response factor
    double available_light;      // Light available at effective crown height
    
    // Constructor
    AdultTree();
    
    // ===========================================
    // Geometry Calculations
    // ===========================================
    
    // Update all geometric properties based on current DBH
    void updateGeometry(const SpeciesProfile& sp);
    
    // Update leaf area density in crown volume
    void updateLeafAreaDensity(const SpeciesProfile& sp);
    
    // Apply foliage loss when crown base height rises (self-pruning)
    // Uses volume-proportional retention: retention = (1 - f)^(2*cs + 1)
    // where f = (HB_new - HB_old) / (H - HB_old)
    void applyPruningLoss(const SpeciesProfile& sp, double HB_old);
    
    // Update crown base height with self-pruning logic
    void updateCrownBase(const SpeciesProfile& sp, double AL_top, double AL_HB);
    
    // ===========================================
    // Biomass Calculations (MODIFIED: removed calcBasalArea)
    // ===========================================
    
    // Calculate aboveground biomass (kg) using allometry
    double calcBiomass(const SpeciesProfile& sp) const;
    
    // NOTE: calcBasalArea() has been REMOVED per requirement 4B
    
    // ===========================================
    // Crown Geometry for Light Calculations
    // ===========================================
    
    // Calculate crown radius at height z
    double getCrownRadiusAtHeight(double z, double cs) const;
    
    // Calculate the height where a vertical ray at distance d enters the crown
    double getEntryHeight(double horizontal_dist, double cs) const;
    
    // Calculate effective crown height for light calculation
    double getEffectiveCrownHeight() const;
    
    // ===========================================
    // Seed Production
    // ===========================================
    
    // Calculate seed production for this tree
    double calcSeedProduction(const SpeciesProfile& sp) const;
    
    // ===========================================
    // State Management (OPTIMIZED: with recovery mechanism)
    // ===========================================
    
    // Update stress counter with recovery mechanism
    // If f_env < stress_threshold: stress_years += 1
    // If f_env >= stress_threshold (recovery):
    //   - If f_env >= GOOD_RECOVERY_THRESHOLD: stress_years -= 2 (fast recovery)
    //   - Else: stress_years -= 1 (slow recovery)
    void updateStressCounter(double current_f_env, double stress_threshold);
    
    // Increment age
    void incrementAge();
    
    // Check if tree should be checked for crown overlap
    bool canOverlapAt(double target_x, double target_y) const;
    
    // ===========================================
    // Crown Area for CC calculation (NEW)
    // ===========================================
    double calcCrownArea() const {
        return M_PI * crown_radius * crown_radius;
    }
};

#endif // ADULT_TREE_H