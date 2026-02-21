#ifndef HERB_LAYER_H
#define HERB_LAYER_H

#include <array>
#include "Config.h"

// Forward declaration
struct SpeciesProfile;

class HerbLayer {
private:
    // NEW: Competition control flag
    // When false, herb layer grows but does not affect tree/sapling competition
    bool enabled_;

public:
    // Grid of biomass values (kg/m² per cell)
    std::array<std::array<double, Config::GRID_DIM>, Config::GRID_DIM> biomass;
    
    // Grid of top-of-herb available light values
    std::array<std::array<double, Config::GRID_DIM>, Config::GRID_DIM> al_top;
    
    // Constructor - initialize to minimum biomass
    HerbLayer();
    
    // Reset to minimum biomass (bare ground or post-fire)
    void reset();
    
    // ===========================================
    // Competition Control (NEW)
    // ===========================================
    
    // Set whether herb layer affects tree/sapling competition
    void setEnabled(bool enabled) { enabled_ = enabled; }
    
    // Check if herb competition is enabled
    bool isEnabled() const { return enabled_; }
    
    // ===========================================
    // LAI and Leaf Calculations
    // ===========================================
    
    // Calculate LAI at a specific cell
    double calcLAI(int u, int v) const;
    
    // Calculate foliage weight at a specific cell
    double calcFoliageWeight(int u, int v) const;
    
    // Get average LAI across all cells
    double getAverageLAI() const;
    
    // Get total biomass (kg)
    double getTotalBiomass() const;
    
    // ===========================================
    // Growth Calculations
    // ===========================================
    
    // Update biomass for a single cell using logistic growth
    void updateCellBiomass(int u, int v, double f_env, const SpeciesProfile& herb_sp);
    
    // Update all cells given spatially varying light conditions
    void updateAllCells(const SpeciesProfile& herb_sp, 
                        double f_temp, double f_drought,
                        const std::array<std::array<double, Config::GRID_DIM>, 
                                         Config::GRID_DIM>& light_grid);
    
    // ===========================================
    // Environmental Response Functions
    // ===========================================
    
    // Light response function (same as tree but using herb shade tolerance)
    static double calcLightResponse(double AL, int shade_tol);
    
    // ===========================================
    // Foliage Weight for Water Competition
    // ===========================================
    
    // Get total foliage weight for water demand calculation
    double getTotalFoliageWeight() const;
    
    // NEW: Get effective foliage weight (returns 0 if disabled)
    // Use this for water competition calculations
    double getEffectiveFoliageWeight() const;
    
    // ===========================================
    // Light Attenuation for Short Saplings
    // ===========================================
    
    // Calculate light available to a short sapling within the herb layer
    // NEW: Returns al_top directly if herb competition is disabled
    double calcLightForShortSapling(int u, int v, double sapling_height, double k_e) const;
    
    // ===========================================
    // Access Methods
    // ===========================================
    
    // Get biomass at cell
    double getBiomass(int u, int v) const;
    
    // Set biomass at cell
    void setBiomass(int u, int v, double value);
    
    // Set top light at cell
    void setTopLight(int u, int v, double value);
    
    // Get top light at cell
    double getTopLight(int u, int v) const;
    
    // Get mean floor light
    double getMeanFloorLight() const;
};

#endif // HERB_LAYER_H
