#ifndef WATER_MODULE_H
#define WATER_MODULE_H

#include <list>
#include <map>
#include "Config.h"
#include "SpeciesParams.h"

class AdultTree;
class Sapling;
class HerbLayer;
struct YearClimate;
struct MonthlyClimate;

// NEW v2.0: Per-species water balance result
struct SpeciesWaterResult {
    double DI;              // Growing-season mean drought index
    double AT_total;        // Total actual transpiration (mm/yr)
    double AT_top;          // Transpiration from top layer (mm/yr)
    double AT_sub;          // Transpiration from sub layer (mm/yr)
    
    SpeciesWaterResult() : DI(0), AT_total(0), AT_top(0), AT_sub(0) {}
};

struct WaterBalanceResult {
    // NEW v2.0: Per-species results
    std::map<int, SpeciesWaterResult> by_species;
    
    // Summary values (backward-compatible)
    double DI_herb;
    double DI_annual_mean;     // Mean DI across all species (weighted)
    double AT_total;           // Total transpiration all species
    
    double soil_moisture_end;
    double water_top_mean;
    double water_sub_mean;
    
    WaterBalanceResult();
    
    // NEW v2.0: Get DI for a specific species (returns 0 if not found)
    double getSpeciesDI(int species_id) const;
};

class WaterModule {
private:
    const SpeciesParamsManager& species_params_;
    double WHC_total_;
    double WHC_top_;
    double WHC_sub_;
    
    // MODIFIED v2.0: Dual-layer independent water tracking
    double water_top_;      // Current water in top layer (mm)
    double water_sub_;      // Current water in sub layer (mm)
    
public:
    explicit WaterModule(const SpeciesParamsManager& sp);
    
    void initialize(double WHC_total);
    void reset();
    
    // Phenology factor (unchanged)
    static double getPhenologyFactor(bool is_evergreen, int month);
    
    // NEW v2.0: Compute ECRB (Effective Competing Root Biomass) per species per layer
    void calcECRB(const std::list<AdultTree>& trees,
                  const std::list<Sapling>& saplings,
                  const HerbLayer& herb_layer,
                  int month,
                  std::map<int, double>& ECRB_top,
                  std::map<int, double>& ECRB_sub) const;
    
    // Retained utility functions
    static double calcCrownCover(const std::list<AdultTree>& trees);
    static double calcTotalLAI(const std::list<AdultTree>& trees,
                               const std::list<Sapling>& saplings,
                               const HerbLayer& herb_layer,
                               const SpeciesParamsManager& sp_mgr);
    
    // MODIFIED v2.0: Monthly water balance with NPPC parallel competition
    void calcMonthlyWaterBalance(const MonthlyClimate& climate,
                                  const std::list<AdultTree>& trees,
                                  const std::list<Sapling>& saplings,
                                  const HerbLayer& herb_layer,
                                  std::map<int, double>& DI_month,
                                  std::map<int, double>& AT_month);
    
    // MODIFIED v2.0: Annual water balance with per-species results
    WaterBalanceResult calcAnnualWaterBalance(const YearClimate& climate,
                                               const std::list<AdultTree>& trees,
                                               const std::list<Sapling>& saplings,
                                               const HerbLayer& herb_layer,
                                               int year);
    
    // Drought response (unchanged)
    static double calcDroughtResponse(double DI, double dr_tol);
    
    // Getters
    double getWaterTop() const { return water_top_; }
    double getWaterSub() const { return water_sub_; }
    double getSoilMoisture() const { return water_top_ + water_sub_; }
    double getWHCTotal() const { return WHC_total_; }
    double getWHCTop() const { return WHC_top_; }
    double getWHCSub() const { return WHC_sub_; }
    void setSoilMoisture(double sm);
};

#endif // WATER_MODULE_H