#ifndef ESTABLISHMENT_MODULE_H
#define ESTABLISHMENT_MODULE_H

#include <vector>
#include <list>
#include <set>
#include <map>
#include <random>
#include <string>
#include "Config.h"
#include "SpeciesParams.h"

class Sapling;
class SeedBank;
class HerbLayer;
struct WaterBalanceResult;  // NEW v2.0: forward declaration

struct EstablishmentDebugInfo {
    int species_id;
    double seeds_in_cell;    // Actual seed count in this cell
    double P_seed;           // Seed availability probability
    double P_env;            // Environmental suitability probability [0,1] continuous
    double P_est;            // Final establishment probability
    // [Phase 2a v1] Changed from bool to double: continuous response factors [0,1]
    double f_winter;         // Winter temperature response
    double f_temp;           // Growing degree-day response (reuses GrowthModule::calcTempResponse)
    double f_drought;        // Drought response (reuses WaterModule::calcDroughtResponse)
    double f_light;          // Light response (reuses LightModule::calcLightResponse)
    int n_established;       // Number established (0 or 1)
    std::string limiting_factor;  // Now records: "factor_name(value)" for min factor
    
    EstablishmentDebugInfo();
};

// ===========================================
// Seed Dispersal Statistics (NEW)
// ===========================================

// Per-species seed statistics for a plot
struct SpeciesSeedStats {
    double seeds_local;      // Seeds from internal + between-plot dispersal
    double seeds_external;   // Seeds from external sources
    double seeds_total;      // Total seeds (local + external)
    double seeds_mean_per_cell;  // Mean seeds per cell
    
    SpeciesSeedStats();
};

// Plot-level seed statistics
struct PlotSeedStats {
    std::map<int, SpeciesSeedStats> by_species;  // Per-species statistics
    double total_seeds_all_species;               // Total seeds across all species
    double total_seeds_local;                     // Total local seeds
    double total_seeds_external;                  // Total external seeds
    
    PlotSeedStats();
};

class EstablishmentModule {
private:
    const SpeciesParamsManager& species_params_;
    int next_sapling_id_;
    
public:
    explicit EstablishmentModule(const SpeciesParamsManager& sp);
    
    // [Phase 2a v1] Continuous environmental response functions replacing boolean filters.
    // Each returns [0, 1] probability instead of pass/fail boolean.
    
    // Winter temperature: smooth linear decay within WINTER_SMOOTH_WIDTH of boundaries
    static double calcWinterResponse(double T_winter, const SpeciesProfile& sp);
    
    // GDD: delegates to GrowthModule::calcTempResponse (already continuous)
    // Light: delegates to LightModule::calcLightResponse (already continuous)
    // Drought: delegates to WaterModule::calcDroughtResponse (already continuous)
    
    // P_env = f_winter × f_temp × f_light × f_drought (Liebig multiplicative)
    static double calcEnvironmentalProb(double f_winter, double f_temp,
                                         double f_drought, double f_light);
    
    bool generateSaplingCoordinates(int u, int v,
                                     const std::list<Sapling>& existing_saplings,
                                     std::vector<std::pair<double, double>>& new_coords_this_cell,
                                     std::mt19937& rng,
                                     double& out_x, double& out_y);
    
    // MODIFIED v2.0: Accept WaterBalanceResult for per-species DI
    void processEstablishment(const SeedBank& seed_bank,
                               const HerbLayer& herb_layer,
                               double GDD_deciduous,
                               double GDD_evergreen,
                               const WaterBalanceResult& water_result,
                               double T_winter,
                               const std::vector<int>& species_ids,
                               std::list<Sapling>& saplings,
                               std::mt19937& rng,
                               bool saturated_mode,
                               std::vector<EstablishmentDebugInfo>* debug_info = nullptr);
    
    int getNextSaplingId() const;
    void setNextSaplingId(int id);
    
    static void summarizeDebugInfo(const std::vector<EstablishmentDebugInfo>& info,
                                   std::map<int, int>& established_by_species,
                                   std::map<std::string, int>& limiting_factors);
};

#endif // ESTABLISHMENT_MODULE_H