#include "EstablishmentModule.h"
#include "Sapling.h"
#include "SeedBank.h"
#include "HerbLayer.h"
#include "GrowthModule.h"
#include "LightModule.h"
#include "WaterModule.h"
#include <cmath>
#include <algorithm>
#include <array>
#include <cstdio>

// [Phase 2b v1] Updated: continuous factors initialized to 0.0 instead of false
EstablishmentDebugInfo::EstablishmentDebugInfo() : species_id(0), seeds_in_cell(0),
    P_seed(0), P_env(0), P_est(0),
    f_winter(0), f_temp(0), f_drought(0), f_light(0), 
    n_established(0), limiting_factor("") {}

SpeciesSeedStats::SpeciesSeedStats() : seeds_local(0), seeds_external(0), 
    seeds_total(0), seeds_mean_per_cell(0) {}

PlotSeedStats::PlotSeedStats() : total_seeds_all_species(0), 
    total_seeds_local(0), total_seeds_external(0) {}

EstablishmentModule::EstablishmentModule(const SpeciesParamsManager& sp) 
    : species_params_(sp), next_sapling_id_(1) {}

// [Phase 2b v1] Continuous winter temperature response.
// Smooth linear decay within buffer zone at species range boundaries.
// Outside range: 0. Well inside range: 1.0.  Near boundary: linear ramp.
//
// Buffer width = 3°C (physical constant: typical interannual temperature variability).
// NOT a species parameter — represents measurement/climate uncertainty at range edges.
//
// Derivation: At T = T_min → f = 0 (same as old boolean).
//             At T = T_min + 3°C → f = 1.0 (full suitability).
//             At T between: linear interpolation.
double EstablishmentModule::calcWinterResponse(double T_winter, const SpeciesProfile& sp) {
    constexpr double WINTER_SMOOTH_WIDTH = 3.0;  // °C
    
    double f_low  = std::min(1.0, std::max(0.0,
                   (T_winter - sp.ngs_tmin_c) / WINTER_SMOOTH_WIDTH));
    double f_high = std::min(1.0, std::max(0.0,
                   (sp.ngs_tmax_c - T_winter) / WINTER_SMOOTH_WIDTH));
    return std::min(f_low, f_high);
}

// [Phase 2b v1] P_env = product of four continuous response factors.
// Replaces the old boolean all-or-nothing gate.
// Each factor ∈ [0, 1]; product ∈ [0, 1].
double EstablishmentModule::calcEnvironmentalProb(double f_winter, double f_temp,
                                                   double f_drought, double f_light) {
    return f_winter * f_temp * f_drought * f_light;
}

bool EstablishmentModule::generateSaplingCoordinates(int u, int v,
                                     const std::list<Sapling>& existing_saplings,
                                     std::vector<std::pair<double, double>>& new_coords_this_cell,
                                     std::mt19937& rng,
                                     double& out_x, double& out_y) {
    std::uniform_real_distribution<double> pos_dist(0.0, Config::CELL_SIZE);
    
    double base_x = u * Config::CELL_SIZE;
    double base_y = v * Config::CELL_SIZE;
    
    int max_attempts = 10;
    for (int attempt = 0; attempt < max_attempts; ++attempt) {
        double x = base_x + pos_dist(rng);
        double y = base_y + pos_dist(rng);
        
        x = std::max(0.0, std::min(x, Config::PLOT_SIZE - Config::EPSILON));
        y = std::max(0.0, std::min(y, Config::PLOT_SIZE - Config::EPSILON));
        
        bool too_close = false;
        
        for (const auto& sapling : existing_saplings) {
            double dx = x - sapling.x;
            double dy = y - sapling.y;
            if (std::sqrt(dx*dx + dy*dy) < Config::MIN_SAPLING_DIST) {
                too_close = true;
                break;
            }
        }
        
        if (!too_close) {
            for (const auto& coord : new_coords_this_cell) {
                double dx = x - coord.first;
                double dy = y - coord.second;
                if (std::sqrt(dx*dx + dy*dy) < Config::MIN_SAPLING_DIST) {
                    too_close = true;
                    break;
                }
            }
        }
        
        if (!too_close) {
            out_x = x;
            out_y = y;
            return true;
        }
    }
    
    return false;
}

// MODIFIED v2.0: Accept WaterBalanceResult for per-species DI lookup
void EstablishmentModule::processEstablishment(const SeedBank& seed_bank,
                               const HerbLayer& herb_layer,
                               double GDD_deciduous,
                               double GDD_evergreen,
                               const WaterBalanceResult& water_result,
                               double T_winter,
                               const std::vector<int>& species_ids,
                               std::list<Sapling>& saplings,
                               std::mt19937& rng,
                               bool saturated_mode,
                               std::vector<EstablishmentDebugInfo>* debug_info) {
    
    std::uniform_real_distribution<double> prob_dist(0.0, 1.0);
    
    std::array<std::array<std::vector<std::pair<double, double>>, Config::GRID_DIM>, Config::GRID_DIM> new_coords;
    
    std::array<std::array<std::set<int>, Config::GRID_DIM>, Config::GRID_DIM> species_in_cell;
    
    // [Phase 4d v2] Cell capacity grid: count all saplings per 1m² cell (cross-species).
    // Physical constant: 1m² cannot support more than 3 sapling stems regardless of species.
    // This prevents density blow-up and provides implicit self-thinning at establishment.
    constexpr int MAX_SAPLINGS_PER_CELL = 3;
    std::array<std::array<int, Config::GRID_DIM>, Config::GRID_DIM> sapling_count{};
    
    for (const auto& sapling : saplings) {
        int u = sapling.getGridU();
        int v = sapling.getGridV();
        if (u >= 0 && u < Config::GRID_DIM && v >= 0 && v < Config::GRID_DIM) {
            species_in_cell[u][v].insert(sapling.species_id);
            sapling_count[u][v]++;
        }
    }
    
    for (int u = 0; u < Config::GRID_DIM; ++u) {
        for (int v = 0; v < Config::GRID_DIM; ++v) {
            // [Phase 4d v2] Cell capacity check — skip saturated cells
            if (sapling_count[u][v] >= MAX_SAPLINGS_PER_CELL) continue;
            const SpeciesProfile& herb_sp = species_params_.hasSpecies(0) ? 
                species_params_.getById(0) : species_params_.getByIndex(0);
            double AL_est = herb_layer.calcLightForShortSapling(
                u, v, Config::INITIAL_SAPLING_HEIGHT, herb_sp.extinction_ke);
            
            for (int species_id : species_ids) {
                if (!species_params_.hasSpecies(species_id)) continue;
                if (species_id == 0) continue;
                
                const SpeciesProfile& sp = species_params_.getById(species_id);
                
                if (species_in_cell[u][v].count(species_id) > 0) continue;
                
                // Get actual seed count from seed bank
                double N_total = seed_bank.getSeeds(u, v, species_id);
                
                double P_seed;
                if (saturated_mode) {
                    P_seed = 1.0;
                } else {
                    P_seed = std::min(N_total / Config::N_UNLIMITED, 1.0);
                }
                
                if (P_seed < Config::EPSILON) continue;
                
                double GDD = sp.is_evergreen ? GDD_evergreen : GDD_deciduous;
                
                // MODIFIED v2.0: Use per-species DI
                // If species has no DI in results (new arrival), fall back to mean DI
                double DI_sp = water_result.getSpeciesDI(species_id);
                if (water_result.by_species.find(species_id) == water_result.by_species.end()) {
                    // Species not yet present on plot Ã¢â‚¬â€ use overall mean as fallback
                    DI_sp = water_result.DI_annual_mean;
                }
                
                // [Phase 2b v1] Continuous environmental response factors
                // Each factor ∈ [0, 1], reusing existing response functions from
                // GrowthModule, LightModule, WaterModule.
                double f_winter  = calcWinterResponse(T_winter, sp);
                double f_temp    = GrowthModule::calcTempResponse(GDD, sp.dd_min,
                                                                   Config::K_GDD_TREE);
                double f_light   = LightModule::calcLightResponse(AL_est, sp.shade_tol);
                double f_drought = WaterModule::calcDroughtResponse(DI_sp, sp.dr_tol);
                
                double P_env = calcEnvironmentalProb(f_winter, f_temp, f_drought, f_light);
                
                double P_est = P_seed * P_env;
                
                if (debug_info) {
                    EstablishmentDebugInfo info;
                    info.species_id = species_id;
                    info.seeds_in_cell = N_total;
                    info.P_seed = P_seed;
                    info.P_env = P_env;
                    info.P_est = P_est;
                    info.f_winter = f_winter;
                    info.f_temp = f_temp;
                    info.f_drought = f_drought;
                    info.f_light = f_light;
                    
                    // [Phase 2b v1] Record the minimum factor as the limiting one.
                    // Format: "FactorName(value)" for easy diagnosis.
                    double min_val = f_winter;
                    std::string min_name = "Winter";
                    if (f_temp < min_val)    { min_val = f_temp;    min_name = "GDD"; }
                    if (f_drought < min_val) { min_val = f_drought; min_name = "Drought"; }
                    if (f_light < min_val)   { min_val = f_light;   min_name = "Light"; }
                    
                    if (P_env > 1.0 - Config::EPSILON) {
                        info.limiting_factor = "None";
                    } else {
                        // Format: "Light(0.342)" — identifies bottleneck at a glance
                        char buf[64];
                        std::snprintf(buf, sizeof(buf), "%s(%.3f)", min_name.c_str(), min_val);
                        info.limiting_factor = buf;
                    }
                    
                    debug_info->push_back(info);
                }
                
                if (prob_dist(rng) < P_est) {
                    // [Phase 4d v2] Re-check capacity before each establishment
                    if (sapling_count[u][v] >= MAX_SAPLINGS_PER_CELL) continue;
                    
                    double x, y;
                    if (generateSaplingCoordinates(u, v, saplings, 
                                                    new_coords[u][v], rng, x, y)) {
                        Sapling new_sapling;
                        new_sapling.id = next_sapling_id_++;
                        new_sapling.species_id = species_id;
                        new_sapling.x = x;
                        new_sapling.y = y;
                        new_sapling.age = 0;
                        new_sapling.stress_years = 0;
                        
                        // Calculate f_env for initialization using per-species DI
                        double f_temp = GrowthModule::calcTempResponse(GDD, sp.dd_min, Config::K_GDD_TREE);
                        double f_light = LightModule::calcLightResponse(AL_est, sp.shade_tol);
                        double f_drought = WaterModule::calcDroughtResponse(DI_sp, sp.dr_tol);
                        double f_env_init = GrowthModule::calcEnvironmentalFactor(f_temp, f_light, f_drought);
                        
                        // Generate independent random variations
                        double d0_variation = GrowthModule::getRandomVariation(rng);
                        double height_variation = GrowthModule::getRandomVariation(rng);
                        
                        // Initialize sapling with d0 and height using random variations
                        new_sapling.initializeFromEstablishment(f_env_init, sp, d0_variation, height_variation);
                        
                        saplings.push_back(new_sapling);
                        new_coords[u][v].push_back({x, y});
                        species_in_cell[u][v].insert(species_id);
                        sapling_count[u][v]++;  // [Phase 4d v2] track new establishment
                        
                        if (debug_info && !debug_info->empty()) {
                            debug_info->back().n_established = 1;
                        }
                    }
                }
            }
        }
    }
}

int EstablishmentModule::getNextSaplingId() const { return next_sapling_id_; }
void EstablishmentModule::setNextSaplingId(int id) { next_sapling_id_ = id; }

void EstablishmentModule::summarizeDebugInfo(const std::vector<EstablishmentDebugInfo>& info,
                               std::map<int, int>& established_by_species,
                               std::map<std::string, int>& limiting_factors) {
    for (const auto& entry : info) {
        if (entry.n_established > 0) {
            established_by_species[entry.species_id] += entry.n_established;
        }
        if (!entry.limiting_factor.empty() && entry.limiting_factor != "None") {
            limiting_factors[entry.limiting_factor]++;
        }
    }
}