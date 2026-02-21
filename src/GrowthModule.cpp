#include "GrowthModule.h"
#include "AdultTree.h"
#include "Sapling.h"
#include "HerbLayer.h"
#include "LightModule.h"
#include "WaterModule.h"
#include "ClimateManager.h"
#include "FloatUtil.h"
#include <cmath>
#include <algorithm>

GrowthModule::GrowthModule(const SpeciesParamsManager& sp) : species_params_(sp) {}

double GrowthModule::calcTempResponse(double GDD, double DD_min, double k_gdd) {
    return std::max(0.0, 1.0 - std::exp((DD_min - GDD) * k_gdd));
}

// [Phase 3c v1] Liebig multiplicative: f_env = f_temp * f_light * f_drought
// (Code uses direct product — the commented-out cbrt was never active.
//  Direct product matches ForClim's Liebig-style minimum resource limitation.)
double GrowthModule::calcEnvironmentalFactor(double f_temp, double f_light, double f_drought) {
    double product = f_temp * f_light * f_drought;
    if (product <= 0.0) return 0.0;
    return std::min(1.0, product);
}

// NEW: Optimal DBH increment using s and g parameters
// ÃŽâ€D_opt = [g Ãƒâ€” dbh Ãƒâ€” (1 - dbhÃƒâ€”height/(dbh_maxÃƒâ€”h_max))] / 
//          [2Ãƒâ€”height + sÃƒâ€”dbhÃƒâ€”exp(-sÃƒâ€”dbh/(h_max-1.37))]
double GrowthModule::calcOptimalDiameterIncrement(double dbh, double height,
                                                const SpeciesProfile& sp) {
    if (sp.dbh_max_cm <= 0 || sp.h_max_m <= Config::BREAST_HEIGHT) return 0;
    if (dbh <= 0) return 0;
    
    // Numerator: g Ãƒâ€” dbh Ãƒâ€” (1 - dbhÃƒâ€”height/(dbh_maxÃƒâ€”h_max))
    double size_factor = 1.0 - (dbh * height) / (sp.dbh_max_cm * sp.h_max_m);
    if (size_factor <= 0) return 0;  // Tree at maximum size
    
    double numerator = sp.g * dbh * size_factor;
    
    // Denominator: 2Ãƒâ€”height + sÃƒâ€”dbhÃƒâ€”exp(-sÃƒâ€”dbh/(h_max-1.37))
    double h_eff = sp.h_max_m - Config::BREAST_HEIGHT;
    if (h_eff <= 0) return 0;
    
    double exp_term = sp.s * dbh * std::exp(-sp.s * dbh / h_eff);
    double denominator = 2.0 * height + exp_term;
    
    if (std::abs(denominator) < Config::EPSILON) return 0;
    
    return std::max(0.0, numerator / denominator);
}

// Height from DBH using mean adult H-D curve (utility / fallback).
// [Phase 5c note] Main growth now uses dual-curve plasticity in growAdultTree.
double GrowthModule::calcHeightFromDBH(double dbh, const SpeciesProfile& sp) {
    // H = adulta * exp(adultb * DBH) + 1.37 - adulta
    // When DBHÃ¢â€ â€™0: H = adulta Ãƒâ€” 1 + 1.37 - adulta = 1.37 (correct boundary)
    double height = sp.adulta * std::exp(sp.adultb * dbh) + Config::BREAST_HEIGHT - sp.adulta;
    
    // Clamp to valid range
    return std::max(Config::BREAST_HEIGHT, std::min(height, sp.h_max_m));
}

// [Phase 5c v2] Dual-curve H-D plasticity + temperature height limit
//
// Replaces the single mean H-D curve with a light-driven envelope:
//   H_low(DBH)  = full-light, stocky form  (from field survey 95th percentile)
//   H_high(DBH) = deep-shade, slender form (from field survey 5th percentile)
//   H_plastic   = H_low + w × (H_high - H_low),  w = 1 - available_light
//
// Large-tree convergence lock: as DBH → dbh_max, height forced toward h_max_eff.
// Temperature height limit: h_max_eff = h_max × f_temp  (treeline dwarfing).
// Monotonic constraint: tree.height can only increase (trees don't shrink).
//
// Zero new species parameters — envelope params are field-survey inputs (5b).
void GrowthModule::growAdultTree(AdultTree& tree, double GDD, double DI_tree) {
    const SpeciesProfile& sp = species_params_.getById(tree.species_id);
    
    // --- Environmental factors (unchanged) ---
    double f_temp = calcTempResponse(GDD, sp.dd_min, Config::K_GDD_TREE);
    double f_light = LightModule::calcLightResponse(tree.available_light, sp.shade_tol);
    double f_drought = WaterModule::calcDroughtResponse(DI_tree, sp.dr_tol);
    double f_env = calcEnvironmentalFactor(f_temp, f_light, f_drought);
    
    // --- DBH growth (unchanged) ---
    double delta_dbh_opt = calcOptimalDiameterIncrement(tree.dbh, tree.height, sp);
    double delta_dbh = delta_dbh_opt * f_env;
    tree.dbh += delta_dbh;
    tree.dbh = std::max(0.0, std::min(tree.dbh, sp.dbh_max_cm));
    
    // --- [Phase 5c v2] Dual-curve H-D plasticity ---
    
    // Step A: Temperature-limited maximum height (treeline dwarfing)
    double h_max_eff = sp.h_max_m * f_temp;
    h_max_eff = std::max(Config::BREAST_HEIGHT, h_max_eff);
    
    // Step B: Dual-curve envelope interpolation
    double H_low  = sp.adulta_low  * std::exp(sp.adultb_low  * tree.dbh)
                    + Config::BREAST_HEIGHT - sp.adulta_low;
    double H_high = sp.adulta_high * std::exp(sp.adultb_high * tree.dbh)
                    + Config::BREAST_HEIGHT - sp.adulta_high;
    
    // Both curves capped by temperature-limited height
    H_low  = std::min(H_low,  h_max_eff);
    H_high = std::min(H_high, h_max_eff);
    
    // Light interpolation weight: dark → slender (w=1), bright → stocky (w=0)
    double w = 1.0 - tree.available_light;
    double H_plastic = H_low + w * (H_high - H_low);
    
    // Step C: Large-tree convergence lock
    // Sigmoid: lock ≈ 0 for DBH < 0.5×dbh_max, ≈ 1 for DBH > 0.9×dbh_max
    double dbh_rel = tree.dbh / sp.dbh_max_cm;
    double lock = 1.0 / (1.0 + std::exp(-10.0 * (dbh_rel - 0.7)));
    double H_target = H_plastic * (1.0 - lock) + h_max_eff * lock;
    
    // Step D: Monotonic constraint — trees never shrink
    tree.height = std::max(tree.height, std::min(H_target, h_max_eff));
    
    // --- Geometry update & stress ---
    tree.updateGeometry(sp);
    tree.updateStressCounter(f_env, sp.stress_threshold);
    tree.incrementAge();
}

void GrowthModule::growAllAdultTrees(std::list<AdultTree>& trees, 
                           double GDD, double DI_tree) {
    for (auto& tree : trees) {
        growAdultTree(tree, GDD, DI_tree);
    }
}

// NEW: Sapling growth using linear d0-height relationship
void GrowthModule::growSapling(Sapling& sapling, double GDD, double DI_sapling) {
    const SpeciesProfile& sp = species_params_.getById(sapling.species_id);
    
    double f_temp = calcTempResponse(GDD, sp.dd_min, Config::K_GDD_TREE);
    double f_light = LightModule::calcLightResponse(sapling.available_light, sp.shade_tol);
    double f_drought = WaterModule::calcDroughtResponse(DI_sapling, sp.dr_tol);
    
    // Liebig multiplicative f_env
    double f_env = calcEnvironmentalFactor(f_temp, f_light, f_drought);
    
    // Use Sapling::grow() which implements the linear d0-height relationship
    sapling.grow(sp, f_env);
    
    // Update stress counter
    sapling.updateStressCounter(sp.stress_threshold);
    
    sapling.incrementAge();
}

void GrowthModule::growAllSaplings(std::list<Sapling>& saplings, 
                         double GDD, double DI_sapling) {
    for (auto& sapling : saplings) {
        growSapling(sapling, GDD, DI_sapling);
    }
}

// Herb growth using Liebig multiplicative f_env
void GrowthModule::growHerbLayer(HerbLayer& herb_layer, double GDD, double DI_herb) {
    if (!species_params_.hasSpecies(0)) return;
    
    const SpeciesProfile& herb_sp = species_params_.getById(0);
    
    double f_temp = calcTempResponse(GDD, herb_sp.dd_min, Config::K_GDD_HERB);
    double f_drought = WaterModule::calcDroughtResponse(DI_herb, herb_sp.dr_tol);
    
    for (int u = 0; u < Config::GRID_DIM; ++u) {
        for (int v = 0; v < Config::GRID_DIM; ++v) {
            double AL = herb_layer.getTopLight(u, v);
            double f_light = HerbLayer::calcLightResponse(AL, herb_sp.shade_tol);
            
            // Liebig multiplicative f_env for herb layer
            double f_env = calcEnvironmentalFactor(f_temp, f_light, f_drought);
            
            herb_layer.updateCellBiomass(u, v, f_env, herb_sp);
        }
    }
}

// Generate random variation factor for recruitment
double GrowthModule::getRandomVariation(std::mt19937& rng) {
    std::uniform_real_distribution<double> dist(
        1.0 - mRecruitmentVariation, 
        1.0 + mRecruitmentVariation
    );
    return dist(rng);
}

// [Phase 3a v1] Convert sapling to adult — smooth transition, no jumps.
//
// OLD PROBLEMS (Bug 3.2):
//   (a) Height recalculated via adult H-D curve → jumps from 1.37m to 3.3m
//   (b) potential_fol_weight = 0 → triggers full allometric init in updateGeometry
//   (c) stress_years = 0 → accumulated stress memory wiped
//   (d) Independent random perturbations on H and DBH → H/D ratio corrupted
//
// NEW DESIGN:
//   (a) Height directly inherited from sapling (continuous at 1.37m)
//   (b) fol_weight inherited from sapling, potential_fol set to allometric value
//       → updateGeometry takes normal increment path, not full-init
//   (c) stress_years inherited from sapling (pressure continuity)
//   (d) Only DBH gets micro-perturbation; height is deterministic inheritance
//
// Zero new parameters.
AdultTree GrowthModule::convertSaplingToAdult(const Sapling& sapling, int new_id, 
                                              std::mt19937& rng) const {
    const SpeciesProfile& sp = species_params_.getById(sapling.species_id);
    
    AdultTree adult;
    adult.id = new_id;
    adult.species_id = sapling.species_id;
    adult.x = sapling.x;
    adult.y = sapling.y;
    adult.age = sapling.age;
    
    // (c) Inherit stress memory from sapling — no more amnesia at recruitment
    adult.stress_years = sapling.stress_years;
    
    // (d) DBH: only micro-perturbation on base value (±5%)
    double base_dbh = sapling.height / sp.sapHD;
    adult.dbh = base_dbh * getRandomVariation(rng);
    adult.dbh = std::max(0.05, std::min(adult.dbh, sp.dbh_max_cm));
    
    // (a) Height: directly inherited — NO recalculation from adult H-D curve
    // This eliminates the 1.37m → 3.3m jump.
    // Height will evolve via dual-curve H-D plasticity in subsequent annual growth cycles.
    adult.height = sapling.height;
    adult.height = std::max(Config::BREAST_HEIGHT, std::min(adult.height, sp.h_max_m));
    
    // Crown base height from open-grown baseline
    adult.crown_base_height = adult.height * (1.0-sp.alpha_base);
    
    // (b) Foliage: inherit from sapling, set potential to allometric ceiling
    // This ensures updateGeometry takes the NORMAL INCREMENT path (line 38-42)
    // instead of the FULL INIT path (line 34-37), preventing instant crown explosion.
    double sapling_fol = sapling.calcFoliageWeight(sp);
    adult.potential_fol_weight = sp.fw1 * std::pow(adult.dbh, sp.fw2);
    adult.fol_weight = std::min(sapling_fol, adult.potential_fol_weight);
    
    // Update derived geometry (crown_radius, leaf_area, leaf_area_density)
    // Note: updateGeometry will see potential_fol_weight > 0 and take normal path.
    // It will recalculate potential_fol_weight from current dbh (same value),
    // compute delta_W = max(0, new - old) ≈ 0, and keep fol_weight unchanged.
    adult.updateGeometry(sp);
    
    return adult;
}

double GrowthModule::calcSpeciesGDD(const YearClimate& climate, bool is_evergreen, int year) {
    return climate.calcGDD(is_evergreen, year);
}