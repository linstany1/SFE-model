#include "WaterModule.h"
#include "AdultTree.h"
#include "Sapling.h"
#include "HerbLayer.h"
#include "ClimateManager.h"
#include "FloatUtil.h"
#include <cmath>
#include <algorithm>

// ============================================================
// WaterBalanceResult
// ============================================================

WaterBalanceResult::WaterBalanceResult() : 
    DI_herb(0), DI_annual_mean(0), AT_total(0),
    soil_moisture_end(0), water_top_mean(0), water_sub_mean(0) {}

double WaterBalanceResult::getSpeciesDI(int species_id) const {
    auto it = by_species.find(species_id);
    if (it != by_species.end()) {
        return it->second.DI;
    }
    return 0.0;
}

// ============================================================
// WaterModule Construction & Initialization
// ============================================================

WaterModule::WaterModule(const SpeciesParamsManager& sp) 
    : species_params_(sp), WHC_total_(0), WHC_top_(0), WHC_sub_(0),
      water_top_(0), water_sub_(0) {}

void WaterModule::initialize(double WHC_total) {
    WHC_total_ = WHC_total;
    WHC_top_ = std::min(Config::WHC_TOP_DEFAULT, WHC_total_);  // min(180, WHC_tot)
    WHC_sub_ = std::max(0.0, WHC_total_ - WHC_top_);
    water_top_ = WHC_top_;    // Start saturated
    water_sub_ = WHC_sub_;    // Start saturated
}

void WaterModule::reset() {
    water_top_ = WHC_top_;
    water_sub_ = WHC_sub_;
}

void WaterModule::setSoilMoisture(double sm) {
    sm = FloatUtil::clamp(sm, 0.0, WHC_total_);
    // Distribute from bottom up (gravity)
    water_sub_ = std::min(sm, WHC_sub_);
    water_top_ = std::min(sm - water_sub_, WHC_top_);
}

// ============================================================
// Phenology Factor (unchanged)
// ============================================================

double WaterModule::getPhenologyFactor(bool is_evergreen, int month) {
    if (is_evergreen) {
        return 1.0;
    } else {
        return (month >= 4 && month <= 10) ? 1.0 : 0.0;
    }
}

// ============================================================
// NEW v2.0: ECRB Computation
// ============================================================
// Computes Effective Competing Root Biomass per species per soil layer.
// ECRB_top_s = sum( foliage_weight_i * delta_leaf_i * k_top_s )
// ECRB_sub_s = sum( foliage_weight_i * delta_leaf_i * k_sub_s )
// where i iterates over individuals of species s.

void WaterModule::calcECRB(const std::list<AdultTree>& trees,
                           const std::list<Sapling>& saplings,
                           const HerbLayer& herb_layer,
                           int month,
                           std::map<int, double>& ECRB_top,
                           std::map<int, double>& ECRB_sub) const {
    ECRB_top.clear();
    ECRB_sub.clear();
    
    // --- Herb layer (species_id = 0) ---
    if (species_params_.hasSpecies(0)) {
        double delta_herb = getPhenologyFactor(false, month);
        double W_herb = herb_layer.getEffectiveFoliageWeight() * delta_herb;
        ECRB_top[0] = W_herb * Config::K_TOP_HERB;   // = W_herb * 1.0
        ECRB_sub[0] = W_herb * Config::K_SUB_HERB;   // = W_herb * 0.0
    }
    
    // --- Saplings (accumulated by species) ---
    for (const auto& sapling : saplings) {
        const SpeciesProfile& sp = species_params_.getById(sapling.species_id);
        double delta = getPhenologyFactor(sp.is_evergreen, month);
        double W = sapling.calcFoliageWeight(sp) * delta;
        double k_top = sp.k_top_sapling;
        ECRB_top[sapling.species_id] += W * k_top;
        ECRB_sub[sapling.species_id] += W * (1.0 - k_top);
    }
    
    // --- Adult trees (accumulated by species) ---
    for (const auto& tree : trees) {
        const SpeciesProfile& sp = species_params_.getById(tree.species_id);
        double delta = getPhenologyFactor(sp.is_evergreen, month);
        double W = tree.fol_weight * delta;
        double k_top = sp.k_top_adult;
        ECRB_top[tree.species_id] += W * k_top;
        ECRB_sub[tree.species_id] += W * (1.0 - k_top);
    }
}

// ============================================================
// Retained Utility Functions (unchanged)
// ============================================================

double WaterModule::calcCrownCover(const std::list<AdultTree>& trees) {
    double total_crown_area = 0;
    for (const auto& tree : trees) {
        total_crown_area += M_PI * tree.crown_radius * tree.crown_radius;
    }
    double plot_area = Config::PLOT_SIZE * Config::PLOT_SIZE;
    return std::min(1.0, total_crown_area / plot_area);
}

double WaterModule::calcTotalLAI(const std::list<AdultTree>& trees,
                               const std::list<Sapling>& saplings,
                               const HerbLayer& herb_layer,
                               const SpeciesParamsManager& sp_mgr) {
    double total_LA = 0;
    
    for (const auto& tree : trees) {
        total_LA += tree.leaf_area;
    }
    
    for (const auto& sapling : saplings) {
        total_LA += sapling.calcLeafArea(sp_mgr.getById(sapling.species_id));
    }
    
    if (herb_layer.isEnabled()) {
        total_LA += herb_layer.getAverageLAI() * Config::PLOT_SIZE * Config::PLOT_SIZE;
    }
    
    double plot_area = Config::PLOT_SIZE * Config::PLOT_SIZE;
    return total_LA / plot_area;
}

// ============================================================
// CORE v2.0: Monthly Water Balance — NPPC Parallel Competition
// ============================================================
//
// Algorithm overview:
//   1. Interception:    P_int, P_soil
//   2. Stand supply:    S_stand = min(C_ET, WHC_tot) * (water_top+water_sub)/WHC_tot
//   3. Stand demand:    D_stand, D_species (by foliage proportion)
//   4. Top supply:      S_top = min(water_top + P_soil, S_stand)
//   5. Phase 1:         Parallel competition for S_top via ECRB_top
//   6. Phase 2:         Deep compensation for unmet demand via ECRB_sub
//   7. End-of-month:    Store remaining water bottom-up; overflow leaves system
//   8. DI per species

void WaterModule::calcMonthlyWaterBalance(
        const MonthlyClimate& climate,
        const std::list<AdultTree>& trees,
        const std::list<Sapling>& saplings,
        const HerbLayer& herb_layer,
        std::map<int, double>& DI_month,
        std::map<int, double>& AT_month) {
    
    int month = climate.month;
    DI_month.clear();
    AT_month.clear();
    
    // ================================================================
    // Step 1: Interception (identical to original LandClim)
    // ================================================================
    double LAI = calcTotalLAI(trees, saplings, herb_layer, species_params_);
    double F_veg = (LAI <= 8.0) ? (0.75 + LAI / 32.0) : 1.0;
    
    double CC = calcCrownCover(trees);
    double alpha_int = Config::INTERCEPTION_MIN + 
                       (Config::INTERCEPTION_MAX - Config::INTERCEPTION_MIN) * 
                       std::min(1.0, CC);
    
    double P_int = std::min(alpha_int * climate.precipitation, F_veg * climate.pet);
    double P_soil = std::max(0.0, climate.precipitation - P_int);
    
    // ================================================================
    // Step 2: Stand-level supply ceiling (Federer 1982)
    // ================================================================
    // S_stand is the maximum amount the entire soil can release this month,
    // scaled by current relative moisture.
    double total_water = water_top_ + water_sub_;
    double S_stand = std::min(Config::C_ET, WHC_total_) * (total_water / WHC_total_);
    
    // ================================================================
    // Step 3: Stand-level demand & per-species demand
    // ================================================================
    double D_stand = std::max(0.0, F_veg * climate.pet - P_int);
    
    // Compute ECRB for this month
    std::map<int, double> ECRB_top, ECRB_sub;
    calcECRB(trees, saplings, herb_layer, month, ECRB_top, ECRB_sub);
    
    // Total effective foliage (W_s = ECRB_top_s + ECRB_sub_s for each species)
    double W_total = 0;
    std::map<int, double> W_species;
    for (auto& p : ECRB_top) {
        int sp_id = p.first;
        double W_sp = ECRB_top[sp_id] + ECRB_sub[sp_id];
        W_species[sp_id] = W_sp;
        W_total += W_sp;
    }
    // Ensure ECRB_sub keys also appear in W_species
    for (auto& p : ECRB_sub) {
        if (W_species.find(p.first) == W_species.end()) {
            W_species[p.first] = p.second;
            W_total += p.second;
        }
    }
    
    // Per-species demand = D_stand × (W_s / W_total)
    std::map<int, double> D_species;
    for (auto& p : W_species) {
        double frac = (W_total > Config::EPSILON) ? p.second / W_total : 0;
        D_species[p.first] = D_stand * frac;
    }
    
    // ================================================================
    // Step 4: Top-layer available supply
    // ================================================================
    // All P_soil enters the top pool as potential supply.
    // Supply is capped by S_stand (system physical limit).
    double S_top = std::min(water_top_ + P_soil, S_stand);
    
    // ================================================================
    // Step 5: Phase 1 — Top-layer parallel competition
    // ================================================================
    double ECRB_top_total = 0;
    for (auto& p : ECRB_top) ECRB_top_total += p.second;
    
    std::map<int, double> AT_top_sp;  // Per-species top-layer uptake
    double AT_phase1_total = 0;
    
    if (ECRB_top_total > Config::EPSILON && S_top > Config::EPSILON) {
        for (auto& p : ECRB_top) {
            int sp_id = p.first;
            double phi = p.second / ECRB_top_total;   // Competition share
            double supply_sp = phi * S_top;            // Supply available
            
            // Species top-layer demand = D_s × k_top_s
            // k_top_s can be recovered from ECRB ratio
            double k_top_sp = (W_species[sp_id] > Config::EPSILON) ?
                              ECRB_top[sp_id] / W_species[sp_id] : 0;
            double D_top_sp = D_species[sp_id] * k_top_sp;
            
            AT_top_sp[sp_id] = std::min(supply_sp, D_top_sp);
            AT_phase1_total += AT_top_sp[sp_id];
        }
    }
    
    // ================================================================
    // Step 6: Phase 2 — Deep-layer compensation competition
    // ================================================================
    // Remaining S_stand budget after Phase 1
    double S_stand_remaining = std::max(0.0, S_stand - AT_phase1_total);
    
    // Deep-layer supply: capped by actual deep water and remaining budget
    double S_sub = std::min(water_sub_, S_stand_remaining);
    
    // Only species with ECRB_sub > 0 participate
    double ECRB_sub_eligible = 0;
    std::map<int, double> D_unmet;
    
    for (auto& p : ECRB_sub) {
        int sp_id = p.first;
        if (p.second > Config::EPSILON) {
            double at_top = AT_top_sp.count(sp_id) ? AT_top_sp[sp_id] : 0;
            double unmet = std::max(0.0, D_species[sp_id] - at_top);
            if (unmet > Config::EPSILON) {
                D_unmet[sp_id] = unmet;
                ECRB_sub_eligible += p.second;
            }
        }
    }
    
    std::map<int, double> AT_sub_sp;
    double AT_phase2_total = 0;
    
    if (ECRB_sub_eligible > Config::EPSILON && S_sub > Config::EPSILON) {
        for (auto& p : D_unmet) {
            int sp_id = p.first;
            double psi = ECRB_sub[sp_id] / ECRB_sub_eligible;
            double supply_sub = psi * S_sub;
            AT_sub_sp[sp_id] = std::min(supply_sub, p.second);
            AT_phase2_total += AT_sub_sp[sp_id];
        }
    }
    
    // ================================================================
    // [Phase 5a v1] Step 6.5: Bare soil evaporation
    // ================================================================
    // When canopy cover is low (fire aftermath, open areas), the exposed
    // soil surface loses water via direct evaporation. This creates the
    // "nurse effect": pioneer species with high drought tolerance colonize
    // first, then shade the soil, reducing evaporation for late-successional
    // species that couldn't survive the bare-soil drought stress.
    //
    // canopy_cover_eff: Beer-Lambert effective cover from total LAI
    // E_soil: potential bare-soil evaporation (0.3 = bare soil PET fraction)
    // E_actual: limited by available top-layer water
    //
    // Dense forest (LAI=6): cover≈0.95 → E_soil≈0.015×PET (negligible)
    // Bare ground (LAI=0): cover=0    → E_soil=0.3×PET (significant)
    //
    // Zero new parameters: 0.3 is physical bare-soil evaporation coefficient.
    {
        double canopy_cover_eff = 1.0 - std::exp(-Config::DEFAULT_EXTINCTION_KE * LAI);
        double E_soil_potential = (1.0 - canopy_cover_eff) * climate.pet * 0.3;
        double E_actual = std::min(E_soil_potential, water_top_);
        water_top_ -= E_actual;
    }
    
    // ================================================================
    // Step 7: End-of-month water storage update (bottom-up fill)
    // ================================================================
    // Total actual transpiration
    double AT_all = AT_phase1_total + AT_phase2_total;
    
    // Remaining water = current storage (post-evaporation) + P_soil - transpiration
    // [Phase 5a v1] Use current water_top_ (reduced by E_actual) instead of stale total_water
    double remaining = (water_top_ + water_sub_) + P_soil - AT_all;
    remaining = std::max(0.0, remaining);
    
    // Overflow: if remaining exceeds WHC_tot, excess leaves the system
    if (remaining > WHC_total_) {
        remaining = WHC_total_;
    }
    
    // Fill from bottom up (gravity: water settles to deepest layer first)
    water_sub_ = std::min(remaining, WHC_sub_);
    water_top_ = std::min(remaining - water_sub_, WHC_top_);
    water_top_ = std::max(0.0, water_top_);
    
    // ================================================================
    // Step 8: Per-species DI calculation
    // ================================================================
    for (auto& p : D_species) {
        int sp_id = p.first;
        double D_sp = p.second;
        double AT_sp = (AT_top_sp.count(sp_id) ? AT_top_sp[sp_id] : 0)
                     + (AT_sub_sp.count(sp_id) ? AT_sub_sp[sp_id] : 0);
        
        AT_month[sp_id] = AT_sp;
        
        if (D_sp > Config::EPSILON) {
            DI_month[sp_id] = FloatUtil::clamp(1.0 - AT_sp / D_sp, 0.0, 1.0);
        } else {
            DI_month[sp_id] = 0.0;
        }
    }
}

// ============================================================
// MODIFIED v2.0: Annual Water Balance
// ============================================================

WaterBalanceResult WaterModule::calcAnnualWaterBalance(
        const YearClimate& climate,
        const std::list<AdultTree>& trees,
        const std::list<Sapling>& saplings,
        const HerbLayer& herb_layer,
        int /*year*/) {
    
    WaterBalanceResult result;
    
    // Accumulators per species
    std::map<int, double> sum_DI;
    std::map<int, double> sum_AT;
    std::map<int, double> sum_AT_top;
    std::map<int, double> sum_AT_sub;
    std::map<int, int>    grow_months;
    double sum_water_top = 0, sum_water_sub = 0;
    
    for (int m = 1; m <= 12; ++m) {
        std::map<int, double> DI_m, AT_m;
        
        calcMonthlyWaterBalance(climate[m], trees, saplings, herb_layer, DI_m, AT_m);
        
        // Accumulate AT
        for (auto& p : AT_m) {
            sum_AT[p.first] += p.second;
        }
        
        // Record layer water levels
        sum_water_top += water_top_;
        sum_water_sub += water_sub_;
        
        // Accumulate DI only in growing season (Apr-Oct, T >= 0)
        if (m >= 4 && m <= 10 && climate[m].temp_mean >= 0) {
            for (auto& p : DI_m) {
                sum_DI[p.first] += p.second;
                grow_months[p.first]++;
            }
        }
    }
    
    // Compute per-species annual results
    double total_AT = 0;
    double sum_all_DI = 0;
    int n_species = 0;
    
    for (auto& p : sum_AT) {
        int sp_id = p.first;
        SpeciesWaterResult swr;
        swr.AT_total = p.second;
        swr.DI = (grow_months.count(sp_id) && grow_months[sp_id] > 0) ? 
                 sum_DI[sp_id] / grow_months[sp_id] : 0;
        // AT_top and AT_sub not tracked monthly for simplicity; can be added if needed
        result.by_species[sp_id] = swr;
        total_AT += swr.AT_total;
        sum_all_DI += swr.DI;
        n_species++;
    }
    
    // Backward-compatible summary fields
    result.DI_herb = result.getSpeciesDI(0);
    result.DI_annual_mean = (n_species > 0) ? sum_all_DI / n_species : 0;
    result.AT_total = total_AT;
    result.soil_moisture_end = water_top_ + water_sub_;
    result.water_top_mean = sum_water_top / 12.0;
    result.water_sub_mean = sum_water_sub / 12.0;
    
    return result;
}

// ============================================================
// Drought Response (unchanged)
// ============================================================

double WaterModule::calcDroughtResponse(double DI, double dr_tol) {
    if (dr_tol <= Config::EPSILON) return 0.0;
    double ratio = 1.0 - DI / dr_tol;
    return std::sqrt(std::max(0.0, ratio));
}