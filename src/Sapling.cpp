#include "Sapling.h"
#include "SpeciesParams.h"
#include <algorithm>

// Constructor
Sapling::Sapling() :
    id(0), species_id(0), x(0), y(0),
    d0(0), height(Config::INITIAL_SAPLING_HEIGHT),
    age(0), stress_years(0), f_env(1.0), available_light(1.0) {}

// NEW: Initialize sapling from establishment with random variations
void Sapling::initializeFromEstablishment(double f_env_val, const SpeciesProfile& sp,
                                          double d0_variation, double height_variation) {
    // d0 = deltaBD_pot Ã— f_env Ã— random_variation
    d0 = sp.deltaBD_pot * f_env_val * d0_variation;
    d0 = std::max(0.01, d0);  // Minimum basal diameter
    
    // height = ya Ã— d0 Ã— random_variation
    height = sp.ya * d0 * height_variation;
    height = std::max(Config::INITIAL_SAPLING_HEIGHT, height);
    
    f_env = f_env_val;
    age = 0;
    stress_years = 0;
}

double Sapling::calcLeafArea(const SpeciesProfile& sp) const {
    // DBH_equiv = H / sapHD (equivalent DBH for leaf area calculation)
    double dbh_equiv = height / sp.sapHD;
    // FolW = fw1 * DBH^fw2
    double fol_weight = sp.fw1 * std::pow(dbh_equiv, sp.fw2);
    // LA = SLA * FolW
    return sp.SLA * fol_weight;
}

double Sapling::calcFoliageWeight(const SpeciesProfile& sp) const {
    double dbh_equiv = height / sp.sapHD;
    return sp.fw1 * std::pow(dbh_equiv, sp.fw2);
}

// NEW: Annual growth using linear d0-height relationship
void Sapling::grow(const SpeciesProfile& sp, double f_env_current) {
    f_env = f_env_current;
    
    // Basal diameter increment: Î”d0 = deltaBD_pot Ã— f_env
    double delta_d0 = sp.deltaBD_pot * f_env;
    d0 += delta_d0;
    
    // Height follows linear relationship: height = ya Ã— d0
    height = sp.ya * d0;
    
    // Apply constraints
    height = std::max(Config::INITIAL_SAPLING_HEIGHT, height);
    // Note: height can exceed 1.37m, recruitment check happens separately
}

bool Sapling::checkRecruitment() const {
    return height >= Config::BREAST_HEIGHT;
}

double Sapling::calcInitialDBH(const SpeciesProfile& sp) const {
    // DBH = height / sapHD
    return height / sp.sapHD;
}

int Sapling::getGridU() const {
    return static_cast<int>(x);
}

int Sapling::getGridV() const {
    return static_cast<int>(y);
}

void Sapling::updateStressCounter(double stress_threshold) {
    // [Phase 1a v1] Stress state update for saplings — gradual decay
    // Aligned with AdultTree::updateStressCounter decay logic.
    // Under stress:   N_{t+1} = N_t + 1
    // Recovery:        N_{t+1} = max(0, N_t - 1)   (was: N_t = 0)
    // Rationale: Saplings have limited carbon reserves, so recovery is slower
    // than adults (who get -2 under excellent conditions). A single good year
    // no longer grants full immunity — it merely reduces accumulated stress by 1.
    if (f_env < stress_threshold) {
        stress_years++;
    } else {
        stress_years = std::max(0, stress_years - 1);
    }
}

bool Sapling::checkStressDeath() const {
    // Sapling dies if 3 consecutive years f_env < stress_threshold
    return stress_years >= Config::STRESS_THRESHOLD_YEARS;
}

void Sapling::incrementAge() {
    age++;
}

bool Sapling::isTallSapling() const {
    return height >= Config::HERB_HEIGHT;
}

bool Sapling::isShortSapling() const {
    return height < Config::HERB_HEIGHT;
}