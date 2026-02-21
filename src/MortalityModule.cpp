#include "MortalityModule.h"
#include "AdultTree.h"
#include "Sapling.h"
#include "FloatUtil.h"
#include <cmath>
#include <algorithm>

MortalityModule::MortalityModule(const SpeciesParamsManager& sp) : species_params_(sp) {}

// OPTIMIZED: Age-dependent natural mortality using shifted Sigmoid
// P_nat = P_base + (1 - P_base) / (1 + exp(k * (inflection - Age_rel)))
// This ensures trees have very low mortality until ~95% of max lifespan
double MortalityModule::calcNaturalMortality(int age, int max_age) {
    if (max_age <= 0) return 1.0;  // Die immediately if no lifespan
    
    double Age_rel = static_cast<double>(age) / static_cast<double>(max_age);
    
    // Shifted Sigmoid: mortality stays low until near inflection point (0.95)
    // k=20 makes the curve steep, inflection=0.95 delays mortality to late life
    double exponent = Config::MORTALITY_K * (Config::MORTALITY_INFLECTION - Age_rel);
    double sigmoid = (1.0 - Config::P_BASE_MORTALITY) / (1.0 + std::exp(exponent));
    
    // Total natural mortality
    double P_nat = Config::P_BASE_MORTALITY + sigmoid;
    
    return FloatUtil::clamp(P_nat, 0.0, 1.0);
}

// OPTIMIZED: Stress mortality with non-linear response and resilience
// P_stress = P_max_stress * (Intensity)^2 * DurationFactor
// Intensity = max(0, 1 - f_env/stress_threshold)  [squared reduces mild stress impact]
// DurationFactor = min(1.0, (stress_years/resilience_limit)^3)  [cubic ramp-up]
// Only triggered when f_env < stress_threshold
double MortalityModule::calcStressMortality(double f_env, int stress_years, double stress_threshold) {
    if (f_env >= stress_threshold) {
        return 0.0;
    }
    
    // Stress intensity factor (how far below threshold)
    // Squared to reduce impact of mild stress
    double intensity = std::max(0.0, 1.0 - f_env / stress_threshold);
    double intensity_squared = intensity * intensity;
    
    // Duration factor using cubic function for gradual ramp-up
    // Reaches full effect at RESILIENCE_LIMIT years (default 5)
    double duration_ratio = static_cast<double>(stress_years) / 
                            static_cast<double>(Config::RESILIENCE_LIMIT);
    double duration_factor = std::min(1.0, duration_ratio * duration_ratio * duration_ratio);
    
    // Combined stress mortality probability
    double P_stress = Config::P_MAX_STRESS * intensity_squared * duration_factor;
    
    return FloatUtil::clamp(P_stress, 0.0, 1.0);
}

// Combined mortality using independent event formula:
// P_mort = 1 - (1 - P_nat) * (1 - P_stress)
double MortalityModule::calcTotalMortality(double P_nat, double P_stress) {
    double P_survive_both = (1.0 - P_nat) * (1.0 - P_stress);
    return FloatUtil::clamp(1.0 - P_survive_both, 0.0, 1.0);
}

bool MortalityModule::checkTreeMortality(const AdultTree& tree, std::mt19937& rng) const {
    const SpeciesProfile& sp = species_params_.getById(tree.species_id);
    
    // Natural mortality (age-dependent)
    double P_nat = calcNaturalMortality(tree.age, sp.max_age_yr);
    
    // Stress mortality (environment-dependent with smoothing)
    double P_stress = calcStressMortality(tree.f_env, tree.stress_years, sp.stress_threshold);
    
    // Total mortality probability
    double P_mort = calcTotalMortality(P_nat, P_stress);
    
    // Random draw to determine if tree dies
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    return dist(rng) < P_mort;
}

// [Phase 1c v1] Sapling mortality: background (intrinsic) + stress, paralleling adult logic.
//
// Two independent mortality channels combined via:
//   P_mort = 1 - (1 - P_base)(1 - P_stress_deterministic)
//
// Channel 1 — Background mortality (stochastic):
//   P_base = 4.0 / max_age_yr   (ForClim classic: ensures ~2% survive to max age)
//   e.g. Abge(450yr) → P_base ≈ 0.89%/yr, Lagr(400yr) → 1.0%/yr
//   Applied every year as a random draw.
//
// Channel 2 — Stress mortality (deterministic, unchanged):
//   Die when stress_years >= STRESS_THRESHOLD_YEARS (default 3).
//   With Phase 1a's decay fix, this now requires ~4 bad years in 5 to trigger,
//   rather than 3 strictly consecutive years.
//
// Zero new parameters: reuses max_age_yr (already in species_params.csv).
bool MortalityModule::checkSaplingMortality(const Sapling& sapling, double /*stress_threshold*/,
                                             int max_age_yr, std::mt19937& rng) const {
    // Channel 2: deterministic stress death (unchanged)
    if (sapling.stress_years >= Config::STRESS_THRESHOLD_YEARS) {
        return true;
    }
    
    // Channel 1: stochastic background mortality
    if (max_age_yr > 0) {
        double P_base = 4.0 / static_cast<double>(max_age_yr);
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        if (dist(rng) < P_base) {
            return true;
        }
    }
    
    return false;
}

std::vector<int> MortalityModule::getDeadTreeIds(const std::list<AdultTree>& trees, 
                                                  std::mt19937& rng) const {
    std::vector<int> dead_ids;
    for (const auto& tree : trees) {
        if (checkTreeMortality(tree, rng)) {
            dead_ids.push_back(tree.id);
        }
    }
    return dead_ids;
}

// [Phase 1c v1] Updated to pass max_age_yr and rng through to checkSaplingMortality.
std::vector<int> MortalityModule::getDeadSaplingIds(const std::list<Sapling>& saplings,
                                                     std::mt19937& rng) const {
    std::vector<int> dead_ids;
    for (const auto& sapling : saplings) {
        const SpeciesProfile& sp = species_params_.getById(sapling.species_id);
        if (checkSaplingMortality(sapling, sp.stress_threshold, sp.max_age_yr, rng)) {
            dead_ids.push_back(sapling.id);
        }
    }
    return dead_ids;
}

void MortalityModule::applyFireMortality(std::list<AdultTree>& trees, std::list<Sapling>& saplings) {
    trees.clear();
    saplings.clear();
}

double MortalityModule::calcExpectedMortalityRate(const std::list<AdultTree>& trees) const {
    if (trees.empty()) return 0;
    
    double total_prob = 0;
    for (const auto& tree : trees) {
        const SpeciesProfile& sp = species_params_.getById(tree.species_id);
        double P_nat = calcNaturalMortality(tree.age, sp.max_age_yr);
        double P_stress = calcStressMortality(tree.f_env, tree.stress_years, sp.stress_threshold);
        total_prob += calcTotalMortality(P_nat, P_stress);
    }
    
    return total_prob / trees.size();
}

int MortalityModule::countStressedTrees(const std::list<AdultTree>& trees, double stress_threshold) {
    int count = 0;
    for (const auto& tree : trees) {
        if (tree.f_env < stress_threshold) {
            count++;
        }
    }
    return count;
}

int MortalityModule::countStressedSaplings(const std::list<Sapling>& saplings) {
    int count = 0;
    for (const auto& sapling : saplings) {
        if (sapling.stress_years >= Config::STRESS_THRESHOLD_YEARS) {
            count++;
        }
    }
    return count;
}