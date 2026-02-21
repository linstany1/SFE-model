#ifndef MORTALITY_MODULE_H
#define MORTALITY_MODULE_H

#include <list>
#include <vector>
#include <random>
#include "Config.h"
#include "SpeciesParams.h"

class AdultTree;
class Sapling;

class MortalityModule {
private:
    const SpeciesParamsManager& species_params_;
    
public:
    explicit MortalityModule(const SpeciesParamsManager& sp);
    
    // OPTIMIZED: Age-dependent natural mortality using shifted Sigmoid
    // P_nat = P_base + (1 - P_base) / (1 + exp(k * (inflection - Age_rel)))
    // Ensures low mortality until ~95% of max lifespan
    static double calcNaturalMortality(int age, int max_age);
    
    // OPTIMIZED: Stress mortality with non-linear response and resilience
    // P_stress = P_max_stress * (Intensity)^2 * DurationFactor
    // Intensity = max(0, 1 - f_env/stress_threshold)  [squared for soft response]
    // DurationFactor = min(1.0, (stress_years/resilience_limit)^3)  [cubic ramp-up]
    static double calcStressMortality(double f_env, int stress_years, double stress_threshold);
    
    // Combined mortality: P_mort = 1 - (1 - P_nat) * (1 - P_stress)
    static double calcTotalMortality(double P_nat, double P_stress);
    
    bool checkTreeMortality(const AdultTree& tree, std::mt19937& rng) const;
    
    // [Phase 1b v1] Sapling mortality: background + stress, paralleling adult logic.
    // Added rng for stochastic background mortality draw.
    // Added max_age_yr for ForClim P_base = 4/Amax calculation.
    bool checkSaplingMortality(const Sapling& sapling, double stress_threshold,
                                int max_age_yr, std::mt19937& rng) const;
    
    std::vector<int> getDeadTreeIds(const std::list<AdultTree>& trees, std::mt19937& rng) const;
    
    // [Phase 1b v1] Updated to pass max_age_yr and rng through to checkSaplingMortality.
    std::vector<int> getDeadSaplingIds(const std::list<Sapling>& saplings,
                                        std::mt19937& rng) const;
    
    static void applyFireMortality(std::list<AdultTree>& trees, std::list<Sapling>& saplings);
    
    double calcExpectedMortalityRate(const std::list<AdultTree>& trees) const;
    static int countStressedTrees(const std::list<AdultTree>& trees, double stress_threshold);
    static int countStressedSaplings(const std::list<Sapling>& saplings);
};

#endif // MORTALITY_MODULE_H