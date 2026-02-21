#ifndef GROWTH_MODULE_H
#define GROWTH_MODULE_H

#include <list>
#include <random>
#include "Config.h"
#include "SpeciesParams.h"

class AdultTree;
class Sapling;
class HerbLayer;
struct YearClimate;

class GrowthModule {
private:
    const SpeciesParamsManager& species_params_;
    
public:
    // ===========================================
    // Recruitment Variation Parameter
    // ===========================================
    // Random variation factor for recruitment (Â±mRecruitmentVariation)
    // value = base_value Ã— uniform_random(1 - mRV, 1 + mRV)
    static constexpr double mRecruitmentVariation = 0.1;  // Â±10% variation
    
    explicit GrowthModule(const SpeciesParamsManager& sp);
    
    // ===========================================
    // Environmental Response
    // ===========================================
    
    // Temperature response
    static double calcTempResponse(double GDD, double DD_min, double k_gdd);
    
    // Combined environmental factor (Liebig multiplicative, clamped to [0,1])
    // f_env = min(f_temp * f_light * f_drought, 1.0)
    static double calcEnvironmentalFactor(double f_temp, double f_light, double f_drought);
    
    // ===========================================
    // Adult Tree Growth
    // ===========================================
    
    // Calculate optimal DBH increment using modified Botkin equation
    // Î”D_opt = [g Ã— dbh Ã— (1 - dbhÃ—height/(dbh_maxÃ—h_max))] / 
    //          [2Ã—height + sÃ—dbhÃ—exp(-sÃ—dbh/(h_max-1.37))]
    static double calcOptimalDiameterIncrement(double dbh, double height, const SpeciesProfile& sp);
    
    // Calculate height from DBH using adult H-D relationship
    // H = adulta Ã— exp(adultb Ã— DBH) + 1.37 - adulta
    static double calcHeightFromDBH(double dbh, const SpeciesProfile& sp);
    
    // [Phase 5c v2] Adult tree growth with dual-curve H-D plasticity
    // H_plastic = H_low + (1-AL) × (H_high - H_low), with convergence lock & temp limit
    void growAdultTree(AdultTree& tree, double GDD, double DI_tree);
    void growAllAdultTrees(std::list<AdultTree>& trees, double GDD, double DI_tree);
    
    // ===========================================
    // Sapling Growth (linear d0-height relationship)
    // ===========================================
    
    void growSapling(Sapling& sapling, double GDD, double DI_sapling);
    void growAllSaplings(std::list<Sapling>& saplings, double GDD, double DI_sapling);
    
    // ===========================================
    // Herb Growth
    // ===========================================
    
    void growHerbLayer(HerbLayer& herb_layer, double GDD, double DI_herb);
    
    // ===========================================
    // Recruitment (sapling â†’ adult)
    // ===========================================
    
    // [Phase 3a v1] Convert sapling to adult — smooth transition
    // Height inherited from sapling (no jump), foliage inherited, stress inherited.
    // Only DBH gets micro-perturbation (±10%).
    AdultTree convertSaplingToAdult(const Sapling& sapling, int new_id, std::mt19937& rng) const;
    
    // ===========================================
    // Utility
    // ===========================================
    
    static double calcSpeciesGDD(const YearClimate& climate, bool is_evergreen, int year);
    
    // Generate random variation factor for recruitment
    static double getRandomVariation(std::mt19937& rng);
};

#endif // GROWTH_MODULE_H