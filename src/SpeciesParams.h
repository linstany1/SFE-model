#ifndef SPECIES_PARAMS_H
#define SPECIES_PARAMS_H

#include <string>
#include <vector>
#include <map>

// Structure holding all species-specific parameters
struct SpeciesProfile {
    // ===========================================
    // Basic Information
    // ===========================================
    int species_id;              // Unique index (0=Herb, 1..N=Trees)
    std::string species_name;    // Species name for logging
    bool is_evergreen;           // 1=evergreen, 0=deciduous
    
    // ===========================================
    // Dispersal Parameters (Tree only)
    // ===========================================
    double lambda1_m;            // Effective dispersal distance (m)
    double lambda2_m;            // Maximum dispersal distance (m)
    double fecundity_f;          // Seed production (seeds/mÂ² leaf area)
    int maturity_age_yr;         // Minimum age for seed production
    
    // ===========================================
    // Growth Parameters
    // ===========================================
    int max_age_yr;              // Maximum physiological lifespan
    double h_max_m;              // Maximum height (m)
    double dbh_max_cm;           // Maximum DBH (cm)
    double r_a;                  // Crown radius allometry parameter A (slope)
    double r_b;                  // Crown radius allometry parameter B (intercept)
    double cs;                   // Crown shape coefficient (1.0=cone)
    double alpha_base;           // Open-grown crown base height ratio
    double SLA;                  // Specific leaf area (mÂ²/kg)
    double s;                    // Adult H-D shape parameter (renamed from opt_s)
    double g;                    // Adult growth rate parameter (renamed from ga)
    double gs;                   // [DEPRECATED] Kept for CSV compatibility, not used
    double fw1;                  // Leaf biomass allometry parameter 1
    double fw2;                  // Leaf biomass allometry parameter 2
    double sapHD;                // Sapling height-diameter ratio (m/cm)
    
    // ===========================================
    // NEW: Sapling Growth Parameters
    // ===========================================
    double deltaBD_pot;          // Sapling max potential basal diameter increment (cm/yr)
    double ya;                   // Sapling height-d0 relationship slope (m/cm)
    
    // ===========================================
    // NEW: Adult H-D Relationship Parameters
    // ===========================================
    double adulta;               // Adult H-D relationship parameter a (m)
    double adultb;               // Adult H-D relationship parameter b (1/cm)
    
    // ===========================================
    // [Phase 5b v1] H-D Plasticity Envelope Parameters
    // ===========================================
    // From field survey quantile fits:
    //   _low  = full-light stocky form (95th percentile H/D)
    //   _high = deep-shade slender form (5th percentile H/D)
    // H(DBH) = adulta_x * exp(adultb_x * DBH) + 1.37 - adulta_x
    double adulta_low;           // H-D lower envelope a (full-light, stocky)
    double adultb_low;           // H-D lower envelope b
    double adulta_high;          // H-D upper envelope a (deep-shade, slender)
    double adultb_high;          // H-D upper envelope b
    
    // ===========================================
    // Biomass Allometry
    // ===========================================
    double B0, B1, B2;           // Large tree biomass params (DBH >= 5cm)
    double b0_small, b1_small, b2_small; // Small tree biomass params (DBH < 5cm)
    
    // ===========================================
    // Environmental Tolerance
    // ===========================================
    int shade_tol;               // Shade tolerance (1=intolerant .. 5=tolerant)
    double dr_tol;               // Drought tolerance threshold
    double dd_min;               // Minimum growing degree days requirement
    double l_min;                // Light compensation point (0-1)
    double l_max;                // Maximum light tolerance (0-1)
    double ngs_tmin_c;           // Non-growing season minimum temp tolerance (Â°C)
    double ngs_tmax_c;           // Non-growing season maximum temp tolerance (Â°C)
    
    // ===========================================
    // Mortality Parameters (NEW - Section 5)
    // ===========================================
    double stress_threshold;     // f_env threshold for stress (default 0.1)
    
    // ===========================================
    // Herb-specific Parameters
    // ===========================================
    double k_max_herb;           // Maximum herb biomass (kg/mÂ²)
    double r_max_herb;           // Herb intrinsic growth rate
    
    // ===========================================
    // Common Parameters
    // ===========================================
    double extinction_ke;        // Extinction coefficient (default 0.5)
    
    // ===========================================
    // NEW v2.0: Root Distribution Parameters
    // ===========================================
    // Fraction of root biomass allocated to the top soil layer.
    // k_sub = 1.0 - k_top (computed automatically, not stored).
    // For herb (species_id=0): k_top = 1.0 (fixed in Config).
    // Values constrained by stable isotope observations (see tech doc Â§3.4.1).
    double k_top_sapling;        // Sapling-stage top-layer root fraction (e.g. 0.58 for Abies)
    double k_top_adult;          // Adult-stage top-layer root fraction  (e.g. 0.35 for Abies)
    
    // Constructor with default values
    SpeciesProfile();
    
    // Check if this is a tree species (not herb)
    bool isTree() const { return species_id > 0; }
    
    // Check if this is the herb species
    bool isHerb() const { return species_id == 0; }
    
    // NEW v2.0: Get k_top for a given life stage
    double getKtop(bool is_adult) const {
        if (isHerb()) return 1.0;  // Herb always 100% top layer
        return is_adult ? k_top_adult : k_top_sapling;
    }
    
    // NEW v2.0: Get k_sub for a given life stage
    double getKsub(bool is_adult) const {
        return 1.0 - getKtop(is_adult);
    }
};

// Class to manage species parameters
class SpeciesParamsManager {
private:
    std::vector<SpeciesProfile> profiles_;
    std::map<int, size_t> id_to_index_;
    
public:
    SpeciesParamsManager() = default;
    
    // Load species parameters from file
    bool loadFromFile(const std::string& filename);
    
    // Get species profile by ID
    const SpeciesProfile& getById(int species_id) const;
    
    // Get species profile by index
    const SpeciesProfile& getByIndex(size_t index) const;
    
    // Get all profiles
    const std::vector<SpeciesProfile>& getAllProfiles() const { return profiles_; }
    
    // Get number of species
    size_t getNumSpecies() const { return profiles_.size(); }
    
    // Get number of tree species (excluding herb)
    size_t getNumTreeSpecies() const;
    
    // Check if species exists
    bool hasSpecies(int species_id) const;
    
    // Get herb profile (species_id = 0)
    const SpeciesProfile& getHerbProfile() const { return getById(0); }
    
    // Get all tree species IDs
    std::vector<int> getTreeSpeciesIds() const;
};

#endif // SPECIES_PARAMS_H