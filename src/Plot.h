#ifndef PLOT_H
#define PLOT_H

#include <list>
#include <vector>
#include <string>
#include <memory>
#include <map>
#include "Config.h"
#include "AdultTree.h"
#include "Sapling.h"
#include "HerbLayer.h"
#include "SeedBank.h"
#include "SpatialHash.h"
#include "CanopyProjectionMap.h"
#include "SpeciesParams.h"

// Forward declaration
class WaterModule;
struct WaterBalanceResult;

// Site information structure
struct SiteInfo {
    int plot_id;
    std::string plot_name;
    std::string region_group;
    int climate_id;
    double elevation;
    double WHC_total;
    
    SiteInfo();
};

class Plot {
public:
    // Site Information
    SiteInfo site_info;
    
    // CRITICAL: Using std::list to prevent pointer invalidation
    // when elements are added/removed during simulation
    std::list<AdultTree> trees;
    std::list<Sapling> saplings;
    
    HerbLayer herb_layer;
    SeedBank seed_bank;
    
    // Spatial Acceleration Structures
    SpatialHash spatial_hash;
    CanopyProjectionMap canopy_map;
    
    // Water Balance State
    std::unique_ptr<WaterModule> water_module;
    
    // State Flags
    bool is_in_fire_recovery;
    int fire_recovery_year;
    
    // ID Counters
    int next_tree_id;
    int next_sapling_id;
    
    // Constructor
    explicit Plot(const SpeciesParamsManager& sp);
    
    // Destructor
    ~Plot();
    
    // Initialization
    void initialize(const SiteInfo& info);
    void reset();
    
    // Fire Event Handling
    void burn(int year);
    
    // Spatial Hash Management
    void rebuildSpatialHash();
    void rebuildCanopyMap();
    void insertTreeToHash(AdultTree* tree);
    void removeTreeFromHash(int tree_id);
    
    // Statistics - Single biomass method (no rough estimates)
    double calcTotalLAI(const SpeciesParamsManager& sp) const;
    double calcTotalBiomass(const SpeciesParamsManager& sp) const;
    double calcBiomassPerHa(const SpeciesParamsManager& sp) const;
    double calcTreeDensityPerHa() const;
    double calcSaplingDensityPerHa() const;
    double calcSourceStrength() const;
    double getMeanFloorLight() const;
    
    // Crown cover calculation (for DataOutput)
    double calcTotalCrownCover() const;
    std::map<int, double> calcCrownCoverBySpecies() const;
    
    // Species Statistics
    std::map<int, int> countTreesBySpecies() const;
    std::map<int, int> countSaplingsBySpecies() const;
    std::map<int, double> calcLeafAreaBySpecies() const;  // For LAI output
    std::map<int, double> calcMeanDBHBySpecies() const;
    std::map<int, double> calcMeanHeightBySpecies() const;
    
    // ID Management
    int getNextTreeId();
    int getNextSaplingId();
    void setNextTreeId(int id);
    void setNextSaplingId(int id);
};

// Site Loader
class SiteLoader {
public:
    static std::vector<SiteInfo> loadFromFile(const std::string& filename);
};

#endif // PLOT_H
