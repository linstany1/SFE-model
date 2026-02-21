#ifndef LIGHT_MODULE_H
#define LIGHT_MODULE_H

#include <list>
#include <vector>
#include "Config.h"
#include "SpeciesParams.h"

class AdultTree;
class Sapling;
class HerbLayer;
class SpatialHash;
class CanopyProjectionMap;

class LightModule {
private:
    const SpeciesParamsManager& species_params_;
    
public:
    explicit LightModule(const SpeciesParamsManager& sp);
    
    // Beer-Lambert light calculation
    double calcAvailableLight(double x, double y, double z,
                               const std::vector<AdultTree*>& nearby_trees,
                               double k_e = Config::DEFAULT_EXTINCTION_KE) const;
    
    // Adult tree light
    void calcAdultTreeLight(std::list<AdultTree>& trees,
                            const SpatialHash& spatial_hash) const;
    
    // Crown base height update (self-pruning)
    void updateCrownBaseHeights(std::list<AdultTree>& trees,
                                 const SpatialHash& spatial_hash) const;
    
    // Herb layer light
    void calcHerbLayerLight(HerbLayer& herb_layer,
                            const CanopyProjectionMap& canopy_map,
                            double k_e = Config::DEFAULT_EXTINCTION_KE) const;
    
    // Sapling light
    void calcSaplingLight(std::list<Sapling>& saplings,
                          const HerbLayer& herb_layer,
                          const CanopyProjectionMap& canopy_map) const;
    
    // Establishment light
    double calcEstablishmentLight(int u, int v, 
                                   const HerbLayer& herb_layer,
                                   double k_e) const;
    
    // Light response function (clamps only at outer level)
    static double calcLightResponse(double AL, int shade_tol);
    
    // Complete light calculation pipeline
    void runLightCalculation(std::list<AdultTree>& trees,
                             std::list<Sapling>& saplings,
                             HerbLayer& herb_layer,
                             SpatialHash& spatial_hash,
                             CanopyProjectionMap& canopy_map);
};

#endif // LIGHT_MODULE_H
