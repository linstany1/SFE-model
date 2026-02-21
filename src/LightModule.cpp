#include "LightModule.h"
#include "AdultTree.h"
#include "Sapling.h"
#include "HerbLayer.h"
#include "SpatialHash.h"
#include "CanopyProjectionMap.h"
#include "FloatUtil.h"
#include <cmath>
#include <algorithm>
#include <array>

LightModule::LightModule(const SpeciesParamsManager& sp) : species_params_(sp) {}

double LightModule::calcAvailableLight(double x, double y, double z,
                               const std::vector<AdultTree*>& nearby_trees,
                               double k_e) const {
    double total_attenuation = 0.0;
    
    for (const AdultTree* tree : nearby_trees) {
        if (!tree) continue;
        
        const SpeciesProfile& sp = species_params_.getById(tree->species_id);
        
        // [Phase 4c(A) v1] Toroidal horizontal distance
        double dx = Config::toroidalDelta(x, tree->x);
        double dy = Config::toroidalDelta(y, tree->y);
        double d_horiz = std::sqrt(dx * dx + dy * dy);
        
        if (d_horiz >= tree->crown_radius) continue;
        
        double d_ratio = d_horiz / tree->crown_radius;
        double Z_in = tree->crown_base_height + 
                      (tree->height - tree->crown_base_height) * 
                      (1.0 - std::pow(d_ratio, 1.0 / sp.cs));
        
        double L_path = std::max(0.0, Z_in - std::max(z, tree->crown_base_height));
        
        total_attenuation += tree->leaf_area_density * L_path;
    }
    
    return std::exp(-k_e * total_attenuation);
}

void LightModule::calcAdultTreeLight(std::list<AdultTree>& trees,
                            const SpatialHash& spatial_hash) const {
    for (auto& tree : trees) {
        const SpeciesProfile& sp = species_params_.getById(tree.species_id);
        
        double Z_eff = tree.getEffectiveCrownHeight();
        
        std::vector<AdultTree*> neighbors = spatial_hash.getNearbyTrees(tree.x, tree.y);
        
        tree.available_light = calcAvailableLight(tree.x, tree.y, Z_eff, 
                                                   neighbors, sp.extinction_ke);
    }
}

void LightModule::updateCrownBaseHeights(std::list<AdultTree>& trees,
                                 const SpatialHash& spatial_hash) const {
    for (auto& tree : trees) {
        const SpeciesProfile& sp = species_params_.getById(tree.species_id);
        
        std::vector<AdultTree*> neighbors = spatial_hash.getNearbyTrees(tree.x, tree.y);
        
        double AL_top = calcAvailableLight(tree.x, tree.y, tree.height, 
                                            neighbors, sp.extinction_ke);
        
        double AL_HB = calcAvailableLight(tree.x, tree.y, tree.crown_base_height, 
                                           neighbors, sp.extinction_ke);
        
        AL_top = std::max(Config::MIN_LIGHT_FOR_LOG, AL_top);
        AL_HB = std::max(Config::MIN_LIGHT_FOR_LOG, AL_HB);
        
        // ============================================================
        // MODIFIED v2.1: Record HB before update for pruning loss calc
        // ============================================================
        double HB_old = tree.crown_base_height;
        
        double HB_open = (1.0 - sp.alpha_base) * tree.height;
        
        double z_star = tree.crown_base_height;
        
        if (AL_HB < sp.l_min && sp.l_min < AL_top) {
            double ln_AL_top = FloatUtil::safeLog(AL_top);
            double ln_AL_HB = FloatUtil::safeLog(AL_HB);
            double ln_L_min = FloatUtil::safeLog(sp.l_min);
            double ln_diff = ln_AL_top - ln_AL_HB;
            
            if (std::abs(ln_diff) > Config::EPSILON) {
                z_star = tree.crown_base_height + 
                         (tree.height - tree.crown_base_height) * 
                         (ln_L_min - ln_AL_HB) / ln_diff;
            }
        }
        
        double new_HB = std::max({tree.crown_base_height, HB_open, z_star});
        tree.crown_base_height = std::min(Config::MAX_CROWN_RATIO * tree.height, new_HB);
        
        // NEW v2.1: Apply foliage pruning loss if crown base moved up
        tree.applyPruningLoss(sp, HB_old);
        
        tree.updateLeafAreaDensity(sp);
    }
}

void LightModule::calcHerbLayerLight(HerbLayer& herb_layer,
                            const CanopyProjectionMap& canopy_map,
                            double k_e) const {
    for (int u = 0; u < Config::GRID_DIM; ++u) {
        for (int v = 0; v < Config::GRID_DIM; ++v) {
            double x = u + 0.5;
            double y = v + 0.5;
            
            const std::vector<AdultTree*>& shaders = canopy_map.getTreesAtCell(u, v);
            
            double AL = calcAvailableLight(x, y, Config::HERB_HEIGHT, shaders, k_e);
            
            herb_layer.setTopLight(u, v, AL);
        }
    }
}

void LightModule::calcSaplingLight(std::list<Sapling>& saplings,
                          const HerbLayer& herb_layer,
                          const CanopyProjectionMap& canopy_map) const {
    // [Phase 4c(B) v1] Pre-generate sapling LAI grid for peer shading.
    // Each sapling shades its cell-mates via Beer-Lambert.
    // Without this, 10+ saplings/m² all get identical light from above,
    // creating "ghost growth" where they don't compete with each other.
    std::array<std::array<double, Config::GRID_DIM>, Config::GRID_DIM> sapling_LAI_grid{};
    
    for (const auto& sapling : saplings) {
        int u = sapling.getGridU();
        int v = sapling.getGridV();
        if (u >= 0 && u < Config::GRID_DIM && v >= 0 && v < Config::GRID_DIM) {
            const SpeciesProfile& sp = species_params_.getById(sapling.species_id);
            sapling_LAI_grid[u][v] += sapling.calcLeafArea(sp);
        }
    }
    
    for (auto& sapling : saplings) {
        const SpeciesProfile& sp = species_params_.getById(sapling.species_id);
        
        if (sapling.isTallSapling()) {
            int u = sapling.getGridU();
            int v = sapling.getGridV();
            const std::vector<AdultTree*>& shaders = canopy_map.getTreesAtCell(u, v);
            
            sapling.available_light = calcAvailableLight(
                sapling.x, sapling.y, sapling.height, shaders, sp.extinction_ke);
        } else {
            int u = sapling.getGridU();
            int v = sapling.getGridV();
            
            sapling.available_light = herb_layer.calcLightForShortSapling(
                u, v, sapling.height, sp.extinction_ke);
        }
        
        // Apply peer shading: subtract self, attenuate by remainder
        int u = sapling.getGridU();
        int v = sapling.getGridV();
        if (u >= 0 && u < Config::GRID_DIM && v >= 0 && v < Config::GRID_DIM) {
            double self_LA = sapling.calcLeafArea(sp);
            double peer_LAI = std::max(0.0, sapling_LAI_grid[u][v] - self_LA);
            
            if (peer_LAI > Config::EPSILON) {
                sapling.available_light *= std::exp(-sp.extinction_ke * peer_LAI);
            }
        }
    }
}

double LightModule::calcEstablishmentLight(int u, int v, 
                                   const HerbLayer& herb_layer,
                                   double k_e) const {
    return herb_layer.calcLightForShortSapling(u, v, 
                                                Config::INITIAL_SAPLING_HEIGHT, k_e);
}

// Light response function (clamps only at outer level)
double LightModule::calcLightResponse(double AL, int shade_tol) {
    AL = std::max(0.0, AL);
    
    // f_intol (ST=1) larix
    double f_intol = 2.24 * (1.0 - std::exp(-1.136 * (AL - 0.08)));
    
    // f_tol (ST=5) abge
    double f_tol = 1.0 - std::exp(-4.68 * (AL - 0.012));
    
    double weight = (shade_tol - 1) / 4.0;
    double f_raw = f_intol + weight * (f_tol - f_intol);
    
    // IMPORTANT: Only clamp at the outer level
    return std::max(0.0, std::min(1.0, f_raw));
}

void LightModule::runLightCalculation(std::list<AdultTree>& trees,
                             std::list<Sapling>& saplings,
                             HerbLayer& herb_layer,
                             SpatialHash& spatial_hash,
                             CanopyProjectionMap& canopy_map) {
    spatial_hash.rebuild(trees);
    canopy_map.build(trees);
    
    calcAdultTreeLight(trees, spatial_hash);
    
    updateCrownBaseHeights(trees, spatial_hash);
    
    canopy_map.build(trees);
    
    double k_e = Config::DEFAULT_EXTINCTION_KE;
    if (species_params_.hasSpecies(0)) {
        k_e = species_params_.getHerbProfile().extinction_ke;
    }
    calcHerbLayerLight(herb_layer, canopy_map, k_e);
    
    calcSaplingLight(saplings, herb_layer, canopy_map);
}