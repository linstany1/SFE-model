#include "SeedDispersalModule.h"
#include "AdultTree.h"
#include "SeedBank.h"
#include "DistanceMatrix.h"
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>

ExternalSourceInfo::ExternalSourceInfo() : plot_id(0), species_id(0), composition_ratio(0) {}
ExternalSourceInfo::ExternalSourceInfo(int p, int s, double r) : plot_id(p), species_id(s), composition_ratio(r) {}

SeedDispersalModule::SeedDispersalModule(const SpeciesParamsManager& sp, const DistanceMatrix& dm)
    : species_params_(sp), distance_matrix_(dm) {}

bool SeedDispersalModule::loadExternalSources(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) return false;
    
    std::string line;
    std::getline(file, line);
    
    external_sources_.clear();
    
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        
        std::istringstream iss(line);
        std::vector<std::string> tokens;
        std::string token;
        
        char delim = (line.find(',') != std::string::npos) ? ',' : '\t';
        while (std::getline(iss, token, delim)) {
            tokens.push_back(token);
        }
        
        if (tokens.size() >= 3) {
            try {
                ExternalSourceInfo src;
                src.plot_id = std::stoi(tokens[0]);
                src.species_id = std::stoi(tokens[1]);
                src.composition_ratio = std::stod(tokens[2]);
                external_sources_.push_back(src);
            } catch (...) {
                continue;
            }
        }
    }
    
    return true;
}

// Dynamic b parameter based on lambda2
double SeedDispersalModule::dispersalKernel(double distance, double lambda1, double lambda2) {
    if (distance <= lambda1) {
        return Config::SEED_NEAR_PROB;  // 0.95
    } else if (distance <= lambda2) {
        // Dynamic b: b = 2.0 if lambda2 < 100, else b = 5.0
        double b = (lambda2 < 100.0) ? 2.0 : 5.0;
        return std::exp(-b * distance / lambda2);
    } else {
        return Config::SEED_BACKGROUND_PROB;  // 0.001
    }
}

double SeedDispersalModule::calcTreeSeedProduction(const AdultTree& tree, const SpeciesProfile& sp) {
    if (tree.age < sp.maturity_age_yr) return 0;
    return sp.fecundity_f * tree.leaf_area;
}

double SeedDispersalModule::calcPlotSourceStrength(const std::list<AdultTree>& trees) {
    double total_LA = 0;
    for (const auto& tree : trees) {
        total_LA += tree.leaf_area;
    }
    
    double area = Config::PLOT_SIZE * Config::PLOT_SIZE;
    return std::min(total_LA / (area * Config::LAI_SAT), 1.0);
}

void SeedDispersalModule::calcWithinPlotDispersal(const std::list<AdultTree>& trees,
                                  SeedBank& seed_bank,
                                  int /*target_plot_id*/) const {
    for (const auto& tree : trees) {
        const SpeciesProfile& sp = species_params_.getById(tree.species_id);
        
        if (tree.age < sp.maturity_age_yr) continue;
        
        double Q_seed = sp.fecundity_f * tree.leaf_area;
        if (Q_seed <= 0) continue;
        
        for (int u = 0; u < Config::GRID_DIM; ++u) {
            for (int v = 0; v < Config::GRID_DIM; ++v) {
                double cx = u + 0.5;
                double cy = v + 0.5;
                
                double d = std::hypot(tree.x - cx, tree.y - cy);
                
                double K = dispersalKernel(d, sp.lambda1_m, sp.lambda2_m);
                
                double arrival = Q_seed * K;
                
                seed_bank.addSeeds(u, v, tree.species_id, arrival);
            }
        }
    }
}

void SeedDispersalModule::calcBetweenPlotDispersal(int source_plot_id,
                                   const std::list<AdultTree>& source_trees,
                                   int target_plot_id,
                                   SeedBank& seed_bank) const {
    double D = distance_matrix_.getDistance(source_plot_id, target_plot_id);
    
    if (DistanceMatrix::isIsolated(D) || source_plot_id == target_plot_id) {
        return;
    }
    
    double S_source = calcPlotSourceStrength(source_trees);
    if (S_source <= 0) return;
    
    for (const auto& tree : source_trees) {
        const SpeciesProfile& sp = species_params_.getById(tree.species_id);
        
        if (tree.age < sp.maturity_age_yr) continue;
        
        double Q_seed = sp.fecundity_f * tree.leaf_area;
        
        double K = dispersalKernel(D, sp.lambda1_m, sp.lambda2_m);
        
        double N_imm = Q_seed * K * S_source / 
                       (Config::PLOT_SIZE * Config::PLOT_SIZE);
        
        for (int u = 0; u < Config::GRID_DIM; ++u) {
            for (int v = 0; v < Config::GRID_DIM; ++v) {
                seed_bank.addSeeds(u, v, tree.species_id, N_imm);
            }
        }
    }
}

// [Phase 2c v1] External seed dispersal: probability arrival model.
//
// OLD MODEL (Bug 1.1):
//   S_potential = ratio × fecundity × LAI_SAT × PLOT_SIZE²  (天文数字)
//   N_ext = S_potential × K / PLOT_SIZE²  →  ~20 seeds/cell/source
//   With 8 sources → P_seed saturates to 1.0 for all species
//
// NEW MODEL:
//   Each external source provides an independent probability of seed arrival.
//   For species s, from all sources j:
//     P_arrive_s = 1 - Π_j(1 - ratio_j × K(D_j))
//   Then inject N_equiv = N_UNLIMITED × P_arrive_s uniformly.
//
//   Physical meaning: ratio_j is the fraction of the source landscape
//   occupied by species s. K(D_j) is the probability that a single seed
//   packet from distance D_j reaches the target plot. Their product is the
//   probability of arrival from source j. Independent sources combine via
//   the complement product.
//
//   Effect: Abge from distant sources (D=200m, K≈0.001) → P_arrive ≈ 0.01
//           Abge from nearby sources (D=80m, K≈0.35)   → P_arrive ≈ 0.60
//           Eliminates seed quantity saturation; pioneer species no longer
//           overwhelmed by fecundity × area scaling.
//
// Zero new parameters: reuses dispersalKernel, composition_ratio, N_UNLIMITED.
void SeedDispersalModule::calcExternalSourceDispersal(int target_plot_id,
                                      SeedBank& seed_bank) const {
    // Step 1: Collect all unique species across external sources
    std::set<int> all_species;
    for (const auto& src : external_sources_) {
        if (src.composition_ratio > 0 && species_params_.hasSpecies(src.species_id)) {
            all_species.insert(src.species_id);
        }
    }
    
    // Step 2: For each species, accumulate independent arrival probability
    // across all external source plots
    for (int sp_id : all_species) {
        const SpeciesProfile& sp = species_params_.getById(sp_id);
        
        // Accumulate survival probability (= probability NO seed arrives)
        double P_no_arrival = 1.0;
        
        for (const auto& src : external_sources_) {
            if (src.species_id != sp_id) continue;
            if (src.composition_ratio <= 0) continue;
            
            double D = distance_matrix_.getDistance(src.plot_id, target_plot_id);
            if (DistanceMatrix::isIsolated(D)) continue;
            
            double K = dispersalKernel(D, sp.lambda1_m, sp.lambda2_m);
            
            // Probability of arrival from this source
            // = composition_ratio × kernel probability
            // Clamped to [0, 1] for safety
            double P_from_source = std::min(1.0, src.composition_ratio * K);
            
            // Complement multiplication: P(no arrival from any) = Π(1 - P_j)
            P_no_arrival *= (1.0 - P_from_source);
        }
        
        double P_arrive = 1.0 - P_no_arrival;
        
        if (P_arrive < Config::EPSILON) continue;
        
        // Inject equivalent seed count: N_equiv = N_UNLIMITED × P_arrive
        // This maps [0,1] arrival probability to [0, N_UNLIMITED] seeds,
        // which then enters the standard P_seed = min(N/N_UNLIMITED, 1) pathway.
        double N_equiv = Config::N_UNLIMITED * P_arrive;
        
        for (int u = 0; u < Config::GRID_DIM; ++u) {
            for (int v = 0; v < Config::GRID_DIM; ++v) {
                seed_bank.addSeeds(u, v, sp_id, N_equiv);
            }
        }
    }
}

void SeedDispersalModule::applySaturatedSeedRain(SeedBank& seed_bank,
                                 const std::vector<int>& species_ids) const {
    for (int species_id : species_ids) {
        if (!species_params_.hasSpecies(species_id)) continue;
        
        for (int u = 0; u < Config::GRID_DIM; ++u) {
            for (int v = 0; v < Config::GRID_DIM; ++v) {
                seed_bank.addSeeds(u, v, species_id, Config::N_UNLIMITED);
            }
        }
    }
}

double SeedDispersalModule::calcSeedProb(double N_total) {
    return std::min(N_total / Config::N_UNLIMITED, 1.0);
}

std::vector<int> SeedDispersalModule::getExternalSpeciesIds() const {
    std::vector<int> ids;
    for (const auto& src : external_sources_) {
        if (std::find(ids.begin(), ids.end(), src.species_id) == ids.end()) {
            ids.push_back(src.species_id);
        }
    }
    return ids;
}

std::vector<int> SeedDispersalModule::getSourcesForTarget(int target_plot_id) const {
    std::vector<int> sources;
    for (int source_id : distance_matrix_.getSourceIds()) {
        if (distance_matrix_.areConnected(source_id, target_plot_id)) {
            sources.push_back(source_id);
        }
    }
    return sources;
}

double SeedDispersalModule::getTotalInternalSeeds(int species_id, 
                                  const SeedBank& seed_bank) const {
    return seed_bank.getTotalSeedsForSpecies(species_id);
}

std::vector<int> SeedDispersalModule::getRegionalSpeciesPool(
    int target_plot_id, 
    const std::vector<int>& target_plot_ids) const {
    
    std::set<int> species_set;  // Use set for automatic deduplication
    
    // 1. Get species from connected external seed sources
    for (const auto& src : external_sources_) {
        // Check if this external source is connected to target plot
        double dist = distance_matrix_.getDistance(src.plot_id, target_plot_id);
        if (!DistanceMatrix::isIsolated(dist) && src.composition_ratio > 0) {
            species_set.insert(src.species_id);
        }
    }
    
    // 2. Get species from connected target plots (they may share species pools)
    for (int other_plot_id : target_plot_ids) {
        if (other_plot_id == target_plot_id) continue;
        
        double dist = distance_matrix_.getDistance(other_plot_id, target_plot_id);
        if (!DistanceMatrix::isIsolated(dist)) {
            // Add species from external sources connected to this other plot
            for (const auto& src : external_sources_) {
                double src_dist = distance_matrix_.getDistance(src.plot_id, other_plot_id);
                if (!DistanceMatrix::isIsolated(src_dist) && src.composition_ratio > 0) {
                    species_set.insert(src.species_id);
                }
            }
        }
    }
    
    return std::vector<int>(species_set.begin(), species_set.end());
}