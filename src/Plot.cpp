#include "Plot.h"
#include "WaterModule.h"
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>

SiteInfo::SiteInfo() : plot_id(0), plot_name(""), region_group(""), 
                        climate_id(0), elevation(0), WHC_total(0) {}

Plot::Plot(const SpeciesParamsManager& sp) 
    : is_in_fire_recovery(false), fire_recovery_year(-1),
      next_tree_id(1), next_sapling_id(1) {
    water_module = std::make_unique<WaterModule>(sp);
}

Plot::~Plot() = default;

void Plot::initialize(const SiteInfo& info) {
    site_info = info;
    water_module->initialize(info.WHC_total);
    reset();
}

void Plot::reset() {
    trees.clear();
    saplings.clear();
    herb_layer.reset();
    seed_bank.clear();
    spatial_hash.clear();
    canopy_map.clear();
    water_module->reset();
    is_in_fire_recovery = false;
}

void Plot::burn(int year) {
    // Clear all vegetation
    trees.clear();
    saplings.clear();
    
    // Reset herb to minimum
    herb_layer.reset();
    
    // Clear seed bank
    seed_bank.clear();
    
    // Clear spatial structures
    spatial_hash.clear();
    canopy_map.clear();
    
    // Reset soil moisture
    water_module->reset();
    
    // Set recovery flag
    is_in_fire_recovery = true;
    fire_recovery_year = year;
}

void Plot::rebuildSpatialHash() {
    spatial_hash.rebuild(trees);
}

void Plot::rebuildCanopyMap() {
    canopy_map.build(trees);
}

void Plot::insertTreeToHash(AdultTree* tree) {
    spatial_hash.insert(tree);
    canopy_map.addTree(tree);
}

void Plot::removeTreeFromHash(int tree_id) {
    spatial_hash.remove(tree_id);
    canopy_map.removeTree(tree_id);
}

// MODIFIED: Conditionally include herb LAI based on enabled state
double Plot::calcTotalLAI(const SpeciesParamsManager& sp) const {
    double total_LA = 0;
    
    // Trees
    for (const auto& tree : trees) {
        total_LA += tree.leaf_area;
    }
    
    // Saplings
    for (const auto& sapling : saplings) {
        total_LA += sapling.calcLeafArea(sp.getById(sapling.species_id));
    }
    
    // Herbs - MODIFIED: Only include if competition is enabled
    if (herb_layer.isEnabled()) {
        for (int u = 0; u < Config::GRID_DIM; ++u) {
            for (int v = 0; v < Config::GRID_DIM; ++v) {
                total_LA += herb_layer.calcLAI(u, v);
            }
        }
    }
    
    double plot_area = Config::PLOT_SIZE * Config::PLOT_SIZE;
    return total_LA / plot_area;
}

// Single biomass calculation method using proper allometry
double Plot::calcTotalBiomass(const SpeciesParamsManager& sp) const {
    double total = 0;
    
    for (const auto& tree : trees) {
        total += tree.calcBiomass(sp.getById(tree.species_id));
    }
    
    // Add herb biomass
    total += herb_layer.getTotalBiomass();
    
    return total;
}

double Plot::calcBiomassPerHa(const SpeciesParamsManager& sp) const {
    double biomass_kg = calcTotalBiomass(sp);
    double plot_area_ha = (Config::PLOT_SIZE * Config::PLOT_SIZE) / 10000.0;
    return (biomass_kg / 1000.0) / plot_area_ha;  // kg to t, then per ha
}

double Plot::calcTreeDensityPerHa() const {
    double plot_area_ha = (Config::PLOT_SIZE * Config::PLOT_SIZE) / 10000.0;
    return trees.size() / plot_area_ha;
}

double Plot::calcSaplingDensityPerHa() const {
    double plot_area_ha = (Config::PLOT_SIZE * Config::PLOT_SIZE) / 10000.0;
    return saplings.size() / plot_area_ha;
}

double Plot::calcSourceStrength() const {
    if (trees.empty()) return 0;
    
    double total_LA = 0;
    for (const auto& tree : trees) {
        total_LA += tree.leaf_area;
    }
    
    double plot_area = Config::PLOT_SIZE * Config::PLOT_SIZE;
    return std::min(total_LA / (plot_area * Config::LAI_SAT), 1.0);
}

double Plot::getMeanFloorLight() const {
    return herb_layer.getMeanFloorLight();
}

double Plot::calcTotalCrownCover() const {
    double total_crown_area = 0;
    for (const auto& tree : trees) {
        total_crown_area += M_PI * tree.crown_radius * tree.crown_radius;
    }
    double plot_area = Config::PLOT_SIZE * Config::PLOT_SIZE;
    return std::min(1.0, total_crown_area / plot_area);
}

std::map<int, double> Plot::calcCrownCoverBySpecies() const {
    std::map<int, double> cc;
    double plot_area = Config::PLOT_SIZE * Config::PLOT_SIZE;
    
    for (const auto& tree : trees) {
        double crown_area = M_PI * tree.crown_radius * tree.crown_radius;
        cc[tree.species_id] += crown_area / plot_area;
    }
    
    return cc;
}

std::map<int, int> Plot::countTreesBySpecies() const {
    std::map<int, int> counts;
    for (const auto& tree : trees) {
        counts[tree.species_id]++;
    }
    return counts;
}

std::map<int, int> Plot::countSaplingsBySpecies() const {
    std::map<int, int> counts;
    for (const auto& sapling : saplings) {
        counts[sapling.species_id]++;
    }
    return counts;
}

std::map<int, double> Plot::calcLeafAreaBySpecies() const {
    std::map<int, double> la;
    for (const auto& tree : trees) {
        la[tree.species_id] += tree.leaf_area;
    }
    return la;
}

std::map<int, double> Plot::calcMeanDBHBySpecies() const {
    std::map<int, double> sum_dbh;
    std::map<int, int> counts;
    
    for (const auto& tree : trees) {
        sum_dbh[tree.species_id] += tree.dbh;
        counts[tree.species_id]++;
    }
    
    std::map<int, double> mean;
    for (const auto& pair : sum_dbh) {
        mean[pair.first] = pair.second / counts[pair.first];
    }
    return mean;
}

std::map<int, double> Plot::calcMeanHeightBySpecies() const {
    std::map<int, double> sum_h;
    std::map<int, int> counts;
    
    for (const auto& tree : trees) {
        sum_h[tree.species_id] += tree.height;
        counts[tree.species_id]++;
    }
    
    std::map<int, double> mean;
    for (const auto& pair : sum_h) {
        mean[pair.first] = pair.second / counts[pair.first];
    }
    return mean;
}

int Plot::getNextTreeId() { return next_tree_id++; }
int Plot::getNextSaplingId() { return next_sapling_id++; }
void Plot::setNextTreeId(int id) { next_tree_id = id; }
void Plot::setNextSaplingId(int id) { next_sapling_id = id; }

// SiteLoader implementation
std::vector<SiteInfo> SiteLoader::loadFromFile(const std::string& filename) {
    std::vector<SiteInfo> sites;
    std::ifstream file(filename);
    
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open site file: " + filename);
    }
    
    std::string line;
    std::getline(file, line);  // Skip header
    
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        
        std::istringstream iss(line);
        std::vector<std::string> tokens;
        std::string token;
        
        char delim = (line.find(',') != std::string::npos) ? ',' : '\t';
        while (std::getline(iss, token, delim)) {
            tokens.push_back(token);
        }
        
        if (tokens.size() >= 6) {
            try {
                SiteInfo site;
                site.plot_id = std::stoi(tokens[0]);
                site.plot_name = tokens[1];
                site.region_group = tokens[2];
                site.climate_id = std::stoi(tokens[3]);
                site.elevation = std::stod(tokens[4]);
                site.WHC_total = std::stod(tokens[5]);
                sites.push_back(site);
            } catch (...) {
                continue;
            }
        }
    }
    
    return sites;
}
