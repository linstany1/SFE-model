/**
 * @file main.cpp
 * @brief SFE Model - Postfire Forest Succession
 * * Study Area: Gongga Mountain SW Slope
 * Climate: 1901-2021
 * Fire Event: 1954
 * Species: Herb + Abge + Pilil + Quaq + Lagr + Jusq
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <map>
#include <cmath>
#include <random>
#include <algorithm>
#include <chrono>

// Use the existing relative paths
#include "core/GlobalConfig.hpp"
#include "core/Types.hpp"
#include "core/Grid.hpp"
#include "core/RandomGenerator.hpp"
#include "entities/Plot.hpp"
#include "entities/AdultTree.hpp"
#include "entities/Sapling.hpp"
#include "io/SpeciesProfile.hpp"
#include "modules/LightEngine.hpp"
#include "modules/Growth.hpp"

using namespace sfe;

// ============================================================================
// Simulation Configuration
// ============================================================================
namespace SimConfig {
    // Spin-up Phase
    constexpr int SPIN_UP_YEARS = 500;        // Years for spin-up
    constexpr int SATURATED_YEARS = 50;       // Years of saturated seed rain
    constexpr int CLIMATE_WINDOW = 30;        // Climate resampling window (1901-1930)
    
    // Transient Simulation
    constexpr int START_YEAR = 1901;          // Start year
    constexpr int END_YEAR = 2021;            // End year
    
    // Output Control
    constexpr int OUTPUT_INTERVAL = 10;       // Output interval (years)
    constexpr bool VERBOSE_DIAGNOSTICS = true; 
    
    // Random Seed
    constexpr uint64_t RANDOM_SEED = 12345;
}

// ============================================================================
// Data Structures
// ============================================================================
struct MonthlyClimate {
    double temp = 0, temp_min = 0, temp_max = 0, prec = 0, pet = 0;
};

struct PlotConfig {
    int plot_id;
    std::string plot_name;
    int climate_id;
    std::string region;
    double elevation, latitude, whc, slope, aspect;
};

struct SeedSource {
    int plot_id;
    std::map<int, double> species_ratios;
};

// Global Containers
std::map<int, std::map<int, std::vector<MonthlyClimate>>> g_climate;
std::vector<PlotConfig> g_sites;
std::map<int, std::vector<int>> g_fire_events;  // year -> plot_ids
std::vector<std::vector<double>> g_dist_matrix;
std::vector<SeedSource> g_seed_sources;

// ============================================================================
// Data Loading Functions
// ============================================================================
bool loadSpeciesParams(const std::string& filename, SpeciesManager& mgr) {
    std::ifstream file(filename);
    if (!file.is_open()) { std::cerr << "Cannot open: " << filename << "\n"; return false; }
    
    auto parseD = [](const std::string& s, double d = 0.0) {
        if (s == "NA" || s.empty()) return d;
        try { return std::stod(s); } catch(...) { return d; }
    };
    auto parseI = [](const std::string& s, int d = 0) {
        if (s == "NA" || s.empty()) return d;
        try { return std::stoi(s); } catch(...) { return d; }
    };
    
    std::string line;
    std::getline(file, line); // header
    
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        std::istringstream iss(line);
        std::vector<std::string> f;
        std::string field;
        while (std::getline(iss, field, '\t')) f.push_back(field);
        if (f.size() < 44) continue;
        
        SpeciesProfile sp;
        sp.species_id = parseI(f[0]);
        sp.species_name = f[1];
        sp.is_evergreen = (parseI(f[2]) == 1);
        sp.lambda1_m = parseD(f[3], 30);
        sp.lambda2_m = parseD(f[4], 100);
        sp.k_ldd = parseD(f[5], 0.1);
        sp.fecundity_f = parseD(f[6], 10);
        sp.maturity_age_yr = parseI(f[7], 20);
        sp.ms_seed_kg = parseD(f[8], 0.0001);
        sp.h_max_m = parseD(f[9], 30);
        sp.dbh_max_cm = parseD(f[10], 100);
        sp.r_a = parseD(f[11], 0.5);
        sp.r_b = parseD(f[12], 0.3);
        sp.cs = parseD(f[13], 0.5);
        sp.alpha_base = parseD(f[14], 0.3);
        sp.ka1 = parseD(f[15], 0.1);
        sp.ka2 = parseD(f[16], 1.5);
        sp.kc2 = parseD(f[17], 6);
        sp.opt_s = parseD(f[18], 1);
        sp.ga = parseD(f[19], 1);
        sp.gs = parseD(f[20], 0.5);
        sp.dw = parseD(f[21], 0.4);
        sp.fw1 = parseD(f[22], 0.1);
        sp.fw2 = parseD(f[23], 1.5);
        sp.sa = parseD(f[24], 0.1);
        sp.sb = parseD(f[25], 2);
        sp.b0 = parseD(f[26], 0.1);
        sp.b1 = parseD(f[27], 2);
        sp.b2 = parseD(f[28], 0.5);
        sp.b0_small = parseD(f[29], 0.1);
        sp.b1_small = parseD(f[30], 1.5);
        sp.b2_small = parseD(f[31], 0.5);
        sp.shade_tol = parseI(f[32], 3);
        sp.dr_tol = parseD(f[33], 0.3);
        sp.dd_min = parseD(f[34], 500);
        sp.l_min = parseD(f[35], 0.05);
        sp.l_max = parseD(f[36], 1.0);
        sp.ngs_tmin_c = parseD(f[37], -30);
        sp.ngs_tmax_c = parseD(f[38], 10);
        sp.k_max_herb = parseD(f[39], 0.8);
        sp.r_max_herb = parseD(f[40], 1.2);
        sp.max_age_yr = parseI(f[41], 300);
        sp.extinction_ke = parseD(f[42], 0.5);
        mgr.addSpecies(sp);
    }
    return true;
}

bool loadSites(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) return false;
    std::string line;
    std::getline(file, line);
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        std::istringstream iss(line);
        PlotConfig p;
        iss >> p.plot_id >> p.plot_name >> p.climate_id >> p.region
            >> p.elevation >> p.latitude >> p.whc >> p.slope >> p.aspect;
        g_sites.push_back(p);
    }
    return true;
}

bool loadClimate(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) return false;
    std::string line;
    std::getline(file, line);
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        std::istringstream iss(line);
        int cid, year, month;
        MonthlyClimate mc;
        iss >> cid >> year >> month >> mc.temp >> mc.temp_min >> mc.temp_max >> mc.prec >> mc.pet;
        if (g_climate[cid][year].empty()) g_climate[cid][year].resize(12);
        if (month >= 1 && month <= 12) g_climate[cid][year][month-1] = mc;
    }
    return true;
}

bool loadEvents(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) return false;
    std::string line;
    std::getline(file, line);
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        std::istringstream iss(line);
        int year, pid;
        std::string type;
        double intensity;
        iss >> year >> pid >> type >> intensity;
        if (type == "Fire") g_fire_events[year].push_back(pid);
    }
    return true;
}

bool loadDistMatrix(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) return false;
    std::string line;
    std::getline(file, line);
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        std::istringstream iss(line);
        int src;
        iss >> src;
        std::vector<double> row;
        double d;
        while (iss >> d) row.push_back(d);
        g_dist_matrix.push_back(row);
    }
    return true;
}

bool loadSeedSources(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) return false;
    std::map<int, SeedSource> src_map;
    std::string line;
    std::getline(file, line);
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        std::istringstream iss(line);
        int pid, sid;
        double ratio;
        iss >> pid >> sid >> ratio;
        src_map[pid].plot_id = pid;
        src_map[pid].species_ratios[sid] = ratio;
    }
    for (auto& [id, src] : src_map) g_seed_sources.push_back(src);
    return true;
}

// ============================================================================
// Climate Functions
// ============================================================================
double calcGDD(const std::vector<MonthlyClimate>& m, bool evergreen = true) {
    constexpr double T_BASE = 5.5;
    constexpr int DAYS[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    double gdd = 0;
    for (int i = 0; i < 12; ++i) {
        double t = m[i].temp;
        if (evergreen || t >= T_BASE) {
            gdd += DAYS[i] * std::max(0.0, t - T_BASE);
        }
    }
    return gdd;
}

double calcNGSTemp(const std::vector<MonthlyClimate>& m) {
    double sum = 0;
    int cnt = 0;
    for (int i = 0; i < 12; ++i) {
        if (m[i].temp < 5.5) {
            sum += m[i].temp;
            cnt++;
        }
    }
    return cnt > 0 ? sum / cnt : 0;
}

double calcPrecip(const std::vector<MonthlyClimate>& m) {
    double total = 0;
    for (const auto& mc : m) total += mc.prec;
    return total;
}

double calcDroughtIndex(const std::vector<MonthlyClimate>& m) {
    double precip = calcPrecip(m);
    double pet = 0;
    for (const auto& mc : m) pet += mc.pet;
    if (pet <= 0) return 0;
    return std::max(0.0, std::min(1.0, 1.0 - precip / pet));
}

// ============================================================================
// Main Function
// ============================================================================
int main() {
    // Pure ASCII Header - Safe for all terminals
    std::cout << "+----------------------------------------------------------+\n";
    std::cout << "|      SFE Model - Postfire Forest Succession              |\n";
    std::cout << "|      Gongga Mountain Simulation (Phase 9)                |\n";
    std::cout << "+----------------------------------------------------------+\n\n";
    
    auto t_start = std::chrono::high_resolution_clock::now();
    std::string data_dir = "./data/";
    
    // ========================================================================
    // Step 1: Load Data
    // ========================================================================
    std::cout << "> Step 1: Loading Data Files\n" << std::string(50, '-') << "\n";
    
    SpeciesManager species_mgr;
    if (!loadSpeciesParams(data_dir + "species_params.txt", species_mgr)) return 1;
    std::cout << "  * Species params: " << species_mgr.getSpeciesCount() << " species\n";
    
    if (!loadSites(data_dir + "site.txt")) return 1;
    std::cout << "  * Site config: " << g_sites.size() << " plots\n";
    
    if (!loadClimate(data_dir + "climate.txt")) return 1;
    std::cout << "  * Climate data: " << g_climate.size() << " sequences\n";
    
    if (!loadEvents(data_dir + "events.txt")) return 1;
    std::cout << "  * Events: " << g_fire_events.size() << " fire events\n";
    
    if (!loadDistMatrix(data_dir + "distance_matrix.txt")) return 1;
    std::cout << "  * Distance matrix: " << g_dist_matrix.size() << " sources\n";
    
    if (!loadSeedSources(data_dir + "seed_source_composition.txt")) return 1;
    std::cout << "  * Seed sources: " << g_seed_sources.size() << " plots\n";
    
    // Diagnostics
    if (SimConfig::VERBOSE_DIAGNOSTICS) {
        std::cout << "\n  Species Parameters:\n";
        std::cout << "  " << std::string(70, '-') << "\n";
        std::cout << "  " << std::left << std::setw(8) << "ID" 
                  << std::setw(10) << "Name" << std::setw(10) << "dd_min"
                  << std::setw(10) << "ngs_tmin" << std::setw(10) << "ngs_tmax"
                  << std::setw(8) << "l_min" << std::setw(8) << "l_max" << "\n";
        for (size_t i = 1; i < species_mgr.getSpeciesCount(); ++i) {
            const auto& sp = species_mgr.getSpecies(i);
            std::cout << "  " << std::setw(8) << sp.species_id 
                      << std::setw(10) << sp.species_name
                      << std::setw(10) << sp.dd_min
                      << std::setw(10) << sp.ngs_tmin_c
                      << std::setw(10) << sp.ngs_tmax_c
                      << std::setw(8) << sp.l_min
                      << std::setw(8) << sp.l_max << "\n";
        }
    }
    
    if (SimConfig::VERBOSE_DIAGNOSTICS && !g_climate.empty()) {
        std::cout << "\n  Climate Check (1901):\n";
        std::cout << "  " << std::string(50, '-') << "\n";
        for (const auto& [cid, years] : g_climate) {
            if (years.find(1901) != years.end()) {
                const auto& m = years.at(1901);
                double gdd = calcGDD(m);
                double ngs = calcNGSTemp(m);
                double di = calcDroughtIndex(m);
                std::cout << "  Climate " << cid << ": GDD=" << std::fixed << std::setprecision(0) << gdd
                          << ", NGS_Temp=" << std::setprecision(1) << ngs << " C"
                          << ", DroughtIndex=" << std::setprecision(2) << di << "\n";
            }
        }
    }
    
    // ========================================================================
    // Step 2: Create Plots
    // ========================================================================
    std::cout << "\n> Step 2: Create Plots\n" << std::string(50, '-') << "\n";
    
    std::vector<Plot> plots;
    std::vector<int> plot_climate_ids;
    std::vector<std::string> plot_names;
    
    for (const auto& s : g_sites) {
        SiteInfo si(s.plot_id, s.plot_name, s.climate_id, s.region,
                    s.elevation, s.latitude, s.whc, s.slope, s.aspect);
        Plot plot(si);
        plots.push_back(std::move(plot));
        plot_climate_ids.push_back(s.climate_id);
        plot_names.push_back(s.plot_name);
        std::cout << "  * Plot " << s.plot_id << " (" << s.plot_name 
                  << "): Elev=" << s.elevation << "m, climate_id=" << s.climate_id << "\n";
    }
    
    // ========================================================================
    // Step 3: Init Engine
    // ========================================================================
    std::cout << "\n> Step 3: Initialize Engine\n" << std::string(50, '-') << "\n";
    
    LightEngine light_engine;
    RandomGenerator& rng = RandomGenerator::getInstance();
    rng.setSeed(SimConfig::RANDOM_SEED);
    std::cout << "  * LightEngine initialized\n";
    std::cout << "  * RandomGenerator seed=" << SimConfig::RANDOM_SEED << "\n";
    
    // Build climate year list
    std::vector<int> base_years;
    for (const auto& [yr, data] : g_climate.begin()->second) {
        if (static_cast<int>(base_years.size()) < SimConfig::CLIMATE_WINDOW)
            base_years.push_back(yr);
    }
    std::cout << "  * Climate window: " << base_years.front() << "-" << base_years.back() << "\n";
    
    TreeId next_tree_id = 1, next_sap_id = 10000;
    
    // ========================================================================
    // Step 4: Spin-up
    // ========================================================================
    std::cout << "\n> Step 4: Spin-up Phase\n";
    std::cout << "  Total years: " << SimConfig::SPIN_UP_YEARS << " years\n";
    std::cout << "  Saturated seed rain: First " << SimConfig::SATURATED_YEARS << " years\n";
    std::cout << "  Fire disturbance: Disabled\n";
    std::cout << std::string(50, '=') << "\n";
    
    int establishment_count = 0;
    int recruitment_count = 0;
    
    for (int step = 1; step <= SimConfig::SPIN_UP_YEARS; ++step) {
        int rand_year = base_years[rng.getInt(0, base_years.size() - 1)];
        bool saturated = (step <= SimConfig::SATURATED_YEARS);
        
        for (size_t p = 0; p < plots.size(); ++p) {
            auto& plot = plots[p];
            int cid = plot_climate_ids[p];
            
            if (g_climate.find(cid) == g_climate.end()) continue;
            if (g_climate[cid].find(rand_year) == g_climate[cid].end()) continue;
            
            const auto& monthly = g_climate[cid][rand_year];
            double gdd = calcGDD(monthly);
            double ngs_temp = calcNGSTemp(monthly);
            double drought = calcDroughtIndex(monthly);
            
            // Light
            Grid<double> herb_lai(Plot::SIZE_M, Plot::SIZE_M, 0.2);
            light_engine.calculateAllLight(plot.getTrees(), &herb_lai);
            
            for (int y = 0; y < Plot::SIZE_M; ++y) {
                for (int x = 0; x < Plot::SIZE_M; ++x) {
                    plot.getCells().get(x, y).light_at_ground = light_engine.getGroundLightAt(x, y);
                }
            }
            
            // Establishment
            for (int y = 0; y < Plot::SIZE_M; ++y) {
                for (int x = 0; x < Plot::SIZE_M; ++x) {
                    double light = plot.getCells().get(x, y).light_at_ground;
                    
                    for (size_t sp_id = 1; sp_id < species_mgr.getSpeciesCount(); ++sp_id) {
                        const auto& sp = species_mgr.getSpecies(sp_id);
                        
                        bool pass_light = (light >= sp.l_min && light <= sp.l_max);
                        bool pass_ngs = (ngs_temp >= sp.ngs_tmin_c && ngs_temp < sp.ngs_tmax_c);
                        bool pass_gdd = (gdd >= sp.dd_min);
                        bool pass_drought = (drought <= sp.dr_tol);
                        
                        if (!pass_light || !pass_ngs || !pass_gdd || !pass_drought) continue;
                        
                        double p_est = saturated ? 0.2 : 0.05;
                        
                        if (!saturated) {
                            bool has_mature_tree = false;
                            for (const auto& tree : plot.getTrees()) {
                                if (tree.getSpeciesId() == static_cast<int>(sp_id) && 
                                    tree.getAge() >= sp.maturity_age_yr) {
                                    has_mature_tree = true;
                                    break;
                                }
                            }
                            if (!has_mature_tree) p_est = 0;
                        }
                        
                        if (p_est > 0 && rng.getUniform01() < p_est) {
                            Sapling sap(next_sap_id++, sp_id, 
                                        x + rng.getUniform01(), y + rng.getUniform01(), 0.05);
                            plot.getSaplings().push_back(sap);
                            establishment_count++;
                        }
                    }
                }
            }
            
            // Sapling Growth
            for (auto& sap : plot.getSaplings()) {
                if (!sap.isAlive()) continue;
                const auto& sp = species_mgr.getSpecies(sap.getSpeciesId());
                
                int cx = std::clamp(static_cast<int>(sap.getX()), 0, Plot::SIZE_M - 1);
                int cy = std::clamp(static_cast<int>(sap.getY()), 0, Plot::SIZE_M - 1);
                double light = light_engine.calcSaplingLight(sap, herb_lai.get(cx, cy));
                
                double f_light = EnvironmentResponse::calcLightResponse(light, sp.shade_tol);
                double f_temp = EnvironmentResponse::calcTemperatureResponse(gdd, sp.dd_min);
                double f_drought = EnvironmentResponse::calcDroughtResponse(drought, sp.dr_tol);
                double f_env = f_light * f_temp * f_drought;
                
                double dh = SaplingGrowthCalculator::calcHeightIncrement(
                    sap.getHeight(), sp.h_max_m, sp.gs, f_env);
                sap.grow(dh);
                sap.incrementAge();
            }
            
            // Adult Growth
            for (auto& tree : plot.getTrees()) {
                if (!tree.isAlive()) continue;
                const auto& sp = species_mgr.getSpecies(tree.getSpeciesId());
                
                double light = tree.getAvailableLight();
                double f_env = EnvironmentResponse::calcLightResponse(light, sp.shade_tol) *
                               EnvironmentResponse::calcTemperatureResponse(gdd, sp.dd_min) *
                               EnvironmentResponse::calcDroughtResponse(drought, sp.dr_tol);
                
                double ddbh = TreeGrowthCalculator::calcOptimalDbhIncrement(
                    tree.getDbh(), tree.getHeight(), sp.dbh_max_cm, sp.h_max_m, sp.ga, sp.opt_s) * f_env;
                
                tree.setDbh(tree.getDbh() + ddbh);
                tree.setHeight(TreeGrowthCalculator::calcHeightFromDbh(sp.h_max_m, sp.opt_s, tree.getDbh()));
                tree.incrementAge();
                
                double r = std::max(0.5, sp.r_a * std::log(std::max(0.1, tree.getDbh())) + sp.r_b);
                tree.setCrownRadius(r);
                tree.setLeafArea(sp.kc2 * sp.ka1 * std::pow(std::max(0.1, tree.getDbh()), sp.ka2));
            }
            
            // Recruitment
            for (auto it = plot.getSaplings().begin(); it != plot.getSaplings().end(); ) {
                if (it->isAlive() && it->getHeight() >= 1.37) {
                    const auto& sp = species_mgr.getSpecies(it->getSpeciesId());
                    AdultTree tree(next_tree_id++, it->getSpeciesId(), 
                                   it->getX(), it->getY(), 0.1, 1.37, it->getAge());
                    tree.setCrownRadius(std::max(0.5, sp.r_a * std::log(0.1) + sp.r_b));
                    tree.setCrownBase(sp.alpha_base * 1.37);
                    tree.setCrownShape(sp.cs);
                    tree.setLeafArea(sp.kc2 * sp.ka1 * std::pow(0.1, sp.ka2));
                    plot.getTrees().push_back(tree);
                    it = plot.getSaplings().erase(it);
                    recruitment_count++;
                } else {
                    ++it;
                }
            }
            
            // Mortality
            for (auto& sap : plot.getSaplings()) {
                if (sap.isAlive()) {
                    double p_mort = 0.02;
                    if (sap.getAge() > 30) p_mort = 0.05;
                    if (rng.getUniform01() < p_mort) sap.kill();
                }
            }
            
            for (auto& tree : plot.getTrees()) {
                if (tree.isAlive()) {
                    const auto& sp = species_mgr.getSpecies(tree.getSpeciesId());
                    double p_mort = 1.0 - std::pow(0.01, 1.0 / sp.max_age_yr);
                    if (rng.getUniform01() < p_mort) tree.kill();
                }
            }
            
            // Cleanup
            plot.getSaplings().erase(
                std::remove_if(plot.getSaplings().begin(), plot.getSaplings().end(),
                    [](const Sapling& s) { return !s.isAlive(); }),
                plot.getSaplings().end());
            
            plot.getTrees().erase(
                std::remove_if(plot.getTrees().begin(), plot.getTrees().end(),
                    [](const AdultTree& t) { return !t.isAlive(); }),
                plot.getTrees().end());
        }
        
        if (step % 100 == 0 || step == SimConfig::SPIN_UP_YEARS) {
            int tt = 0, ts = 0;
            double bio = 0;
            for (const auto& pl : plots) {
                tt += pl.getTrees().size();
                ts += pl.getSaplings().size();
                for (const auto& tr : pl.getTrees()) {
                    const auto& sp = species_mgr.getSpecies(tr.getSpeciesId());
                    bio += sp.b0 * std::pow(tr.getDbh(), sp.b1) * std::pow(tr.getHeight(), sp.b2);
                }
            }
            std::cout << "[Spin-up] Year " << std::setw(4) << step 
                      << " | Trees: " << std::setw(5) << tt
                      << " | Saplings: " << std::setw(5) << ts
                      << " | Biomass: " << std::fixed << std::setprecision(0) << bio << " kg\n";
        }
    }
    
    std::cout << "\nSpin-up Completed:\n";
    std::cout << "  Total Est.: " << establishment_count << "\n";
    std::cout << "  Total Recr.: " << recruitment_count << "\n";
    
    // ========================================================================
    // Step 5: Transient Simulation
    // ========================================================================
    std::cout << "\n> Step 5: Transient Simulation\n";
    std::cout << "  Time range: " << SimConfig::START_YEAR << " - " << SimConfig::END_YEAR << "\n";
    std::cout << "  Fire events: Enabled (1954)\n";
    std::cout << std::string(50, '=') << "\n";
    
    for (int year = SimConfig::START_YEAR; year <= SimConfig::END_YEAR; ++year) {
        bool fire_this_year = false;
        if (g_fire_events.count(year)) {
            fire_this_year = true;
            for (int pid : g_fire_events[year]) {
                for (size_t p = 0; p < plots.size(); ++p) {
                    if (static_cast<int>(p) == pid || pid == -1) {
                        plots[p].getTrees().clear();
                        plots[p].getSaplings().clear();
                        std::cout << "  [FIRE] Year " << year << ": Plot " << p 
                                  << " (" << plot_names[p] << ") ALL CLEARED!\n";
                    }
                }
            }
        }
        
        for (size_t p = 0; p < plots.size(); ++p) {
            auto& plot = plots[p];
            int cid = plot_climate_ids[p];
            
            if (g_climate.find(cid) == g_climate.end()) continue;
            if (g_climate[cid].find(year) == g_climate[cid].end()) continue;
            
            const auto& monthly = g_climate[cid][year];
            double gdd = calcGDD(monthly);
            double ngs_temp = calcNGSTemp(monthly);
            double drought = calcDroughtIndex(monthly);
            
            Grid<double> herb_lai(Plot::SIZE_M, Plot::SIZE_M, 0.2);
            light_engine.calculateAllLight(plot.getTrees(), &herb_lai);
            
            for (int y = 0; y < Plot::SIZE_M; ++y) {
                for (int x = 0; x < Plot::SIZE_M; ++x) {
                    plot.getCells().get(x, y).light_at_ground = light_engine.getGroundLightAt(x, y);
                }
            }
            
            // External Seeds (Post-fire)
            bool is_open = (plot.getTrees().size() < 10);
            if (is_open) {
                for (const auto& src : g_seed_sources) {
                    if (src.plot_id >= static_cast<int>(g_dist_matrix.size())) continue;
                    if (p >= g_dist_matrix[src.plot_id].size()) continue;
                    double dist = g_dist_matrix[src.plot_id][p];
                    if (dist < 0) continue;
                    
                    for (const auto& [sp_id, ratio] : src.species_ratios) {
                        if (ratio <= 0 || sp_id >= static_cast<int>(species_mgr.getSpeciesCount())) continue;
                        const auto& sp = species_mgr.getSpecies(sp_id);
                        
                        double k1 = (1.0 - sp.k_ldd) * std::exp(-dist / sp.lambda1_m);
                        double k2 = sp.k_ldd * std::exp(-dist / sp.lambda2_m);
                        double kernel = k1 + k2;
                        
                        for (int y = 0; y < Plot::SIZE_M; ++y) {
                            for (int x = 0; x < Plot::SIZE_M; ++x) {
                                double light = plot.getCells().get(x, y).light_at_ground;
                                
                                bool ok = (light >= sp.l_min && light <= sp.l_max) &&
                                          (ngs_temp >= sp.ngs_tmin_c && ngs_temp < sp.ngs_tmax_c) &&
                                          (gdd >= sp.dd_min) && (drought <= sp.dr_tol);
                                
                                double p_est = ratio * kernel * 2.0;
                                if (ok && rng.getUniform01() < p_est) {
                                    Sapling sap(next_sap_id++, sp_id,
                                                x + rng.getUniform01(), y + rng.getUniform01(), 0.05);
                                    plot.getSaplings().push_back(sap);
                                }
                            }
                        }
                    }
                }
            }
            
            // Internal Seeds
            for (const auto& tree : plot.getTrees()) {
                if (!tree.isAlive()) continue;
                const auto& sp = species_mgr.getSpecies(tree.getSpeciesId());
                if (tree.getAge() < sp.maturity_age_yr) continue;
                
                double n_seeds = sp.fecundity_f * tree.getLeafArea() * 0.1;
                
                for (int i = 0; i < static_cast<int>(n_seeds); ++i) {
                    double angle = rng.getUniform01() * 2.0 * 3.14159265359;
                    double u = rng.getUniform01();
                    double dist = (u < (1.0 - sp.k_ldd)) ?
                        -sp.lambda1_m * std::log(1.0 - rng.getUniform01()) :
                        -sp.lambda2_m * std::log(1.0 - rng.getUniform01());
                    
                    int cx = static_cast<int>(tree.getX() + dist * std::cos(angle));
                    int cy = static_cast<int>(tree.getY() + dist * std::sin(angle));
                    
                    if (cx < 0 || cx >= Plot::SIZE_M || cy < 0 || cy >= Plot::SIZE_M) continue;
                    
                    double light = plot.getCells().get(cx, cy).light_at_ground;
                    bool ok = (light >= sp.l_min && light <= sp.l_max) &&
                              (ngs_temp >= sp.ngs_tmin_c && ngs_temp < sp.ngs_tmax_c) &&
                              (gdd >= sp.dd_min) && (drought <= sp.dr_tol);
                    
                    if (ok && rng.getUniform01() < 0.1) {
                        Sapling sap(next_sap_id++, tree.getSpeciesId(),
                                    cx + rng.getUniform01(), cy + rng.getUniform01(), 0.05);
                        plot.getSaplings().push_back(sap);
                    }
                }
            }
            
            // Growth, Recruitment, Mortality
            for (auto& sap : plot.getSaplings()) {
                if (!sap.isAlive()) continue;
                const auto& sp = species_mgr.getSpecies(sap.getSpeciesId());
                int cx = std::clamp(static_cast<int>(sap.getX()), 0, Plot::SIZE_M - 1);
                int cy = std::clamp(static_cast<int>(sap.getY()), 0, Plot::SIZE_M - 1);
                double light = light_engine.calcSaplingLight(sap, herb_lai.get(cx, cy));
                double f_env = EnvironmentResponse::calcLightResponse(light, sp.shade_tol) *
                               EnvironmentResponse::calcTemperatureResponse(gdd, sp.dd_min) *
                               EnvironmentResponse::calcDroughtResponse(drought, sp.dr_tol);
                sap.grow(SaplingGrowthCalculator::calcHeightIncrement(sap.getHeight(), sp.h_max_m, sp.gs, f_env));
                sap.incrementAge();
            }
            
            for (auto& tree : plot.getTrees()) {
                if (!tree.isAlive()) continue;
                const auto& sp = species_mgr.getSpecies(tree.getSpeciesId());
                double f_env = EnvironmentResponse::calcLightResponse(tree.getAvailableLight(), sp.shade_tol) *
                               EnvironmentResponse::calcTemperatureResponse(gdd, sp.dd_min) *
                               EnvironmentResponse::calcDroughtResponse(drought, sp.dr_tol);
                double ddbh = TreeGrowthCalculator::calcOptimalDbhIncrement(
                    tree.getDbh(), tree.getHeight(), sp.dbh_max_cm, sp.h_max_m, sp.ga, sp.opt_s) * f_env;
                tree.setDbh(tree.getDbh() + ddbh);
                tree.setHeight(TreeGrowthCalculator::calcHeightFromDbh(sp.h_max_m, sp.opt_s, tree.getDbh()));
                tree.incrementAge();
                tree.setCrownRadius(std::max(0.5, sp.r_a * std::log(std::max(0.1, tree.getDbh())) + sp.r_b));
                tree.setLeafArea(sp.kc2 * sp.ka1 * std::pow(std::max(0.1, tree.getDbh()), sp.ka2));
            }
            
            for (auto it = plot.getSaplings().begin(); it != plot.getSaplings().end(); ) {
                if (it->isAlive() && it->getHeight() >= 1.37) {
                    const auto& sp = species_mgr.getSpecies(it->getSpeciesId());
                    AdultTree tree(next_tree_id++, it->getSpeciesId(), it->getX(), it->getY(), 0.1, 1.37, it->getAge());
                    tree.setCrownRadius(std::max(0.5, sp.r_a * std::log(0.1) + sp.r_b));
                    tree.setCrownBase(sp.alpha_base * 1.37);
                    tree.setCrownShape(sp.cs);
                    tree.setLeafArea(sp.kc2 * sp.ka1 * std::pow(0.1, sp.ka2));
                    plot.getTrees().push_back(tree);
                    it = plot.getSaplings().erase(it);
                } else {
                    ++it;
                }
            }
            
            for (auto& sap : plot.getSaplings())
                if (sap.isAlive() && rng.getUniform01() < 0.02) sap.kill();
            
            for (auto& tree : plot.getTrees()) {
                if (tree.isAlive()) {
                    const auto& sp = species_mgr.getSpecies(tree.getSpeciesId());
                    if (rng.getUniform01() < (1.0 - std::pow(0.01, 1.0 / sp.max_age_yr))) tree.kill();
                }
            }
            
            plot.getSaplings().erase(
                std::remove_if(plot.getSaplings().begin(), plot.getSaplings().end(),
                    [](const Sapling& s) { return !s.isAlive(); }),
                plot.getSaplings().end());
            plot.getTrees().erase(
                std::remove_if(plot.getTrees().begin(), plot.getTrees().end(),
                    [](const AdultTree& t) { return !t.isAlive(); }),
                plot.getTrees().end());
        }
        
        if (year % SimConfig::OUTPUT_INTERVAL == 0 || year == SimConfig::END_YEAR || fire_this_year) {
            int tt = 0, ts = 0;
            double bio = 0;
            for (const auto& pl : plots) {
                tt += pl.getTrees().size();
                ts += pl.getSaplings().size();
                for (const auto& tr : pl.getTrees()) {
                    const auto& sp = species_mgr.getSpecies(tr.getSpeciesId());
                    bio += sp.b0 * std::pow(tr.getDbh(), sp.b1) * std::pow(tr.getHeight(), sp.b2);
                }
            }
            std::cout << "Year " << std::setw(4) << year
                      << " | Trees: " << std::setw(5) << tt
                      << " | Saplings: " << std::setw(5) << ts
                      << " | Biomass: " << std::fixed << std::setprecision(0) << bio << " kg\n";
        }
    }
    
    // ========================================================================
    // Step 6: Final Results
    // ========================================================================
    std::cout << "\n" << std::string(50, '=') << "\n";
    std::cout << "> Simulation Summary\n";
    std::cout << std::string(50, '=') << "\n\n";
    
    for (size_t p = 0; p < plots.size(); ++p) {
        std::cout << "[Plot] " << p << " (" << plot_names[p] << "):\n";
        std::cout << "   Adult trees: " << plots[p].getTrees().size() << "\n";
        std::cout << "   Saplings: " << plots[p].getSaplings().size() << "\n";
        
        std::map<int, int> sp_cnt;
        double total_bio = 0;
        for (const auto& tr : plots[p].getTrees()) {
            sp_cnt[tr.getSpeciesId()]++;
            const auto& sp = species_mgr.getSpecies(tr.getSpeciesId());
            total_bio += sp.b0 * std::pow(tr.getDbh(), sp.b1) * std::pow(tr.getHeight(), sp.b2);
        }
        
        std::cout << "   Total biomass: " << std::fixed << std::setprecision(0) << total_bio << " kg\n";
        std::cout << "   Species composition:\n";
        for (const auto& [sid, cnt] : sp_cnt) {
            std::cout << "     - " << species_mgr.getSpecies(sid).species_name << ": " << cnt << " trees\n";
        }
        std::cout << "\n";
    }
    
    auto t_end = std::chrono::high_resolution_clock::now();
    auto dur = std::chrono::duration_cast<std::chrono::seconds>(t_end - t_start);
    
    std::cout << "Total runtime: " << dur.count() << " seconds\n";
    std::cout << "\n+----------------------------------------------------------+\n";
    std::cout << "|                    Simulation Complete                   |\n";
    std::cout << "+----------------------------------------------------------+\n";
    
    return 0;
}