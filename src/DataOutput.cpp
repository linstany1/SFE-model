#include "DataOutput.h"
#include "Plot.h"
#include "WaterModule.h"
#include "EstablishmentModule.h"
#include <iomanip>
#include <set>

// ===========================================
// AnnualStats Constructor
// ===========================================
AnnualStats::AnnualStats()
    : GDD_deciduous(0), GDD_evergreen(0), T_winter(0),
      trees_died(0), saplings_died(0), mortality_rate_pct(0),
      tree_age_max(0), tree_age_mean(0), sapling_age_max(0), sapling_age_mean(0),
      herb_biomass_kg(0), herb_lai_mean(0),
      saplings_established(0), trees_recruited(0),
      trees_stressed(0), saplings_stressed(0),
      f_env_tree_mean(0), f_env_sapling_mean(0) {}

// ===========================================
// SpeciesAnnualStats Constructor
// ===========================================
SpeciesAnnualStats::SpeciesAnnualStats()
    : trees_died(0), saplings_died(0), saplings_established(0),
      tree_age_max(0), tree_age_mean(0), sapling_age_max(0), sapling_age_mean(0),
      trees_stressed(0), saplings_stressed(0),
      f_env_tree_mean(0), f_env_sapling_mean(0) {}

DataOutput::DataOutput(const std::string& output_dir, bool debug_mode)
    : output_dir_(output_dir), debug_mode_(debug_mode) {}

DataOutput::~DataOutput() {
    closeAllFiles();
}

void DataOutput::initializeOutputFiles() {
    // Plot year summary - UPDATED: with extended diagnostics
    plot_summary_file_.open(output_dir_ + "/plot_year_summary.txt");
    plot_summary_file_ << "year_idx,plot_id,plot_name,region_group,"
                       << "lai_total,cc_total,biomass_ag_t_ha,adult_density_ha,sapling_density_ha,"
                       << "sw_top_mean,sw_sub_mean,al_floor_mean,di_annual_mean,"
                       << "GDD_deciduous,GDD_evergreen,T_winter,"
                       << "trees_died,saplings_died,mortality_rate_pct,"
                       << "tree_age_max,tree_age_mean,sapling_age_max,sapling_age_mean,"
                       << "herb_biomass_kg,herb_lai_mean,"
                       << "saplings_established,trees_recruited,"
                       << "trees_stressed,saplings_stressed,f_env_tree_mean,f_env_sapling_mean\n";
    
    // Species composition - UPDATED: with extended diagnostics
    species_composition_file_.open(output_dir_ + "/species_composition_plot.txt");
    species_composition_file_ << "year_idx,plot_id,plot_name,species_id,"
                              << "adult_count,sapling_count,lai,cc,"
                              << "biomass_sum_kg,mean_dbh_cm,mean_height_m,"
                              << "recruitment_success,saplings_established,"
                              << "trees_died,saplings_died,"
                              << "tree_age_max,tree_age_mean,sapling_age_max,sapling_age_mean,"
                              << "trees_stressed,saplings_stressed,f_env_tree_mean,f_env_sapling_mean\n";
    
    // Disturbance log
    disturbance_log_file_.open(output_dir_ + "/disturbance_log.txt");
    disturbance_log_file_ << "year_idx,plot_id,plot_name,disturbance_type,"
                           << "severity,trees_killed,saplings_killed,biomass_removed\n";
    
    if (debug_mode_) {
        debug_seed_file_.open(output_dir_ + "/debug_seed_regeneration.txt");
        debug_seed_file_ << "year_idx,plot_id,species_id,"
                         << "seeds_local,seeds_external,seeds_total,seeds_mean_per_cell,"
                         << "p_seed_mean,p_env_mean,p_est_mean,"
                         << "n_new_saplings,limiting_factor\n";
    }
}

void DataOutput::closeAllFiles() {
    if (plot_summary_file_.is_open()) plot_summary_file_.close();
    if (species_composition_file_.is_open()) species_composition_file_.close();
    if (disturbance_log_file_.is_open()) disturbance_log_file_.close();
    if (debug_seed_file_.is_open()) debug_seed_file_.close();
}

void DataOutput::writePlotSummary(int year_idx, const Plot& plot, 
                          const SpeciesParamsManager& sp,
                          const WaterBalanceResult& water_result,
                          const AnnualStats& stats) {
    if (!plot_summary_file_.is_open()) return;
    
    plot_summary_file_ << std::fixed << std::setprecision(3);
    plot_summary_file_ << year_idx << ","
                       << plot.site_info.plot_id << ","
                       << "\"" << plot.site_info.plot_name << "\","
                       << "\"" << plot.site_info.region_group << "\","
                       // Basic metrics
                       << plot.calcTotalLAI(sp) << ","
                       << plot.calcTotalCrownCover() << ","
                       << plot.calcBiomassPerHa(sp) << ","
                       << plot.calcTreeDensityPerHa() << ","
                       << plot.calcSaplingDensityPerHa() << ","
                       // Water balance
                       << water_result.water_top_mean << ","
                       << water_result.water_sub_mean << ","
                       << plot.getMeanFloorLight() << ","
                       << water_result.DI_annual_mean << ","
                       // Climate indicators
                       << stats.GDD_deciduous << ","
                       << stats.GDD_evergreen << ","
                       << stats.T_winter << ","
                       // Mortality
                       << stats.trees_died << ","
                       << stats.saplings_died << ","
                       << stats.mortality_rate_pct << ","
                       // Age statistics
                       << stats.tree_age_max << ","
                       << stats.tree_age_mean << ","
                       << stats.sapling_age_max << ","
                       << stats.sapling_age_mean << ","
                       // Herb layer
                       << stats.herb_biomass_kg << ","
                       << stats.herb_lai_mean << ","
                       // Establishment & Recruitment
                       << stats.saplings_established << ","
                       << stats.trees_recruited << ","
                       // Stress statistics
                       << stats.trees_stressed << ","
                       << stats.saplings_stressed << ","
                       << stats.f_env_tree_mean << ","
                       << stats.f_env_sapling_mean << "\n";
    
    plot_summary_file_.flush();
}

void DataOutput::writeSpeciesComposition(int year_idx, const Plot& plot,
                                  const SpeciesParamsManager& sp,
                                  const std::map<int, int>& recruitment_counts,
                                  const std::map<int, int>& establishment_counts,
                                  const std::map<int, SpeciesAnnualStats>& species_stats) {
    if (!species_composition_file_.is_open()) return;
    
    auto tree_counts = plot.countTreesBySpecies();
    auto sapling_counts = plot.countSaplingsBySpecies();
    auto leaf_areas = plot.calcLeafAreaBySpecies();
    auto crown_covers = plot.calcCrownCoverBySpecies();
    auto mean_dbh = plot.calcMeanDBHBySpecies();
    auto mean_height = plot.calcMeanHeightBySpecies();
    
    std::set<int> all_species;
    for (const auto& p : tree_counts) all_species.insert(p.first);
    for (const auto& p : sapling_counts) all_species.insert(p.first);
    for (const auto& p : species_stats) all_species.insert(p.first);
    
    double plot_area = Config::PLOT_SIZE * Config::PLOT_SIZE;
    
    species_composition_file_ << std::fixed << std::setprecision(3);
    
    for (int species_id : all_species) {
        int adult_count = tree_counts.count(species_id) ? tree_counts[species_id] : 0;
        int sapling_count = sapling_counts.count(species_id) ? sapling_counts[species_id] : 0;
        double lai = leaf_areas.count(species_id) ? leaf_areas[species_id] / plot_area : 0;
        double cc = crown_covers.count(species_id) ? crown_covers[species_id] : 0;
        double dbh = mean_dbh.count(species_id) ? mean_dbh[species_id] : 0;
        double h = mean_height.count(species_id) ? mean_height[species_id] : 0;
        int recruit = recruitment_counts.count(species_id) ? recruitment_counts.at(species_id) : 0;
        int established = establishment_counts.count(species_id) ? establishment_counts.at(species_id) : 0;
        
        // Get species stats (use defaults if not present)
        SpeciesAnnualStats stats;
        if (species_stats.count(species_id)) {
            stats = species_stats.at(species_id);
        }
        
        double biomass = 0;
        for (const auto& tree : plot.trees) {
            if (tree.species_id == species_id) {
                biomass += tree.calcBiomass(sp.getById(species_id));
            }
        }
        
        species_composition_file_ << year_idx << ","
                                  << plot.site_info.plot_id << ","
                                  << "\"" << plot.site_info.plot_name << "\","
                                  << species_id << ","
                                  // Basic counts
                                  << adult_count << ","
                                  << sapling_count << ","
                                  << lai << ","
                                  << cc << ","
                                  << biomass << ","
                                  << dbh << ","
                                  << h << ","
                                  // Recruitment & Establishment
                                  << recruit << ","
                                  << established << ","
                                  // Mortality
                                  << stats.trees_died << ","
                                  << stats.saplings_died << ","
                                  // Age statistics
                                  << stats.tree_age_max << ","
                                  << stats.tree_age_mean << ","
                                  << stats.sapling_age_max << ","
                                  << stats.sapling_age_mean << ","
                                  // Stress statistics
                                  << stats.trees_stressed << ","
                                  << stats.saplings_stressed << ","
                                  << stats.f_env_tree_mean << ","
                                  << stats.f_env_sapling_mean << "\n";
    }
    
    species_composition_file_.flush();
}

void DataOutput::writeIndividualTreeStatus(int year_idx, const Plot& plot,
                                    const SpeciesParamsManager& sp) {
    std::string filename = output_dir_ + "/individual_adulttree_status.txt";
    
    static bool header_written = false;
    std::ofstream file;
    
    if (!header_written) {
        file.open(filename);
        file << "year_idx,plot_id,plot_name,tree_id,species_id,"
             << "x_coord,y_coord,age,dbh_cm,height_m,"
             << "biomass_ag_kg,crown_radius_m\n";
        header_written = true;
    } else {
        file.open(filename, std::ios::app);
    }
    
    file << std::fixed << std::setprecision(3);
    
    for (const auto& tree : plot.trees) {
        file << year_idx << ","
             << plot.site_info.plot_id << ","
             << "\"" << plot.site_info.plot_name << "\","
             << tree.id << ","
             << tree.species_id << ","
             << tree.x << ","
             << tree.y << ","
             << tree.age << ","
             << tree.dbh << ","
             << tree.height << ","
             << tree.calcBiomass(sp.getById(tree.species_id)) << ","
             << tree.crown_radius << "\n";
    }
    
    file.close();
}

void DataOutput::writeDisturbanceEvent(int year_idx, const Plot& plot,
                                const std::string& event_type,
                                double severity,
                                int trees_killed, int saplings_killed,
                                double biomass_removed) {
    if (!disturbance_log_file_.is_open()) return;
    
    disturbance_log_file_ << std::fixed << std::setprecision(3);
    disturbance_log_file_ << year_idx << ","
                          << plot.site_info.plot_id << ","
                          << "\"" << plot.site_info.plot_name << "\","
                          << "\"" << event_type << "\","
                          << severity << ","
                          << trees_killed << ","
                          << saplings_killed << ","
                          << biomass_removed << "\n";
    
    disturbance_log_file_.flush();
}

void DataOutput::writeDebugSeedInfo(int year_idx, int plot_id,
                            const std::vector<EstablishmentDebugInfo>& debug_info,
                            const PlotSeedStats& seed_stats) {
    if (!debug_mode_ || !debug_seed_file_.is_open()) return;
    
    // Aggregate debug_info by species
    std::map<int, double> p_seed_sum;
    std::map<int, double> p_env_sum;
    std::map<int, double> p_est_sum;
    std::map<int, int> count;
    std::map<int, int> n_established;
    std::map<int, std::string> limiting_factor;
    
    int total_established = 0;
    
    for (const auto& info : debug_info) {
        p_seed_sum[info.species_id] += info.P_seed;
        p_env_sum[info.species_id] += info.P_env;
        p_est_sum[info.species_id] += info.P_est;
        count[info.species_id]++;
        n_established[info.species_id] += info.n_established;
        total_established += info.n_established;
        // Keep the most common limiting factor (simplified: just take the last non-None)
        if (!info.limiting_factor.empty() && info.limiting_factor != "None") {
            limiting_factor[info.species_id] = info.limiting_factor;
        }
    }
    
    debug_seed_file_ << std::fixed << std::setprecision(4);
    
    // Output per-species statistics
    for (const auto& pair : seed_stats.by_species) {
        int sp_id = pair.first;
        const SpeciesSeedStats& ss = pair.second;
        
        double p_seed_mean = count.count(sp_id) && count[sp_id] > 0 ? 
                             p_seed_sum[sp_id] / count[sp_id] : 0.0;
        double p_env_mean = count.count(sp_id) && count[sp_id] > 0 ? 
                            p_env_sum[sp_id] / count[sp_id] : 0.0;
        double p_est_mean = count.count(sp_id) && count[sp_id] > 0 ? 
                            p_est_sum[sp_id] / count[sp_id] : 0.0;
        int n_est = n_established.count(sp_id) ? n_established[sp_id] : 0;
        std::string lim_factor = limiting_factor.count(sp_id) ? limiting_factor[sp_id] : "None";
        
        debug_seed_file_ << year_idx << ","
                         << plot_id << ","
                         << sp_id << ","
                         << ss.seeds_local << ","
                         << ss.seeds_external << ","
                         << ss.seeds_total << ","
                         << ss.seeds_mean_per_cell << ","
                         << p_seed_mean << ","
                         << p_env_mean << ","
                         << p_est_mean << ","
                         << n_est << ","
                         << "\"" << lim_factor << "\"\n";
    }
    
    // Output plot-level total (species_id = -1 indicates TOTAL row)
    double num_cells = Config::GRID_DIM * Config::GRID_DIM;
    debug_seed_file_ << year_idx << ","
                     << plot_id << ","
                     << -1 << ","  // species_id = -1 for TOTAL
                     << seed_stats.total_seeds_local << ","
                     << seed_stats.total_seeds_external << ","
                     << seed_stats.total_seeds_all_species << ","
                     << seed_stats.total_seeds_all_species / num_cells << ","
                     << 0.0 << ","  // p_seed_mean N/A for total
                     << 0.0 << ","  // p_env_mean N/A for total
                     << 0.0 << ","  // p_est_mean N/A for total
                     << total_established << ","
                     << "\"PLOT_TOTAL\"\n";
    
    debug_seed_file_.flush();
}

void DataOutput::writeFinalSummary(const std::vector<std::unique_ptr<Plot>>& plots,
                           const SpeciesParamsManager& sp,
                           int total_years) {
    std::string filename = output_dir_ + "/simulation_summary.txt";
    std::ofstream file(filename);
    
    file << "SFE-Model Simulation Summary\n";
    file << "============================\n\n";
    file << "Total simulation years: " << total_years << "\n";
    file << "Number of plots: " << plots.size() << "\n\n";
    
    for (const auto& plot : plots) {
        file << "Plot " << plot->site_info.plot_id << " (" 
             << plot->site_info.plot_name << "):\n";
        file << "  Adult trees: " << plot->trees.size() << "\n";
        file << "  Saplings: " << plot->saplings.size() << "\n";
        file << "  LAI: " << std::fixed << std::setprecision(2) 
             << plot->calcTotalLAI(sp) << "\n";
        file << "  Crown cover: " << plot->calcTotalCrownCover() << "\n";
        file << "  Fire recovery: " << (plot->is_in_fire_recovery ? "Yes" : "No") << "\n";
        file << "\n";
    }
    
    file.close();
}

bool DataOutput::isDebugMode() const { return debug_mode_; }
std::string DataOutput::getOutputDir() const { return output_dir_; }