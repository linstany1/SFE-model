#ifndef DATA_OUTPUT_H
#define DATA_OUTPUT_H

#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <memory>
#include "Config.h"
#include "SpeciesParams.h"

class Plot;
struct WaterBalanceResult;
struct EstablishmentDebugInfo;
struct PlotSeedStats;  // NEW: Forward declaration

// ===========================================
// Annual Statistics for Plot Summary
// ===========================================
struct AnnualStats {
    // Climate indicators
    double GDD_deciduous;
    double GDD_evergreen;
    double T_winter;
    
    // Mortality statistics
    int trees_died;
    int saplings_died;
    double mortality_rate_pct;
    
    // Age statistics
    int tree_age_max;
    double tree_age_mean;
    int sapling_age_max;
    double sapling_age_mean;
    
    // Herb layer
    double herb_biomass_kg;
    double herb_lai_mean;
    
    // Establishment & Recruitment
    int saplings_established;
    int trees_recruited;
    
    // Stress statistics (currently under stress)
    int trees_stressed;
    int saplings_stressed;
    double f_env_tree_mean;
    double f_env_sapling_mean;
    
    AnnualStats();
};

// ===========================================
// Species-level Annual Statistics
// ===========================================
struct SpeciesAnnualStats {
    int trees_died;
    int saplings_died;
    int saplings_established;
    int tree_age_max;
    double tree_age_mean;
    int sapling_age_max;
    double sapling_age_mean;
    int trees_stressed;
    int saplings_stressed;
    double f_env_tree_mean;
    double f_env_sapling_mean;
    
    SpeciesAnnualStats();
};

class DataOutput {
private:
    std::string output_dir_;
    bool debug_mode_;
    
    std::ofstream plot_summary_file_;
    std::ofstream species_composition_file_;
    std::ofstream disturbance_log_file_;
    std::ofstream debug_seed_file_;
    
public:
    DataOutput(const std::string& output_dir, bool debug_mode = false);
    ~DataOutput();
    
    void initializeOutputFiles();
    void closeAllFiles();
    
    // UPDATED: Plot summary with extended diagnostics
    void writePlotSummary(int year_idx, const Plot& plot, 
                          const SpeciesParamsManager& sp,
                          const WaterBalanceResult& water_result,
                          const AnnualStats& stats);
    
    // UPDATED: Species composition with extended diagnostics
    void writeSpeciesComposition(int year_idx, const Plot& plot,
                                  const SpeciesParamsManager& sp,
                                  const std::map<int, int>& recruitment_counts,
                                  const std::map<int, int>& establishment_counts,
                                  const std::map<int, SpeciesAnnualStats>& species_stats);
    
    void writeIndividualTreeStatus(int year_idx, const Plot& plot,
                                    const SpeciesParamsManager& sp);
    
    void writeDisturbanceEvent(int year_idx, const Plot& plot,
                                const std::string& event_type,
                                double severity,
                                int trees_killed, int saplings_killed,
                                double biomass_removed);
    
    // UPDATED: Debug seed info with complete seed statistics
    void writeDebugSeedInfo(int year_idx, int plot_id,
                            const std::vector<EstablishmentDebugInfo>& debug_info,
                            const PlotSeedStats& seed_stats);
    
    void writeFinalSummary(const std::vector<std::unique_ptr<Plot>>& plots,
                           const SpeciesParamsManager& sp,
                           int total_years);
    
    bool isDebugMode() const;
    std::string getOutputDir() const;
};

#endif // DATA_OUTPUT_H