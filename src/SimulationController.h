#ifndef SIMULATION_CONTROLLER_H
#define SIMULATION_CONTROLLER_H

#include <vector>
#include <memory>
#include <random>
#include <string>
#include <map>
#include <list>
#include "Config.h"
#include "Plot.h"
#include "SpeciesParams.h"
#include "ClimateManager.h"
#include "DistanceMatrix.h"
#include "EventManager.h"
#include "LightModule.h"
#include "WaterModule.h"
#include "SeedDispersalModule.h"
#include "EstablishmentModule.h"
#include "GrowthModule.h"
#include "MortalityModule.h"
#include "DataOutput.h"

// Simulation configuration
struct SimulationConfig {
    int spin_up_years;
    int saturation_years;
    int transient_start_year;
    int transient_end_year;
    int climate_base_window;
    int output_interval;
    unsigned int random_seed;
    bool debug_mode;
    std::string output_dir;
    
    // Herb competition control
    // When true (default), herb layer affects sapling light and water competition
    // When false, herb layer grows but does not compete with trees/saplings
    bool enable_herb_competition;
    
    // NEW: Spin-up output control
    // When true, output data during spin-up phase at specified interval
    bool enable_spinup_output;
    int spinup_output_interval;  // Output interval for spin-up phase (default: 50)
    
    SimulationConfig();
};

class SimulationController {
private:
    // Data Managers
    SpeciesParamsManager species_params_;
    ClimateManager climate_manager_;
    DistanceMatrix distance_matrix_;
    EventManager event_manager_;
    
    // Plots
    std::vector<std::unique_ptr<Plot>> plots_;
    std::map<int, size_t> plot_id_to_index_;
    
    // Modules
    std::unique_ptr<LightModule> light_module_;
    std::unique_ptr<SeedDispersalModule> seed_module_;
    std::unique_ptr<EstablishmentModule> estab_module_;
    std::unique_ptr<GrowthModule> growth_module_;
    std::unique_ptr<MortalityModule> mortality_module_;
    std::unique_ptr<DataOutput> data_output_;
    
    // RNG
    std::mt19937 rng_;
    
    // Configuration
    SimulationConfig config_;
    
    // State Tracking
    int current_year_;
    int current_calendar_year_;
    bool is_spinup_phase_;
    std::vector<int> tree_species_ids_;
    
    // Regional species pools for each plot (used in saturation phase)
    std::map<int, std::vector<int>> regional_species_pools_;
    
    // Output state
    WaterBalanceResult last_water_result_;
    std::map<int, int> last_recruitment_counts_;
    
    // NEW: Extended diagnostic statistics
    AnnualStats last_annual_stats_;
    std::map<int, int> last_establishment_counts_;
    std::map<int, SpeciesAnnualStats> last_species_stats_;
    
    // Private methods
    void runSpinup();
    void runTransient();
    void runAnnualCycle(Plot& plot, const YearClimate& climate, 
                        const YearClimate& prev_year_climate,
                        int year_idx, int /*calendar_year*/,
                        bool saturated_seed_mode, bool external_seed_enabled);
    void processRecruitment(Plot& plot, std::map<int, int>& recruitment_counts);
    void processMortality(Plot& plot, std::map<int, int>& trees_died_by_species,
                          std::map<int, int>& saplings_died_by_species);
    void collectStatistics(Plot& plot, std::map<int, SpeciesAnnualStats>& species_stats);
    void outputAnnualData(Plot& plot, int year_idx);
    
    // NEW: Calculate seed statistics from SeedBank
    std::map<int, double> calcSeedTotalsBySpecies(const SeedBank& seed_bank, 
                                                   const std::vector<int>& species_ids) const;
    
public:
    SimulationController();
    ~SimulationController() = default;
    
    // Initialization
    bool initialize(const SimulationConfig& config,
                    const std::string& species_file,
                    const std::string& site_file,
                    const std::string& climate_file,
                    const std::string& distance_file,
                    const std::string& events_file,
                    const std::string& seed_source_file);
    
    // Main simulation
    void run();
    
    // Accessors
    const std::vector<std::unique_ptr<Plot>>& getPlots() const;
    const SpeciesParamsManager& getSpeciesParams() const;
    int getCurrentYear() const;
};

#endif // SIMULATION_CONTROLLER_H