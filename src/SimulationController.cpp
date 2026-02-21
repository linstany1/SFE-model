#include "SimulationController.h"
#include <iostream>
#include <algorithm>
#include <filesystem>
#include <set>

namespace fs = std::filesystem;

// ===========================================
// SimulationConfig Implementation
// ===========================================

SimulationConfig::SimulationConfig()
    : spin_up_years(Config::SPIN_UP_YEARS),
      saturation_years(Config::SATURATION_YEARS),
      transient_start_year(1901),
      transient_end_year(2020),
      climate_base_window(Config::CLIMATE_BASE_WINDOW),
      output_interval(Config::DEFAULT_OUTPUT_INTERVAL),
      random_seed(42),
      debug_mode(false),
      output_dir("./output"),
      enable_herb_competition(true),
      enable_spinup_output(false),      // NEW: default to false (no spin-up output)
      spinup_output_interval(50) {}     // NEW: default interval 50 years

// ===========================================
// SimulationController Implementation
// ===========================================

SimulationController::SimulationController()
    : current_year_(0), current_calendar_year_(0), is_spinup_phase_(true) {}

bool SimulationController::initialize(const SimulationConfig& config,
                                       const std::string& species_file,
                                       const std::string& site_file,
                                       const std::string& climate_file,
                                       const std::string& distance_file,
                                       const std::string& events_file,
                                       const std::string& seed_source_file) {
    config_ = config;
    
    // Initialize RNG
    rng_.seed(config_.random_seed);
    
    std::cout << "Loading input files...\n";
    
    // Load species parameters
    if (!species_params_.loadFromFile(species_file)) {
        std::cerr << "Failed to load species parameters\n";
        return false;
    }
    std::cout << "  Loaded " << species_params_.getNumSpecies() << " species\n";
    
    // Get tree species IDs
    tree_species_ids_ = species_params_.getTreeSpeciesIds();
    
    // Load climate data
    if (!climate_manager_.loadFromFile(climate_file)) {
        std::cerr << "Failed to load climate data\n";
        return false;
    }
    std::cout << "  Loaded climate data for " << climate_manager_.getClimateIds().size() << " stations\n";
    
    // Load distance matrix
    if (!distance_matrix_.loadFromFile(distance_file)) {
        std::cerr << "Failed to load distance matrix\n";
        return false;
    }
    std::cout << "  Loaded distance matrix: " << distance_matrix_.getNumSources() 
              << " sources x " << distance_matrix_.getNumTargets() << " targets\n";
    
    // Load events
    event_manager_.loadFromFile(events_file);
    std::cout << "  Loaded " << event_manager_.getNumEvents() << " disturbance events\n";
    
    // Load site information and create plots
    auto sites = SiteLoader::loadFromFile(site_file);
    std::cout << "  Loaded " << sites.size() << " plots\n";
    
    for (const auto& site : sites) {
        auto plot = std::make_unique<Plot>(species_params_);
        plot->initialize(site);
        
        // NEW: Set herb competition mode for each plot's herb layer
        plot->herb_layer.setEnabled(config_.enable_herb_competition);
        
        plot_id_to_index_[site.plot_id] = plots_.size();
        plots_.push_back(std::move(plot));
    }
    
    // Initialize modules
    light_module_ = std::make_unique<LightModule>(species_params_);
    seed_module_ = std::make_unique<SeedDispersalModule>(species_params_, distance_matrix_);
    seed_module_->loadExternalSources(seed_source_file);
    estab_module_ = std::make_unique<EstablishmentModule>(species_params_);
    growth_module_ = std::make_unique<GrowthModule>(species_params_);
    mortality_module_ = std::make_unique<MortalityModule>(species_params_);
    
    // Initialize output
    data_output_ = std::make_unique<DataOutput>(config_.output_dir, config_.debug_mode);
    data_output_->initializeOutputFiles();
    
    // Calculate regional species pools for each plot
    std::vector<int> target_plot_ids;
    for (const auto& plot : plots_) {
        target_plot_ids.push_back(plot->site_info.plot_id);
    }
    
    std::cout << "  Calculating regional species pools...\n";
    for (const auto& plot : plots_) {
        regional_species_pools_[plot->site_info.plot_id] = 
            seed_module_->getRegionalSpeciesPool(plot->site_info.plot_id, target_plot_ids);
        
        std::cout << "    Plot " << plot->site_info.plot_id 
                  << " (" << plot->site_info.plot_name << "): " 
                  << regional_species_pools_[plot->site_info.plot_id].size() 
                  << " species in regional pool\n";
    }
    
    std::cout << "Initialization complete.\n\n";
    return true;
}

void SimulationController::run() {
    std::cout << "Starting simulation...\n";
    
    // Phase I: Spin-up
    runSpinup();
    
    // Phase II: Transient
    runTransient();
    
    // Final output
    data_output_->writeFinalSummary(plots_, species_params_, 
                                     config_.spin_up_years + 
                                     (config_.transient_end_year - config_.transient_start_year + 1));
    
    std::cout << "\nSimulation complete.\n";
}

// ===========================================
// Spin-up Phase
// ===========================================

void SimulationController::runSpinup() {
    std::cout << "Phase I: Spin-up (" << config_.spin_up_years << " years)\n";
    if (config_.enable_spinup_output) {
        std::cout << "  Spin-up output enabled (interval: " << config_.spinup_output_interval << " years)\n";
    }
    is_spinup_phase_ = true;
    
    for (int year = 1; year <= config_.spin_up_years; ++year) {
        current_year_ = year;
        
        if (year % 50 == 0) {
            std::cout << "  Spin-up year " << year << "/" << config_.spin_up_years << "\n";
        }
        
        // MODIFIED: Saturated mode disabled, external seed always enabled
        // Previously: saturated_mode = (year <= config_.saturation_years)
        // Now: Always use normal dispersal with external seed sources
        bool saturated_mode = false;
        bool external_seed_enabled = true;
        
        // Process each plot
        for (auto& plot : plots_) {
            // Get random climate year for spin-up
            const YearClimate& climate = climate_manager_.getRandomYearClimate(
                plot->site_info.climate_id, rng_, config_.climate_base_window);
            
            // Get previous year climate for NGS calculation
            const YearClimate& prev_climate = climate_manager_.getPreviousYearClimate(
                plot->site_info.climate_id, climate.year, rng_);
            
            // Run annual cycle with external seed dispersal enabled
            runAnnualCycle(*plot, climate, prev_climate, year, 0, saturated_mode, external_seed_enabled);
            
            // Output during spin-up phase if enabled
            if (config_.enable_spinup_output) {
                bool should_output = (year % config_.spinup_output_interval == 0) ||
                                    (year == config_.spin_up_years);
                
                if (should_output) {
                    outputAnnualData(*plot, year);
                }
            }
        }
    }
    
    std::cout << "  Spin-up complete.\n\n";
}

// ===========================================
// Transient Phase
// ===========================================

void SimulationController::runTransient() {
    int total_years = config_.transient_end_year - config_.transient_start_year + 1;
    std::cout << "Phase II: Transient simulation (" << total_years << " years)\n";
    is_spinup_phase_ = false;
    
    for (int calendar_year = config_.transient_start_year; 
         calendar_year <= config_.transient_end_year; ++calendar_year) {
        
        current_year_ = config_.spin_up_years + 
                       (calendar_year - config_.transient_start_year + 1);
        current_calendar_year_ = calendar_year;
        
        if ((calendar_year - config_.transient_start_year) % 10 == 0) {
            std::cout << "  Year " << calendar_year << "\n";
        }
        
        // Process each plot
        for (auto& plot : plots_) {
            // Check for disturbance events
            const DisturbanceEvent* event = event_manager_.getEvent(
                calendar_year, plot->site_info.plot_id);
            
            if (event && event->isFire()) {
                // Record pre-fire stats
                int trees_before = static_cast<int>(plot->trees.size());
                int saplings_before = static_cast<int>(plot->saplings.size());
                double biomass_before = plot->calcTotalBiomass(species_params_);
                
                // Apply fire
                plot->burn(calendar_year);
                
                // Log disturbance
                data_output_->writeDisturbanceEvent(current_year_, *plot,
                    "Fire", event->intensity, trees_before, saplings_before, biomass_before);
                
                std::cout << "    Fire at plot " << plot->site_info.plot_id << "\n";
            }
            
            // Get climate for this year
            const ClimateSeries& series = climate_manager_.getSeries(plot->site_info.climate_id);
            
            // Handle case where climate year might not exist
            const YearClimate* climate_ptr = nullptr;
            
            if (series.hasYear(calendar_year)) {
                climate_ptr = &series.getYear(calendar_year);
            } else {
                // Use random year from base window if calendar year not available
                climate_ptr = &climate_manager_.getRandomYearClimate(
                    plot->site_info.climate_id, rng_, config_.climate_base_window);
            }
            
            // Get previous year climate for NGS calculation
            const YearClimate& prev_climate = climate_manager_.getPreviousYearClimate(
                plot->site_info.climate_id, climate_ptr->year, rng_);
            
            // Run annual cycle
            // MODIFIED: External seed always enabled (not just during fire recovery)
            // Previously: external_seed_enabled = plot->is_in_fire_recovery
            // Now: Always true
            bool external_seed_enabled = true;
            runAnnualCycle(*plot, *climate_ptr, prev_climate,
                          current_year_, calendar_year, false, 
                          external_seed_enabled);
            
            // Output
            bool should_output = (current_year_ % config_.output_interval == 0) ||
                                (calendar_year == config_.transient_end_year);
            
            if (should_output) {
                outputAnnualData(*plot, current_year_);
            }
        }
    }
    
    std::cout << "  Transient simulation complete.\n";
}

// ===========================================
// Annual Simulation Cycle
// ===========================================

void SimulationController::runAnnualCycle(Plot& plot, const YearClimate& climate, 
                                           const YearClimate& prev_year_climate,
                                           int year_idx, int /*calendar_year*/,
                                           bool saturated_seed_mode, bool external_seed_enabled) {
    
    // ===========================================
    // Step 1: Rebuild spatial structures
    // ===========================================
    plot.rebuildSpatialHash();
    plot.rebuildCanopyMap();
    
    // ===========================================
    // Step 2: Calculate physical environment
    // ===========================================
    
    // 2a. Light calculation (top-down)
    light_module_->runLightCalculation(plot.trees, plot.saplings, 
                                        plot.herb_layer, plot.spatial_hash, plot.canopy_map);
    
    // 2b. Water balance (bottom-up) with leap year support
    WaterBalanceResult water_result = plot.water_module->calcAnnualWaterBalance(
        climate, plot.trees, plot.saplings, plot.herb_layer, climate.year);
    
    // 2c. GDD calculation (with leap year support via climate year)
    double GDD_deciduous = climate.calcGDD_Deciduous(climate.year);
    double GDD_evergreen = climate.calcGDD_Evergreen(climate.year);
    
    // 2d. Winter temperature calculation using previous year
    double T_winter = ClimateManager::calcWinterTemp(prev_year_climate, climate);
    
    // ===========================================
    // Step 3: Seed dispersal and establishment
    // ===========================================
    
    // 3a. Annual seed bank update
    plot.seed_bank.annualUpdate();
    
    // NEW: Record seed counts at different stages for debug output
    std::map<int, double> seeds_after_annual_update;
    std::map<int, double> seeds_after_local;
    std::map<int, double> seeds_after_external;
    
    if (config_.debug_mode) {
        seeds_after_annual_update = calcSeedTotalsBySpecies(plot.seed_bank, tree_species_ids_);
    }
    
    // 3b. Seed dispersal
    if (saturated_seed_mode) {
        // Mode C: Saturated seed rain (using regional species pool)
        const auto& regional_pool = regional_species_pools_[plot.site_info.plot_id];
        seed_module_->applySaturatedSeedRain(plot.seed_bank, regional_pool);
        
        if (config_.debug_mode) {
            seeds_after_local = calcSeedTotalsBySpecies(plot.seed_bank, tree_species_ids_);
            seeds_after_external = seeds_after_local;  // No external in saturated mode
        }
    } else {
        // Mode A: Internal dispersal (always on)
        seed_module_->calcWithinPlotDispersal(plot.trees, plot.seed_bank, 
                                               plot.site_info.plot_id);
        
        // Between-target dispersal
        for (const auto& other_plot : plots_) {
            if (other_plot->site_info.plot_id != plot.site_info.plot_id) {
                seed_module_->calcBetweenPlotDispersal(
                    other_plot->site_info.plot_id, other_plot->trees,
                    plot.site_info.plot_id, plot.seed_bank);
            }
        }
        
        if (config_.debug_mode) {
            seeds_after_local = calcSeedTotalsBySpecies(plot.seed_bank, tree_species_ids_);
        }
        
        // Mode B: External dispersal (always enabled now)
        if (external_seed_enabled) {
            seed_module_->calcExternalSourceDispersal(plot.site_info.plot_id, 
                                                      plot.seed_bank);
        }
        
        if (config_.debug_mode) {
            seeds_after_external = calcSeedTotalsBySpecies(plot.seed_bank, tree_species_ids_);
        }
    }
    
    // Record sapling count before establishment for statistics
    std::map<int, int> saplings_before;
    for (const auto& sap : plot.saplings) {
        saplings_before[sap.species_id]++;
    }
    
    // 3c. Establishment with NGS temperature check
    std::vector<EstablishmentDebugInfo> debug_info;
    
    // Use regional species pool in saturation mode, otherwise use global species list
    const auto& species_for_establishment = saturated_seed_mode ? 
        regional_species_pools_[plot.site_info.plot_id] : tree_species_ids_;
    
    // MODIFIED v2.0: Pass full water_result for per-species DI lookup
    estab_module_->processEstablishment(
        plot.seed_bank, plot.herb_layer,
        GDD_deciduous, GDD_evergreen, water_result,
        T_winter, species_for_establishment,
        plot.saplings, rng_, saturated_seed_mode,
        config_.debug_mode ? &debug_info : nullptr);
    
    // Calculate established saplings by species
    std::map<int, int> establishment_counts;
    for (const auto& sap : plot.saplings) {
        int count_before = saplings_before.count(sap.species_id) ? saplings_before[sap.species_id] : 0;
        int current_count = 0;
        for (const auto& s : plot.saplings) {
            if (s.species_id == sap.species_id) current_count++;
        }
        establishment_counts[sap.species_id] = current_count - count_before;
    }
    // Correct: count new saplings properly
    std::map<int, int> saplings_after;
    for (const auto& sap : plot.saplings) {
        saplings_after[sap.species_id]++;
    }
    establishment_counts.clear();
    for (const auto& pair : saplings_after) {
        int before = saplings_before.count(pair.first) ? saplings_before[pair.first] : 0;
        if (pair.second > before) {
            establishment_counts[pair.first] = pair.second - before;
        }
    }
    
    if (config_.debug_mode) {
        // Calculate PlotSeedStats from recorded seed counts
        PlotSeedStats seed_stats;
        double num_cells = Config::GRID_DIM * Config::GRID_DIM;
        
        for (int species_id : tree_species_ids_) {
            double seeds_base = seeds_after_annual_update.count(species_id) ? 
                               seeds_after_annual_update[species_id] : 0.0;
            double seeds_local_total = seeds_after_local.count(species_id) ? 
                                       seeds_after_local[species_id] : 0.0;
            double seeds_external_total = seeds_after_external.count(species_id) ? 
                                          seeds_after_external[species_id] : 0.0;
            
            SpeciesSeedStats ss;
            ss.seeds_local = seeds_local_total - seeds_base;  // Local contribution
            ss.seeds_external = seeds_external_total - seeds_local_total;  // External contribution
            ss.seeds_total = seeds_external_total;  // Total after all dispersal
            ss.seeds_mean_per_cell = seeds_external_total / num_cells;
            
            // Only include species with seeds
            if (ss.seeds_total > 0) {
                seed_stats.by_species[species_id] = ss;
                seed_stats.total_seeds_local += ss.seeds_local;
                seed_stats.total_seeds_external += ss.seeds_external;
                seed_stats.total_seeds_all_species += ss.seeds_total;
            }
        }
        
        data_output_->writeDebugSeedInfo(year_idx, plot.site_info.plot_id, debug_info, seed_stats);
    }
    
    // ===========================================
    // Step 4: Growth and dynamics
    // ===========================================
    
    // 4a. Herb growth Ã¢â‚¬â€ use herb-specific DI
    growth_module_->growHerbLayer(plot.herb_layer, GDD_deciduous, water_result.getSpeciesDI(0));
    
    // 4b. Sapling growth Ã¢â‚¬â€ MODIFIED v2.0: per-species DI
    for (auto& sapling : plot.saplings) {
        double DI_sp = water_result.getSpeciesDI(sapling.species_id);
        double GDD_sp = species_params_.getById(sapling.species_id).is_evergreen ?
                        GDD_evergreen : GDD_deciduous;
        growth_module_->growSapling(sapling, GDD_sp, DI_sp);
    }
    
    // 4c. Adult tree growth Ã¢â‚¬â€ MODIFIED v2.0: per-species DI
    for (auto& tree : plot.trees) {
        double DI_sp = water_result.getSpeciesDI(tree.species_id);
        double GDD_sp = species_params_.getById(tree.species_id).is_evergreen ?
                        GDD_evergreen : GDD_deciduous;
        growth_module_->growAdultTree(tree, GDD_sp, DI_sp);
    }
    
    // 4d. Recruitment (sapling -> adult)
    std::map<int, int> recruitment_counts;
    processRecruitment(plot, recruitment_counts);
    
    // 4e. Mortality (with new formulas) - now returns death counts
    std::map<int, int> trees_died_by_species;
    std::map<int, int> saplings_died_by_species;
    processMortality(plot, trees_died_by_species, saplings_died_by_species);
    
    // ===========================================
    // Step 5: Collect statistics
    // ===========================================
    
    // Collect species-level statistics
    std::map<int, SpeciesAnnualStats> species_stats;
    collectStatistics(plot, species_stats);
    
    // Add mortality counts to species stats
    for (const auto& pair : trees_died_by_species) {
        species_stats[pair.first].trees_died = pair.second;
    }
    for (const auto& pair : saplings_died_by_species) {
        species_stats[pair.first].saplings_died = pair.second;
    }
    // Add establishment counts to species stats
    for (const auto& pair : establishment_counts) {
        species_stats[pair.first].saplings_established = pair.second;
    }
    
    // Calculate plot-level statistics
    AnnualStats stats;
    
    // Climate indicators
    stats.GDD_deciduous = GDD_deciduous;
    stats.GDD_evergreen = GDD_evergreen;
    stats.T_winter = T_winter;
    
    // Mortality statistics
    int total_trees_died = 0;
    int total_saplings_died = 0;
    for (const auto& pair : trees_died_by_species) total_trees_died += pair.second;
    for (const auto& pair : saplings_died_by_species) total_saplings_died += pair.second;
    stats.trees_died = total_trees_died;
    stats.saplings_died = total_saplings_died;
    
    // Mortality rate (trees that died / trees at start of mortality phase)
    int trees_before_mortality = static_cast<int>(plot.trees.size()) + total_trees_died;
    stats.mortality_rate_pct = (trees_before_mortality > 0) ? 
        (100.0 * total_trees_died / trees_before_mortality) : 0.0;
    
    // Age statistics
    int tree_age_max = 0;
    double tree_age_sum = 0;
    for (const auto& tree : plot.trees) {
        tree_age_max = std::max(tree_age_max, tree.age);
        tree_age_sum += tree.age;
    }
    stats.tree_age_max = tree_age_max;
    stats.tree_age_mean = plot.trees.empty() ? 0 : tree_age_sum / plot.trees.size();
    
    int sapling_age_max = 0;
    double sapling_age_sum = 0;
    for (const auto& sap : plot.saplings) {
        sapling_age_max = std::max(sapling_age_max, sap.age);
        sapling_age_sum += sap.age;
    }
    stats.sapling_age_max = sapling_age_max;
    stats.sapling_age_mean = plot.saplings.empty() ? 0 : sapling_age_sum / plot.saplings.size();
    
    // Herb layer
    stats.herb_biomass_kg = plot.herb_layer.getTotalBiomass();
    stats.herb_lai_mean = plot.herb_layer.getAverageLAI();
    
    // Establishment & Recruitment
    int total_established = 0;
    int total_recruited = 0;
    for (const auto& pair : establishment_counts) total_established += pair.second;
    for (const auto& pair : recruitment_counts) total_recruited += pair.second;
    stats.saplings_established = total_established;
    stats.trees_recruited = total_recruited;
    
    // Stress statistics
    int trees_stressed = 0;
    double f_env_tree_sum = 0;
    for (const auto& tree : plot.trees) {
        const SpeciesProfile& sp = species_params_.getById(tree.species_id);
        if (tree.f_env < sp.stress_threshold) trees_stressed++;
        f_env_tree_sum += tree.f_env;
    }
    stats.trees_stressed = trees_stressed;
    stats.f_env_tree_mean = plot.trees.empty() ? 0 : f_env_tree_sum / plot.trees.size();
    
    int saplings_stressed = 0;
    double f_env_sapling_sum = 0;
    for (const auto& sap : plot.saplings) {
        const SpeciesProfile& sp = species_params_.getById(sap.species_id);
        if (sap.f_env < sp.stress_threshold) saplings_stressed++;
        f_env_sapling_sum += sap.f_env;
    }
    stats.saplings_stressed = saplings_stressed;
    stats.f_env_sapling_mean = plot.saplings.empty() ? 0 : f_env_sapling_sum / plot.saplings.size();
    
    // Store results for output
    last_water_result_ = water_result;
    last_recruitment_counts_ = recruitment_counts;
    last_annual_stats_ = stats;
    last_establishment_counts_ = establishment_counts;
    last_species_stats_ = species_stats;
}

// ===========================================
// Recruitment Processing
// ===========================================

void SimulationController::processRecruitment(Plot& plot, std::map<int, int>& recruitment_counts) {
    // Find saplings ready for recruitment
    auto it = plot.saplings.begin();
    while (it != plot.saplings.end()) {
        if (it->checkRecruitment()) {
            // Convert to adult with random variation
            AdultTree new_tree = growth_module_->convertSaplingToAdult(
                *it, plot.getNextTreeId(), rng_);
            
            // Add to trees (std::list - no invalidation)
            plot.trees.push_back(new_tree);
            
            // Incremental spatial hash update
            plot.insertTreeToHash(&plot.trees.back());
            
            // Track recruitment
            recruitment_counts[it->species_id]++;
            
            // Remove sapling (std::list - iterator safe erase)
            it = plot.saplings.erase(it);
        } else {
            ++it;
        }
    }
}

// ===========================================
// Mortality Processing
// ===========================================

void SimulationController::processMortality(Plot& plot, 
                                            std::map<int, int>& trees_died_by_species,
                                            std::map<int, int>& saplings_died_by_species) {
    trees_died_by_species.clear();
    saplings_died_by_species.clear();
    
    // Adult tree mortality with new formulas
    auto tree_it = plot.trees.begin();
    while (tree_it != plot.trees.end()) {
        if (mortality_module_->checkTreeMortality(*tree_it, rng_)) {
            // Track death by species
            trees_died_by_species[tree_it->species_id]++;
            // Remove from spatial hash first
            plot.removeTreeFromHash(tree_it->id);
            // Remove from list (std::list - iterator safe erase)
            tree_it = plot.trees.erase(tree_it);
        } else {
            ++tree_it;
        }
    }
    
    // [Phase 1d v1] Sapling mortality: background + stress (paralleling adult logic)
    auto sap_it = plot.saplings.begin();
    while (sap_it != plot.saplings.end()) {
        const SpeciesProfile& sp = species_params_.getById(sap_it->species_id);
        if (mortality_module_->checkSaplingMortality(*sap_it, sp.stress_threshold,
                                                      sp.max_age_yr, rng_)) {
            // Track death by species
            saplings_died_by_species[sap_it->species_id]++;
            // Remove from list (std::list - iterator safe erase)
            sap_it = plot.saplings.erase(sap_it);
        } else {
            ++sap_it;
        }
    }
}

// ===========================================
// Statistics Collection
// ===========================================

void SimulationController::collectStatistics(Plot& plot, 
                                             std::map<int, SpeciesAnnualStats>& species_stats) {
    species_stats.clear();
    
    // Collect tree statistics by species
    std::map<int, std::vector<int>> tree_ages;
    std::map<int, std::vector<double>> tree_f_envs;
    std::map<int, int> trees_stressed_count;
    
    for (const auto& tree : plot.trees) {
        int sp_id = tree.species_id;
        tree_ages[sp_id].push_back(tree.age);
        tree_f_envs[sp_id].push_back(tree.f_env);
        
        const SpeciesProfile& sp = species_params_.getById(sp_id);
        if (tree.f_env < sp.stress_threshold) {
            trees_stressed_count[sp_id]++;
        }
    }
    
    // Collect sapling statistics by species
    std::map<int, std::vector<int>> sapling_ages;
    std::map<int, std::vector<double>> sapling_f_envs;
    std::map<int, int> saplings_stressed_count;
    
    for (const auto& sap : plot.saplings) {
        int sp_id = sap.species_id;
        sapling_ages[sp_id].push_back(sap.age);
        sapling_f_envs[sp_id].push_back(sap.f_env);
        
        const SpeciesProfile& sp = species_params_.getById(sp_id);
        if (sap.f_env < sp.stress_threshold) {
            saplings_stressed_count[sp_id]++;
        }
    }
    
    // Combine all species
    std::set<int> all_species;
    for (const auto& pair : tree_ages) all_species.insert(pair.first);
    for (const auto& pair : sapling_ages) all_species.insert(pair.first);
    
    // Calculate statistics for each species
    for (int sp_id : all_species) {
        SpeciesAnnualStats stats;
        
        // Tree age statistics
        if (tree_ages.count(sp_id) && !tree_ages[sp_id].empty()) {
            const auto& ages = tree_ages[sp_id];
            stats.tree_age_max = *std::max_element(ages.begin(), ages.end());
            double sum = 0;
            for (int a : ages) sum += a;
            stats.tree_age_mean = sum / ages.size();
        }
        
        // Tree f_env statistics
        if (tree_f_envs.count(sp_id) && !tree_f_envs[sp_id].empty()) {
            const auto& f_envs = tree_f_envs[sp_id];
            double sum = 0;
            for (double f : f_envs) sum += f;
            stats.f_env_tree_mean = sum / f_envs.size();
        }
        
        // Trees stressed
        stats.trees_stressed = trees_stressed_count.count(sp_id) ? trees_stressed_count[sp_id] : 0;
        
        // Sapling age statistics
        if (sapling_ages.count(sp_id) && !sapling_ages[sp_id].empty()) {
            const auto& ages = sapling_ages[sp_id];
            stats.sapling_age_max = *std::max_element(ages.begin(), ages.end());
            double sum = 0;
            for (int a : ages) sum += a;
            stats.sapling_age_mean = sum / ages.size();
        }
        
        // Sapling f_env statistics
        if (sapling_f_envs.count(sp_id) && !sapling_f_envs[sp_id].empty()) {
            const auto& f_envs = sapling_f_envs[sp_id];
            double sum = 0;
            for (double f : f_envs) sum += f;
            stats.f_env_sapling_mean = sum / f_envs.size();
        }
        
        // Saplings stressed
        stats.saplings_stressed = saplings_stressed_count.count(sp_id) ? saplings_stressed_count[sp_id] : 0;
        
        species_stats[sp_id] = stats;
    }
}

// ===========================================
// Output
// ===========================================

void SimulationController::outputAnnualData(Plot& plot, int year_idx) {
    data_output_->writePlotSummary(year_idx, plot, species_params_, 
                                    last_water_result_, last_annual_stats_);
    data_output_->writeSpeciesComposition(year_idx, plot, species_params_, 
                                           last_recruitment_counts_,
                                           last_establishment_counts_,
                                           last_species_stats_);
    data_output_->writeIndividualTreeStatus(year_idx, plot, species_params_);
}

// ===========================================
// Seed Statistics Helper
// ===========================================

std::map<int, double> SimulationController::calcSeedTotalsBySpecies(
    const SeedBank& seed_bank, 
    const std::vector<int>& species_ids) const {
    
    std::map<int, double> totals;
    
    for (int species_id : species_ids) {
        double total = 0.0;
        for (int u = 0; u < Config::GRID_DIM; ++u) {
            for (int v = 0; v < Config::GRID_DIM; ++v) {
                total += seed_bank.getSeeds(u, v, species_id);
            }
        }
        totals[species_id] = total;
    }
    
    return totals;
}

// ===========================================
// Accessors
// ===========================================

const std::vector<std::unique_ptr<Plot>>& SimulationController::getPlots() const { 
    return plots_; 
}

const SpeciesParamsManager& SimulationController::getSpeciesParams() const { 
    return species_params_; 
}

int SimulationController::getCurrentYear() const { 
    return current_year_; 
}