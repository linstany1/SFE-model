/**
 * @file SimulationController.cpp
 * @brief 模拟控制器实现
 * 
 * 亚高山森林演替模型 (SFE-Model)
 * Phase 8: 两阶段模拟策略 (Spin-up + Transient)
 */

#include "simulation/SimulationController.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <chrono>

namespace sfe {

// ============================================================================
// 初始化
// ============================================================================

bool SimulationController::initialize(const SimulationConfig& config) {
    config_ = config;
    all_stats_.clear();
    
    // 初始化随机数生成器
    auto seed = std::chrono::steady_clock::now().time_since_epoch().count();
    rng_.seed(static_cast<unsigned int>(seed));
    
    return true;
}

// ============================================================================
// Phase 8: 两阶段运行
// ============================================================================

bool SimulationController::run() {
    // 阶段 I: 预热
    if (config_.run_spin_up) {
        if (!runSpinUp()) {
            std::cerr << "Spin-up phase failed!\n";
            return false;
        }
    }
    
    // 阶段 II: 瞬态模拟
    if (!runTransient()) {
        std::cerr << "Transient phase failed!\n";
        return false;
    }
    
    return true;
}

bool SimulationController::runSpinUp() {
    if (plots_.empty()) {
        std::cerr << "Error: No plots to simulate.\n";
        return false;
    }
    
    current_phase_ = SimulationPhase::SPIN_UP;
    
    if (config_.verbose) {
        std::cout << "\n=== Starting Spin-up Phase ===\n";
        std::cout << "Years: 1 - " << config_.spin_up_years << "\n";
        std::cout << "Saturated seed rain: Years 1 - " << config_.spin_up_saturated_years << "\n";
        std::cout << "Climate resampling window: " << config_.climate_resample_window << " years\n";
        std::cout << "Disturbance: DISABLED\n\n";
    }
    
    for (int step = 1; step <= config_.spin_up_years; ++step) {
        // 预热阶段不推进真实年份，使用随机重采样
        // 这里使用 step 作为虚拟年份
        auto year_stats = runYear(step, SimulationPhase::SPIN_UP, step);
        
        // 仅在配置允许时输出预热数据
        if (config_.output_spin_up) {
            for (const auto& stats : year_stats) {
                all_stats_.push_back(stats);
                if (output_callback_) {
                    output_callback_(stats);
                }
            }
        }
        
        // 进度输出
        if (config_.verbose && (step % 50 == 0 || step == config_.spin_up_years)) {
            if (!year_stats.empty()) {
                std::cout << "[Spin-up] Step " << std::setw(4) << step 
                          << " | Trees: " << year_stats[0].tree_count
                          << " | Saplings: " << year_stats[0].sapling_count
                          << " | LAI: " << std::fixed << std::setprecision(1) 
                          << year_stats[0].total_lai << "\n";
            }
        }
    }
    
    if (config_.verbose) {
        std::cout << "\n=== Spin-up Phase Complete ===\n";
        std::cout << "Vegetation state preserved for transient simulation.\n\n";
    }
    
    return true;
}

bool SimulationController::runTransient() {
    if (plots_.empty()) {
        std::cerr << "Error: No plots to simulate.\n";
        return false;
    }
    
    if (!species_mgr_) {
        std::cerr << "Error: Species manager not set.\n";
        return false;
    }
    
    current_phase_ = SimulationPhase::TRANSIENT;
    
    if (config_.verbose) {
        std::cout << "\n=== Starting Transient Simulation ===\n";
        std::cout << "Years: " << config_.start_year << " - " << config_.end_year << "\n";
        std::cout << "Plots: " << plots_.size() << "\n";
        std::cout << "Fire enabled: " << (config_.enable_fire ? "Yes" : "No") << "\n\n";
    }
    
    for (int year = config_.start_year; year <= config_.end_year; ++year) {
        auto year_stats = runYear(year, SimulationPhase::TRANSIENT, 0);
        
        for (const auto& stats : year_stats) {
            all_stats_.push_back(stats);
            if (output_callback_) {
                output_callback_(stats);
            }
        }
        
        // 控制台输出
        if (config_.verbose && !year_stats.empty()) {
            std::cout << OutputWriter::formatConsoleOutput(year_stats[0]) << "\n";
        }
    }
    
    if (config_.verbose) {
        std::cout << "\n=== Transient Simulation Complete ===\n";
    }
    
    return true;
}

// ============================================================================
// 年度处理
// ============================================================================

std::vector<YearlyStats> SimulationController::runYear(int year, SimulationPhase phase, int spin_up_step) {
    current_year_ = year;
    std::vector<YearlyStats> results;
    
    for (auto& plot : plots_) {
        auto stats = processPlotYear(plot, year, phase, spin_up_step);
        results.push_back(stats);
    }
    
    return results;
}

YearlyStats SimulationController::processPlotYear(Plot& plot, int year, 
                                                   SimulationPhase phase, int spin_up_step) {
    // 获取气候 ID
    int climate_id = plot.getSiteInfo().climate_id;
    
    // Step 1: 处理气候
    processClimate(year, climate_id, phase, spin_up_step);
    
    // Step 2: 干扰检查（仅瞬态阶段）
    bool fire_occurred = false;
    double fire_biomass = 0.0;
    if (phase == SimulationPhase::TRANSIENT) {
        fire_occurred = processDisturbance(plot, year, phase);
        if (fire_occurred) {
            fire_biomass = plot.getLastFireTreesKilled() * 100.0; // 估算
        }
    }
    
    // Step 3: 光照计算
    processLight(plot);
    
    // Step 4: 确定种子散布模式并执行
    SeedDispersalMode seed_mode = determineSeedMode(plot, phase, spin_up_step);
    processSeedDispersal(plot, seed_mode);
    
    // Step 5: 定居
    int new_seedlings = processEstablishment(plot);
    
    // Step 6: 生长
    std::vector<GrowthResponse> sapling_resp, tree_resp;
    processGrowth(plot, sapling_resp, tree_resp);
    
    // Step 7: 死亡
    auto [dead_saplings, dead_trees] = processMortality(plot, sapling_resp, tree_resp);
    
    // Step 8: 晋升
    int recruited = processRecruitment(plot);
    
    // 收集统计
    YearlyStats stats = collectStats(plot, year);
    stats.new_seedlings = new_seedlings;
    stats.recruited_trees = recruited;
    stats.dead_trees = dead_trees;
    stats.dead_saplings = dead_saplings;
    stats.fire_occurred = fire_occurred;
    stats.fire_biomass_removed = fire_biomass;
    
    return stats;
}

// ============================================================================
// 种子散布模式确定
// ============================================================================

SeedDispersalMode SimulationController::determineSeedMode(const Plot& plot, 
                                                           SimulationPhase phase, 
                                                           int spin_up_step) {
    if (phase == SimulationPhase::SPIN_UP) {
        // 预热阶段：前 N 年饱和，之后内部种源
        if (spin_up_step <= config_.spin_up_saturated_years) {
            return SeedDispersalMode::MODE_C_SATURATED;
        } else {
            return SeedDispersalMode::MODE_A_INTERNAL;
        }
    } else {
        // 瞬态阶段：检查火后恢复状态
        if (plot.isInFireRecovery()) {
            // 火后恢复期：内部 + 外部种源
            return SeedDispersalMode::MODE_B_EXTERNAL;
        } else {
            return SeedDispersalMode::MODE_A_INTERNAL;
        }
    }
}

Year SimulationController::getResampledClimateYear(int climate_id) {
    if (!climate_mgr_) return config_.start_year;
    
    auto years = climate_mgr_->getYearList(climate_id);
    if (years.empty()) return config_.start_year;
    
    // 限制在重采样窗口内
    int window_size = std::min(config_.climate_resample_window, static_cast<int>(years.size()));
    std::uniform_int_distribution<int> dist(0, window_size - 1);
    return years[dist(rng_)];
}

// ============================================================================
// 干扰处理
// ============================================================================

bool SimulationController::processDisturbance(Plot& plot, int year, SimulationPhase phase) {
    // 预热阶段强制关闭干扰
    if (phase == SimulationPhase::SPIN_UP) {
        return false;
    }
    
    if (!config_.enable_fire) {
        return false;
    }
    
    // 检查 EventManager
    if (event_mgr_ && event_mgr_->hasFireEvent(year, plot.getPlotId())) {
        FireDisturbance::applyFire(plot, year);
        return true;
    }
    
    // 随机火灾（如果配置了概率）
    if (config_.fire_probability > 0.0 && 
        FireDisturbance::checkFireOccurrence(config_.fire_probability)) {
        FireDisturbance::applyFire(plot, year);
        return true;
    }
    
    return false;
}

// ============================================================================
// 气候处理
// ============================================================================

void SimulationController::processClimate(int year, int climate_id, 
                                           SimulationPhase phase, 
                                           [[maybe_unused]] int spin_up_step) {
    current_climate_id_ = climate_id;
    
    if (!climate_mgr_) {
        // 使用默认值
        current_gdd_ = 1200.0;
        current_ngs_temp_ = -5.0;
        current_drought_index_ = 0.3;
        current_annual_precip_ = 800.0;
        return;
    }
    
    // 确定要使用的气候年份
    Year climate_year;
    if (phase == SimulationPhase::SPIN_UP) {
        // 随机重采样
        climate_year = getResampledClimateYear(climate_id);
    } else {
        // 时序驱动
        climate_year = year;
    }
    
    // 计算气候指标
    current_gdd_ = 0.0;
    current_annual_precip_ = 0.0;
    double ngs_temp_sum = 0.0;
    int ngs_months = 0;
    
    for (int month = 1; month <= 12; ++month) {
        const auto& mc = climate_mgr_->getMonthlyClimate(climate_id, climate_year, month);
        
        current_annual_precip_ += mc.prec;
        
        if (mc.temp > 5.5) {
            current_gdd_ += (mc.temp - 5.5) * 30.0;
        } else {
            ngs_temp_sum += mc.temp;
            ngs_months++;
        }
    }
    
    current_ngs_temp_ = ngs_months > 0 ? ngs_temp_sum / ngs_months : 0.0;
}

// ============================================================================
// 光照处理
// ============================================================================

void SimulationController::processLight(Plot& plot) {
    auto& trees = plot.getTrees();
    
    if (trees.empty()) {
        // 无树木时，地面光照为 1.0
        for (int y = 0; y < Plot::SIZE_M; ++y) {
            for (int x = 0; x < Plot::SIZE_M; ++x) {
                plot.getCells().get(x, y).light_at_ground = 1.0;
            }
        }
        return;
    }
    
    // 更新空间哈希
    light_engine_.rebuildSpatialIndex(trees);
    
    // 计算全部光照（地表和树木）
    light_engine_.calculateAllLight(trees, nullptr);
    
    // 复制地表光照到样方网格
    const auto& ground_light = light_engine_.getGroundLight();
    for (int y = 0; y < Plot::SIZE_M; ++y) {
        for (int x = 0; x < Plot::SIZE_M; ++x) {
            plot.getCells().get(x, y).light_at_ground = ground_light.get(x, y);
        }
    }
}

// ============================================================================
// 种子散布
// ============================================================================

void SimulationController::processSeedDispersal(Plot& plot, SeedDispersalMode mode) {
    if (!species_mgr_) return;
    
    const auto& trees = plot.getTrees();
    int num_species = static_cast<int>(species_mgr_->getSpeciesCount());
    
    // 初始化种子网格
    for (int y = 0; y < Plot::SIZE_M; ++y) {
        for (int x = 0; x < Plot::SIZE_M; ++x) {
            CellState& cell = plot.getCells().get(x, y);
            cell.seed_by_species.assign(num_species, 0.0);
        }
    }
    
    // 根据模式处理
    if (mode == SeedDispersalMode::MODE_C_SATURATED) {
        // 模式 C：饱和种子雨，强制 P_seed = 1.0
        for (int sp = 1; sp < num_species; ++sp) {
            for (int y = 0; y < Plot::SIZE_M; ++y) {
                for (int x = 0; x < Plot::SIZE_M; ++x) {
                    // 饱和阈值（100 seeds/m²）
                    plot.getCells().get(x, y).seed_by_species[sp] = 100.0;
                }
            }
        }
        return;
    }
    
    // 模式 A/B：计算实际种子散布
    for (int sp = 1; sp < num_species; ++sp) {
        if (!species_mgr_->hasSpecies(sp)) continue;
        const auto& species = species_mgr_->getSpecies(sp);
        
        // 样方内散布
        for (const auto& tree : trees) {
            if (!tree.isAlive() || tree.getSpeciesId() != sp) continue;
            if (tree.getAge() < species.maturity_age_yr) continue;
            
            double seed_output = tree.getLeafArea() * species.fecundity_f;
            
            for (int gy = 0; gy < Plot::SIZE_M; ++gy) {
                for (int gx = 0; gx < Plot::SIZE_M; ++gx) {
                    double dist = std::sqrt(
                        std::pow((gx + 0.5) - tree.getX(), 2) +
                        std::pow((gy + 0.5) - tree.getY(), 2)
                    );
                    
                    if (dist < 0.1) dist = 0.1;
                    
                    double kernel_val = SeedKernel::evaluate(dist, species);
                    double seeds = seed_output * kernel_val;
                    
                    if (sp < static_cast<int>(plot.getCells().get(gx, gy).seed_by_species.size())) {
                        plot.getCells().get(gx, gy).seed_by_species[sp] += seeds;
                    }
                }
            }
        }
        
        // 样方间散布（如果有距离矩阵）
        if (dist_matrix_ && mode == SeedDispersalMode::MODE_A_INTERNAL) {
            for (const auto& other_plot : plots_) {
                if (other_plot.getPlotId() == plot.getPlotId()) continue;
                
                double distance = dist_matrix_->getDistance(
                    other_plot.getPlotId(), plot.getPlotId()
                );
                
                // 检查隔离（负距离表示物理隔离）
                if (distance < 0) continue;
                
                // 计算来自其他样方的种子输入
                double total_seed_from_other = 0.0;
                for (const auto& tree : other_plot.getTrees()) {
                    if (!tree.isAlive() || tree.getSpeciesId() != sp) continue;
                    if (tree.getAge() < species.maturity_age_yr) continue;
                    total_seed_from_other += tree.getLeafArea() * species.fecundity_f;
                }
                
                double kernel_val = SeedKernel::evaluate(distance, species);
                double immigration = total_seed_from_other * kernel_val / Plot::AREA_M2;
                
                // 均匀分配到所有网格
                for (int y = 0; y < Plot::SIZE_M; ++y) {
                    for (int x = 0; x < Plot::SIZE_M; ++x) {
                        if (sp < static_cast<int>(plot.getCells().get(x, y).seed_by_species.size())) {
                            plot.getCells().get(x, y).seed_by_species[sp] += immigration;
                        }
                    }
                }
            }
        }
        
        // 模式 B：外部种源（火后恢复）
        if (mode == SeedDispersalMode::MODE_B_EXTERNAL) {
            // 假设外部种源提供背景种子雨
            double external_input = species.fecundity_f * 0.01;  // 背景值
            for (int y = 0; y < Plot::SIZE_M; ++y) {
                for (int x = 0; x < Plot::SIZE_M; ++x) {
                    if (sp < static_cast<int>(plot.getCells().get(x, y).seed_by_species.size())) {
                        plot.getCells().get(x, y).seed_by_species[sp] += external_input;
                    }
                }
            }
        }
    }
}

// ============================================================================
// 定居
// ============================================================================

int SimulationController::processEstablishment(Plot& plot) {
    if (!species_mgr_) return 0;
    
    auto& saplings = plot.getSaplings();
    int new_seedlings = 0;
    int num_species = static_cast<int>(species_mgr_->getSpeciesCount());
    
    // 记录每个网格已有的幼苗
    Grid<std::set<SpeciesId>> occupied(Plot::SIZE_M, Plot::SIZE_M);
    for (const auto& sapling : saplings) {
        if (sapling.isAlive()) {
            int cx = static_cast<int>(sapling.getX());
            int cy = static_cast<int>(sapling.getY());
            if (cx >= 0 && cx < Plot::SIZE_M && cy >= 0 && cy < Plot::SIZE_M) {
                occupied.get(cx, cy).insert(sapling.getSpeciesId());
            }
        }
    }
    
    auto& rng = RandomGenerator::getInstance();
    
    for (int sp = 1; sp < num_species; ++sp) {
        if (!species_mgr_->hasSpecies(sp)) continue;
        const auto& species = species_mgr_->getSpecies(sp);
        
        for (int y = 0; y < Plot::SIZE_M; ++y) {
            for (int x = 0; x < Plot::SIZE_M; ++x) {
                // 检查是否已有该物种幼苗
                if (occupied.get(x, y).count(sp) > 0) continue;
                
                CellState& cell = plot.getCells().get(x, y);
                
                // 种子限制
                double seed_density = 0.0;
                if (sp < static_cast<int>(cell.seed_by_species.size())) {
                    seed_density = cell.seed_by_species[sp];
                }
                double p_seed = std::min(seed_density / 100.0, 1.0);
                if (p_seed <= 0.0) continue;
                
                // 环境过滤
                bool f_temp = true;
                if (climate_mgr_) {
                    f_temp = climate_mgr_->checkNGSTemperature(
                        current_year_, species.ngs_tmin_c, species.ngs_tmax_c
                    );
                }
                
                bool f_gdd = current_gdd_ >= species.dd_min;
                bool f_drought = current_drought_index_ <= species.dr_tol;
                bool f_light = cell.light_at_ground >= species.l_min &&
                               cell.light_at_ground <= species.l_max;
                
                if (!f_temp || !f_gdd || !f_drought || !f_light) continue;
                
                // 定居概率
                if (rng.getUniform01() < p_seed) {
                    double px = x + rng.getUniform01();
                    double py = y + rng.getUniform01();
                    
                    Sapling new_sapling(next_sapling_id_++, sp, px, py, 
                                        GlobalConfig::SEEDLING_INIT_HEIGHT_M);
                    saplings.push_back(std::move(new_sapling));
                    
                    occupied.get(x, y).insert(sp);
                    new_seedlings++;
                }
            }
        }
    }
    
    return new_seedlings;
}

// ============================================================================
// 生长
// ============================================================================

void SimulationController::processGrowth(Plot& plot,
                                          std::vector<GrowthResponse>& sapling_responses,
                                          std::vector<GrowthResponse>& tree_responses) {
    if (!species_mgr_) return;
    
    sapling_responses.clear();
    tree_responses.clear();
    
    // 幼苗生长（Phase 9: 使用新的 calcSaplingLight 接口）
    for (auto& sapling : plot.getSaplings()) {
        if (!sapling.isAlive()) continue;
        
        int sp_id = sapling.getSpeciesId();
        if (!species_mgr_->hasSpecies(sp_id)) continue;
        const auto& species = species_mgr_->getSpecies(sp_id);
        
        // Phase 9: 使用 LightEngine::calcSaplingLight 计算幼树光照
        // 分层逻辑：
        // - H >= 1.0m: 仅受成树遮阴
        // - H < 1.0m: 受成树+草本双重遮阴
        int cx = static_cast<int>(sapling.getX());
        int cy = static_cast<int>(sapling.getY());
        cx = std::max(0, std::min(Plot::SIZE_M - 1, cx));
        cy = std::max(0, std::min(Plot::SIZE_M - 1, cy));
        
        double herb_lai = plot.getHerbLayer().getLAI(cx, cy);
        double light = light_engine_.calcSaplingLight(sapling, herb_lai);
        
        // 计算环境因子
        double f_light = EnvironmentResponse::calcLightResponse(light, species.shade_tol);
        double f_temp = EnvironmentResponse::calcTemperatureResponse(current_gdd_, species.dd_min);
        double f_drought = EnvironmentResponse::calcDroughtResponse(current_drought_index_, species.dr_tol);
        double f_env = f_light * f_temp * f_drought;
        
        // Bertalanffy 生长
        double delta_h = SaplingGrowthCalculator::calcHeightIncrement(
            sapling.getHeight(), species.h_max_m, species.gs, f_env
        );
        
        sapling.grow(delta_h);
        sapling.incrementAge();
        
        GrowthResponse resp;
        resp.entity_id = sapling.getId();
        resp.f_env = f_env;
        resp.growth_rate = delta_h;
        sapling_responses.push_back(resp);
    }
    
    // 成树生长
    for (auto& tree : plot.getTrees()) {
        if (!tree.isAlive()) continue;
        
        int sp_id = tree.getSpeciesId();
        if (!species_mgr_->hasSpecies(sp_id)) continue;
        const auto& species = species_mgr_->getSpecies(sp_id);
        
        double light = tree.getAvailableLight();
        
        double f_light = EnvironmentResponse::calcLightResponse(light, species.shade_tol);
        double f_temp = EnvironmentResponse::calcTemperatureResponse(current_gdd_, species.dd_min);
        double f_drought = EnvironmentResponse::calcDroughtResponse(current_drought_index_, species.dr_tol);
        double f_env = f_light * f_temp * f_drought;
        
        double delta_dbh = TreeGrowthCalculator::calcOptimalDbhIncrement(
            tree.getDbh(), tree.getHeight(), 
            species.dbh_max_cm, species.h_max_m, 
            species.ga, species.opt_s
        ) * f_env;
        
        double new_dbh = tree.getDbh() + delta_dbh;
        double new_height = TreeGrowthCalculator::calcHeightFromDbh(
            new_dbh, species.h_max_m, species.opt_s
        );
        
        tree.setDbh(new_dbh);
        tree.setHeight(new_height);
        tree.incrementAge();
        
        // 更新几何
        tree.setLeafArea(TreeGrowthCalculator::calcLeafArea(
            new_dbh, species.kc2, species.ka1, species.ka2
        ));
        tree.setCrownRadius(TreeGrowthCalculator::calcCrownRadius(
            new_dbh, species.r_a, species.r_b
        ));
        
        // 更新生物量
        double biomass = tree.calcBiomass(species.b0, species.b1, species.b2);
        tree.setBiomass(biomass);
        
        GrowthResponse resp;
        resp.entity_id = tree.getId();
        resp.f_env = f_env;
        resp.growth_rate = delta_dbh;
        tree_responses.push_back(resp);
    }
    
    // 草本生长
    double herb_light = 0.0;
    int count = 0;
    for (int y = 0; y < Plot::SIZE_M; ++y) {
        for (int x = 0; x < Plot::SIZE_M; ++x) {
            herb_light += plot.getCells().get(x, y).light_at_ground;
            count++;
        }
    }
    herb_light /= count;
    
    plot.getHerbLayer().grow(herb_light, current_gdd_, current_drought_index_);
}

// ============================================================================
// 死亡
// ============================================================================

std::pair<int, int> SimulationController::processMortality(
    Plot& plot,
    const std::vector<GrowthResponse>& sapling_responses,
    const std::vector<GrowthResponse>& tree_responses) {
    
    int dead_saplings = mortality_engine_.processSaplingMortality(
        plot.getSaplings(), sapling_responses
    );
    
    int dead_trees = mortality_engine_.processTreeMortality(
        plot.getTrees(), tree_responses
    );
    
    // 清理死亡个体
    auto& saplings = plot.getSaplings();
    saplings.erase(
        std::remove_if(saplings.begin(), saplings.end(),
                       [](const Sapling& s) { return !s.isAlive(); }),
        saplings.end()
    );
    
    auto& trees = plot.getTrees();
    trees.erase(
        std::remove_if(trees.begin(), trees.end(),
                       [](const AdultTree& t) { return !t.isAlive(); }),
        trees.end()
    );
    
    return {dead_saplings, dead_trees};
}

// ============================================================================
// 晋升
// ============================================================================

int SimulationController::processRecruitment(Plot& plot) {
    if (!species_mgr_) return 0;
    
    auto& saplings = plot.getSaplings();
    auto& trees = plot.getTrees();
    int recruited = 0;
    
    for (auto it = saplings.begin(); it != saplings.end(); ) {
        if (it->isAlive() && it->getHeight() >= GlobalConfig::HEIGHT_THRESHOLD_M) {
            int sp_id = it->getSpeciesId();
            if (species_mgr_->hasSpecies(sp_id)) {
                const auto& species = species_mgr_->getSpecies(sp_id);
                
                // Phase 8: 初始化 DBH = 0.1 cm（禁止 DBH=0）
                double init_dbh = 0.1;
                
                AdultTree new_tree(
                    next_tree_id_++,
                    sp_id,
                    it->getX(),
                    it->getY(),
                    init_dbh,
                    it->getHeight(),
                    it->getAge()
                );
                
                new_tree.setLeafArea(TreeGrowthCalculator::calcLeafArea(
                    init_dbh, species.kc2, species.ka1, species.ka2
                ));
                new_tree.setCrownRadius(TreeGrowthCalculator::calcCrownRadius(
                    init_dbh, species.r_a, species.r_b
                ));
                
                trees.push_back(std::move(new_tree));
                recruited++;
            }
            
            it = saplings.erase(it);
        } else {
            ++it;
        }
    }
    
    return recruited;
}

// ============================================================================
// 统计收集
// ============================================================================

YearlyStats SimulationController::collectStats(const Plot& plot, int year) {
    YearlyStats stats;
    stats.year = year;
    stats.plot_id = plot.getPlotId();
    stats.plot_name = plot.getSiteInfo().plot_name;
    stats.region_group = plot.getSiteInfo().region_group;
    
    stats.annual_gdd = current_gdd_;
    stats.annual_precip = current_annual_precip_;
    stats.ngs_temperature = current_ngs_temp_;
    stats.annual_drought_index = current_drought_index_;
    
    for (const auto& tree : plot.getTrees()) {
        if (tree.isAlive()) {
            stats.tree_count++;
            stats.tree_biomass_kg += tree.getBiomass();
            stats.tree_lai += tree.getLeafArea() / Plot::AREA_M2;
        }
    }
    
    for (const auto& sapling : plot.getSaplings()) {
        if (sapling.isAlive()) {
            stats.sapling_count++;
            stats.sapling_biomass_kg += sapling.getBiomass();
        }
    }
    
    stats.herb_biomass_kg = plot.getHerbLayer().getTotalBiomass();
    stats.herb_lai = plot.getHerbLayer().getMeanLAI();
    
    stats.total_biomass_kg = stats.tree_biomass_kg + stats.sapling_biomass_kg + stats.herb_biomass_kg;
    stats.total_lai = stats.tree_lai + stats.herb_lai;
    
    return stats;
}

// ============================================================================
// OutputWriter
// ============================================================================

void OutputWriter::writeSummaryHeader() {
    std::ofstream file(output_dir_ + "simulation_summary.csv");
    if (file.is_open()) {
        file << "year,plot_id,plot_name,region_group,tree_count,sapling_count,"
             << "tree_biomass_kg,sapling_biomass_kg,herb_biomass_kg,total_biomass_kg,"
             << "tree_lai,herb_lai,total_lai,new_seedlings,recruited_trees,"
             << "dead_trees,dead_saplings,fire_occurred,annual_gdd,annual_precip\n";
    }
}

void OutputWriter::appendSummary(const YearlyStats& stats) {
    std::ofstream file(output_dir_ + "simulation_summary.csv", std::ios::app);
    if (file.is_open()) {
        file << stats.year << ","
             << stats.plot_id << ","
             << stats.plot_name << ","
             << stats.region_group << ","
             << stats.tree_count << ","
             << stats.sapling_count << ","
             << std::fixed << std::setprecision(2) << stats.tree_biomass_kg << ","
             << stats.sapling_biomass_kg << ","
             << stats.herb_biomass_kg << ","
             << stats.total_biomass_kg << ","
             << std::setprecision(3) << stats.tree_lai << ","
             << stats.herb_lai << ","
             << stats.total_lai << ","
             << stats.new_seedlings << ","
             << stats.recruited_trees << ","
             << stats.dead_trees << ","
             << stats.dead_saplings << ","
             << (stats.fire_occurred ? 1 : 0) << ","
             << std::setprecision(1) << stats.annual_gdd << ","
             << stats.annual_precip << "\n";
    }
}

} // namespace sfe
