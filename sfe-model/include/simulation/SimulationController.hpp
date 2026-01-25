/**
 * @file SimulationController.hpp
 * @brief 模拟控制器
 * 
 * 亚高山森林演替模型 (SFE-Model)
 * Phase 8: 两阶段模拟策略 (Spin-up + Transient)
 */

#ifndef SFE_SIMULATION_CONTROLLER_HPP
#define SFE_SIMULATION_CONTROLLER_HPP

#include <vector>
#include <string>
#include <functional>
#include <iostream>
#include <iomanip>
#include <set>
#include <cmath>
#include <random>
#include <memory>

#include "core/GlobalConfig.hpp"
#include "core/Types.hpp"
#include "core/RandomGenerator.hpp"
#include "entities/Plot.hpp"
#include "io/SpeciesProfile.hpp"
#include "io/ClimateParser.hpp"
#include "modules/LightEngine.hpp"
#include "modules/SeedDispersal.hpp"
#include "modules/Establishment.hpp"
#include "modules/Mortality.hpp"
#include "modules/Growth.hpp"
#include "simulation/EventManager.hpp"

namespace sfe {

// ============================================================================
// 枚举与配置结构
// ============================================================================

/**
 * @enum SeedDispersalMode
 * @brief Phase 8: 种子散布模式
 */
enum class SeedDispersalMode {
    MODE_A_INTERNAL,    ///< 模式 A: 内部种源传播（样方内+样方间）
    MODE_B_EXTERNAL,    ///< 模式 B: 外部种源传播（来自 Source Plots）
    MODE_C_SATURATED    ///< 模式 C: 饱和种子雨（P_seed = 1.0）
};

/**
 * @enum SimulationPhase
 * @brief Phase 8: 模拟阶段
 */
enum class SimulationPhase {
    SPIN_UP,    ///< 预热阶段
    TRANSIENT   ///< 瞬态模拟阶段
};

/**
 * @struct SimulationConfig
 * @brief Phase 8: 模拟配置参数
 */
struct SimulationConfig {
    // === 时间控制 ===
    int start_year = 2000;              ///< 瞬态模拟起始年
    int end_year = 2100;                ///< 瞬态模拟结束年
    
    // === Phase 8: 预热阶段配置 ===
    int spin_up_years = 500;            ///< 预热阶段年数
    int spin_up_saturated_years = 50;   ///< 饱和种子雨年数（预热前期）
    int climate_resample_window = 30;   ///< 气候随机重采样窗口（年）
    bool run_spin_up = true;            ///< 是否运行预热阶段
    
    // === 干扰控制 ===
    bool enable_fire = false;           ///< 是否启用火灾（瞬态阶段）
    double fire_probability = 0.0;      ///< 随机火灾概率（若不使用 events.txt）
    
    // === 输出控制 ===
    bool verbose = true;                ///< 详细输出
    std::string output_dir = "./output/";
    int output_interval = 5;            ///< 详细输出间隔（年）
    bool output_spin_up = false;        ///< 是否输出预热阶段数据
};

/**
 * @struct YearlyStats
 * @brief 年度统计数据
 */
struct YearlyStats {
    int year = 0;
    int plot_id = 0;
    std::string plot_name = "";         ///< Phase 8: 样方名称
    std::string region_group = "";      ///< Phase 8: 区域分组
    
    double annual_gdd = 0.0;
    double annual_precip = 0.0;
    double ngs_temperature = 0.0;
    double annual_drought_index = 0.0;
    
    int tree_count = 0;
    int sapling_count = 0;
    
    double tree_biomass_kg = 0.0;
    double sapling_biomass_kg = 0.0;
    double herb_biomass_kg = 0.0;
    double total_biomass_kg = 0.0;
    
    double tree_lai = 0.0;
    double herb_lai = 0.0;
    double total_lai = 0.0;
    
    int new_seedlings = 0;
    int recruited_trees = 0;
    int dead_trees = 0;
    int dead_saplings = 0;
    
    bool fire_occurred = false;
    double fire_biomass_removed = 0.0;  ///< Phase 8: 火灾移除的生物量
};

/**
 * @class FireDisturbance
 * @brief 火灾干扰管理器
 */
class FireDisturbance {
public:
    static bool checkFireOccurrence(double probability) {
        if (probability <= 0.0) return false;
        return RandomGenerator::getInstance().getUniform01() < probability;
    }
    
    static double applyFire(Plot& plot, Year year) {
        return plot.burn(year);
    }
};

// ============================================================================
// 模拟控制器
// ============================================================================

/**
 * @class SimulationController
 * @brief Phase 8: 两阶段模拟控制器
 * 
 * 实现预热-瞬态两阶段模拟策略：
 * - 阶段 I (Spin-up): 随机重采样气候，饱和/内部种子模式
 * - 阶段 II (Transient): 时序气候，事件驱动干扰
 */
class SimulationController {
public:
    SimulationController() = default;
    
    // ========================================================================
    // 初始化
    // ========================================================================
    
    void addPlot(Plot&& plot) {
        plots_.push_back(std::move(plot));
    }
    
    void setSpeciesManager(const SpeciesManager* mgr) {
        species_mgr_ = mgr;
    }
    
    void setClimateManager(const ClimateManager* mgr) {
        climate_mgr_ = mgr;
    }
    
    void setEventManager(const EventManager* mgr) {
        event_mgr_ = mgr;
    }
    
    void setDistanceMatrix(const DistanceMatrix* matrix) {
        dist_matrix_ = matrix;
    }
    
    bool initialize(const SimulationConfig& config);
    
    // ========================================================================
    // Phase 8: 两阶段运行接口
    // ========================================================================
    
    /**
     * @brief 运行完整模拟（预热 + 瞬态）
     */
    bool run();
    
    /**
     * @brief 阶段 I: 预热循环
     * 
     * 特点：
     * - 气候：随机重采样
     * - 干扰：强制关闭
     * - 种子：前 N 年饱和，之后内部种源
     */
    bool runSpinUp();
    
    /**
     * @brief 阶段 II: 瞬态模拟循环
     * 
     * 特点：
     * - 气候：时序驱动
     * - 干扰：根据 events.txt 触发
     * - 种子：内部 + 火后外部种源
     */
    bool runTransient();
    
    /**
     * @brief 运行单一年份（内部方法）
     * @param year 年份
     * @param phase 当前阶段
     * @param spin_up_step 预热步数（仅预热阶段使用）
     * @return 年度统计列表
     */
    std::vector<YearlyStats> runYear(int year, SimulationPhase phase, int spin_up_step = 0);
    
    // ========================================================================
    // 输出与回调
    // ========================================================================
    
    void setOutputCallback(std::function<void(const YearlyStats&)> callback) {
        output_callback_ = callback;
    }
    
    const std::vector<YearlyStats>& getAllStats() const { return all_stats_; }
    
    // ========================================================================
    // 访问器
    // ========================================================================
    
    std::vector<Plot>& getPlots() { return plots_; }
    const std::vector<Plot>& getPlots() const { return plots_; }
    
    SimulationPhase getCurrentPhase() const { return current_phase_; }
    int getCurrentYear() const { return current_year_; }

private:
    // 内部处理方法
    YearlyStats processPlotYear(Plot& plot, int year, SimulationPhase phase, int spin_up_step);
    bool processDisturbance(Plot& plot, int year, SimulationPhase phase);
    void processClimate(int year, int climate_id, SimulationPhase phase, int spin_up_step);
    void processLight(Plot& plot);
    void processSeedDispersal(Plot& plot, SeedDispersalMode mode);
    int processEstablishment(Plot& plot);
    void processGrowth(Plot& plot, std::vector<GrowthResponse>& sapling_resp, 
                       std::vector<GrowthResponse>& tree_resp);
    std::pair<int, int> processMortality(Plot& plot, 
                                          const std::vector<GrowthResponse>& sapling_resp,
                                          const std::vector<GrowthResponse>& tree_resp);
    int processRecruitment(Plot& plot);
    YearlyStats collectStats(const Plot& plot, int year);
    
    /**
     * @brief 确定种子散布模式
     */
    SeedDispersalMode determineSeedMode(const Plot& plot, SimulationPhase phase, int spin_up_step);
    
    /**
     * @brief 获取随机重采样的气候年份
     */
    Year getResampledClimateYear(int climate_id);
    
    // 配置
    SimulationConfig config_;
    
    // 数据管理器（外部持有）
    const SpeciesManager* species_mgr_ = nullptr;
    const ClimateManager* climate_mgr_ = nullptr;
    const EventManager* event_mgr_ = nullptr;
    const DistanceMatrix* dist_matrix_ = nullptr;
    
    // 样方列表
    std::vector<Plot> plots_;
    
    // 处理引擎
    LightEngine light_engine_;
    GrowthEngine growth_engine_;
    MortalityEngine mortality_engine_;
    
    // 当前状态
    SimulationPhase current_phase_ = SimulationPhase::SPIN_UP;
    int current_year_ = 0;
    int current_climate_id_ = 0;
    
    // 气候缓存
    double current_gdd_ = 1200.0;
    double current_ngs_temp_ = -5.0;
    double current_drought_index_ = 0.3;
    double current_annual_precip_ = 800.0;
    
    // 统计
    std::vector<YearlyStats> all_stats_;
    std::function<void(const YearlyStats&)> output_callback_;
    
    // ID 计数器
    TreeId next_tree_id_ = 1;
    TreeId next_sapling_id_ = 10000;
    
    // 随机数生成器
    std::mt19937 rng_;
};

// ============================================================================
// 输出格式化工具
// ============================================================================

/**
 * @class OutputWriter  
 * @brief 输出写入器
 */
class OutputWriter {
public:
    explicit OutputWriter(const std::string& output_dir) : output_dir_(output_dir) {}
    
    void writeSummaryHeader();
    void appendSummary(const YearlyStats& stats);
    
    static std::string formatConsoleOutput(const YearlyStats& stats) {
        std::ostringstream oss;
        oss << "Year " << std::setw(4) << stats.year
            << " | Trees: " << std::setw(4) << stats.tree_count
            << " | Saplings: " << std::setw(4) << stats.sapling_count
            << " | LAI: " << std::fixed << std::setprecision(1) << stats.total_lai
            << " | Bio: " << std::setprecision(1) << stats.total_biomass_kg << " kg";
        if (stats.fire_occurred) {
            oss << " | FIRE";
        }
        return oss.str();
    }
    
private:
    std::string output_dir_;
};

} // namespace sfe

#endif // SFE_SIMULATION_CONTROLLER_HPP
