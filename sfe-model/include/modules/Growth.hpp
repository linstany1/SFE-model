/**
 * @file Growth.hpp
 * @brief 生长模块
 * 
 * 亚高山森林演替模型 (SFE-Model)
 * 
 * Phase 6: 过程模块 - 生长
 * 
 * 包含：
 * - Task 6.1: 草本生长 (Logistic方程)
 * - Task 6.2: 幼树生长 (Bertalanffy方程)
 * - Task 6.3: 成树生长 (Botkin/Leemans异速生长)
 * 
 * 关键公式：
 * - 光响应 f_light (公式35-37)
 * - 温度响应 f_temperature (公式78-79)
 * - 干旱响应 f_drought (公式73)
 * - Logistic生长 (公式101)
 * - Bertalanffy生长 (公式86, 106)
 * - 最优胸径增量 (公式105)
 */

#ifndef SFE_GROWTH_HPP
#define SFE_GROWTH_HPP

#include <vector>
#include <cmath>
#include <algorithm>
#include "core/Grid.hpp"
#include "core/GlobalConfig.hpp"
#include "core/Types.hpp"
#include "io/SpeciesProfile.hpp"
#include "modules/Mortality.hpp"

namespace sfe {

// 前向声明
class AdultTree;
class Sapling;
class HerbLayer;
class Plot;

// ============================================================================
// 环境响应函数
// ============================================================================

/**
 * @class EnvironmentResponse
 * @brief 环境响应函数计算器
 * 
 * 计算光照、温度、干旱对生长的限制
 */
class EnvironmentResponse {
public:
    // ========================================================================
    // 光照响应 (公式35-37)
    // ========================================================================
    
    /**
     * @brief 计算耐阴种光响应 (ST=5)
     * @param available_light 可用光照 [0, 1]
     * @return 光响应值 [0, 1]
     * 
     * 公式(35): f_tol = 1.0106 × (1 - e^(-4.6 × (AL - 0.03)))
     */
    static double calcLightResponseTolerant(double available_light) {
        return 1.0106 * (1.0 - std::exp(-4.6 * (available_light - 0.03)));
    }
    
    /**
     * @brief 计算喜光种光响应 (ST=1)
     * @param available_light 可用光照 [0, 1]
     * @return 光响应值 [0, 1]
     * 
     * 公式(36): f_intol = 1.60 × (1 - e^(-1.16 × (AL - 0.15)))
     */
    static double calcLightResponseIntolerant(double available_light) {
        return 1.60 * (1.0 - std::exp(-1.16 * (available_light - 0.15)));
    }
    
    /**
     * @brief 计算光响应（按耐阴等级插值）
     * @param available_light 可用光照 [0, 1]
     * @param shade_tolerance 耐阴等级 [1, 5]
     * @return 光响应值 [0, 1]
     * 
     * 公式(37): f_light = max(0, min(1, f_intol + (ST-1)/4 × (f_tol - f_intol)))
     */
    static double calcLightResponse(double available_light, double shade_tolerance) {
        double f_tol = calcLightResponseTolerant(available_light);
        double f_intol = calcLightResponseIntolerant(available_light);
        
        double f_light = f_intol + (shade_tolerance - 1.0) / 4.0 * (f_tol - f_intol);
        return std::max(0.0, std::min(1.0, f_light));
    }
    
    // ========================================================================
    // 温度响应 (公式78-79)
    // ========================================================================
    
    /**
     * @brief 计算树木温度响应
     * @param gdd 生长度日
     * @param dd_min 最小生长度日需求
     * @return 温度响应值 [0, 1]
     * 
     * 公式(78): f_temperature = max(1 - e^((DD_min - DD) × 0.0013), 0)
     */
    static double calcTemperatureResponse(double gdd, double dd_min) {
        return std::max(1.0 - std::exp((dd_min - gdd) * 0.0013), 0.0);
    }
    
    /**
     * @brief 计算草本温度响应
     * @param gdd 生长度日
     * @param dd_min 最小生长度日需求
     * @return 温度响应值 [0, 1]
     * 
     * 公式(79/96): f_temperature_herb = max(1 - e^((DD_min - DD) × 0.0027), 0)
     */
    static double calcTemperatureResponseHerb(double gdd, double dd_min) {
        return std::max(1.0 - std::exp((dd_min - gdd) * 0.0027), 0.0);
    }
    
    // ========================================================================
    // 干旱响应 (公式73, 97)
    // ========================================================================
    
    /**
     * @brief 计算干旱响应
     * @param drought_index 干旱指数 DI [0, 1]
     * @param drought_tolerance 耐旱阈值 DrTol [0, 1]
     * @return 干旱响应值 [0, 1]
     * 
     * 公式(73/97): f_drought = sqrt(max(1 - DI/DrTol, 0))
     */
    static double calcDroughtResponse(double drought_index, double drought_tolerance) {
        if (drought_tolerance <= 0.0) return 0.0;
        double ratio = drought_index / drought_tolerance;
        return std::sqrt(std::max(1.0 - ratio, 0.0));
    }
    
    // ========================================================================
    // 综合环境响应
    // ========================================================================
    
    /**
     * @brief 计算综合环境响应
     * @param f_light 光响应
     * @param f_temperature 温度响应
     * @param f_drought 干旱响应
     * @return 综合响应 f_env = f_light × f_temperature × f_drought
     */
    static double calcCombinedResponse(double f_light, double f_temperature, double f_drought) {
        return f_light * f_temperature * f_drought;
    }
    
    /**
     * @brief 一次性计算所有环境响应
     */
    static GrowthResponse calcAllResponses(double available_light, double shade_tolerance,
                                            double gdd, double dd_min,
                                            double drought_index, double drought_tolerance) {
        GrowthResponse response;
        response.f_light = calcLightResponse(available_light, shade_tolerance);
        response.f_temperature = calcTemperatureResponse(gdd, dd_min);
        response.f_drought = calcDroughtResponse(drought_index, drought_tolerance);
        return response;
    }
};

// ============================================================================
// Task 6.1: 草本生长
// ============================================================================

/**
 * @struct HerbGrowthParams
 * @brief 草本生长参数
 */
struct HerbGrowthParams {
    double k_max = 0.8;        ///< 最大生物量 (kg/m²)
    double r_max = 1.0;        ///< 最大生长速率
    double dd_min = 200.0;     ///< 最小生长度日
    double dr_tol = 0.5;       ///< 耐旱阈值
    double shade_tol = 3.0;    ///< 耐阴等级
    double min_biomass = 0.02; ///< 最小生物量 (kg/m²)
};

/**
 * @class HerbGrowthCalculator
 * @brief 草本生长计算器
 * 
 * 使用 Logistic 生长方程
 */
class HerbGrowthCalculator {
public:
    static constexpr double CELL_AREA = 1.0;  ///< 单元面积 (m²)
    static constexpr double MIN_BIOMASS = 0.02;  ///< 最小生物量 (kg/m²)
    
    /**
     * @brief 计算环境限制后的容纳量 K
     * @param k_max 最大容纳量 (kg/m²)
     * @param f_env 环境响应因子
     * @return 有效容纳量 (kg)
     * 
     * 公式(99): K = (K_max × A_cell / 1000) × f_env
     * 注意：K_max 单位为 kg/m²，结果为 kg
     */
    static double calcEffectiveK(double k_max, double f_env) {
        // K_max 已经是 kg/m²，A_cell = 1 m²
        return k_max * CELL_AREA * f_env;
    }
    
    /**
     * @brief 计算有效生长速率
     * @param r_max 最大生长速率
     * @param f_env 环境响应因子
     * @return 有效生长速率
     * 
     * 公式(100): r = r_max × f_env
     */
    static double calcEffectiveR(double r_max, double f_env) {
        return r_max * f_env;
    }
    
    /**
     * @brief 计算单网格草本生长
     * @param current_biomass 当前生物量 (kg/m²)
     * @param available_light 可用光照 [0, 1]
     * @param drought_index 干旱指数 [0, 1]
     * @param gdd 生长度日
     * @param params 生长参数
     * @return 新的生物量 (kg/m²)
     * 
     * 公式(101): B_{Y+1} = B_Y + r × (1 - B/K) × B
     * 公式(102): B_{Y+1} = max(B_{Y+1}, 0.02)
     */
    static double calcGrowth(double current_biomass, 
                              double available_light,
                              double drought_index,
                              double gdd,
                              const HerbGrowthParams& params);
    
    /**
     * @brief 使用物种参数计算生长
     */
    static double calcGrowth(double current_biomass,
                              double available_light,
                              double drought_index,
                              double gdd,
                              const SpeciesProfile& herb_species);
};

// ============================================================================
// Task 6.2: 幼树生长
// ============================================================================

/**
 * @class SaplingGrowthCalculator
 * @brief 幼树生长计算器
 * 
 * 使用 Bertalanffy 方程
 */
class SaplingGrowthCalculator {
public:
    static constexpr double INITIAL_HEIGHT = 0.05;  ///< 幼苗初始高度 (m)
    static constexpr double RECRUIT_HEIGHT = 1.37;  ///< 晋升高度阈值 (m)
    static constexpr double SAPLING_DROUGHT_SENSITIVITY = 1.5;  ///< 幼树干旱敏感系数 γ_sapling
    
    /**
     * @brief 计算幼树潜在高度
     * @param current_height 当前高度 (m)
     * @param h_max 最大高度 (m)
     * @param gs 生长参数
     * @return 潜在新高度 (m)
     * 
     * 公式(86/106): H_{t+1} = H_max × [1 - (1 - (H_t/H_max)^(1/3)) × e^(-gs)]^3
     */
    static double calcPotentialHeight(double current_height, double h_max, double gs) {
        if (h_max <= 0.0 || current_height < 0.0) return current_height;
        
        double ratio = current_height / h_max;
        ratio = std::max(0.001, std::min(ratio, 0.999));  // 避免边界问题
        
        double term = 1.0 - std::pow(ratio, 1.0/3.0);
        double new_ratio = 1.0 - term * std::exp(-gs);
        
        return h_max * std::pow(new_ratio, 3.0);
    }
    
    /**
     * @brief 计算高度增量
     * @param current_height 当前高度 (m)
     * @param h_max 最大高度 (m)
     * @param gs 生长参数
     * @param f_env 环境修正因子 [0, 1]
     * @return 高度增量 (m)
     */
    static double calcHeightIncrement(double current_height, double h_max, double gs, double f_env) {
        double potential = calcPotentialHeight(current_height, h_max, gs);
        double delta = potential - current_height;
        return std::max(0.0, delta * f_env);
    }
    
    /**
     * @brief 计算幼树干旱指数（放大）
     * @param drought_index 原始干旱指数
     * @return 放大后的干旱指数
     * 
     * 公式(88): DI_sapling = min(1, γ_sapling × DI)
     */
    static double calcSaplingDroughtIndex(double drought_index) {
        return std::min(1.0, SAPLING_DROUGHT_SENSITIVITY * drought_index);
    }
    
    /**
     * @brief 计算幼树实际生长
     * @param sapling 幼树对象
     * @param available_light 可用光照 [0, 1]
     * @param drought_index 干旱指数 [0, 1]
     * @param gdd 生长度日
     * @param species 物种参数
     * @return 生长响应和新高度
     */
    static std::pair<GrowthResponse, double> calcGrowth(
        const Sapling& sapling,
        double available_light,
        double drought_index,
        double gdd,
        const SpeciesProfile& species
    );
    
    /**
     * @brief 检查是否可以晋升为成树
     */
    static bool canRecruit(double height) {
        return height >= RECRUIT_HEIGHT;
    }
};

// ============================================================================
// Task 6.3: 成树生长
// ============================================================================

/**
 * @class TreeGrowthCalculator
 * @brief 成树生长计算器
 * 
 * 使用 Botkin/Leemans 异速生长方程
 */
class TreeGrowthCalculator {
public:
    /**
     * @brief 计算最优胸径增量
     * @param dbh 当前胸径 (cm)
     * @param height 当前高度 (m)
     * @param dbh_max 最大胸径 (cm)
     * @param h_max 最大高度 (m)
     * @param ga 生长参数 ga
     * @param opt_s 生长参数 opt_s
     * @return 最优胸径增量 (cm)
     * 
     * 公式(105): ΔDBH_opt = ga × DBH × (1 - DBH×H/(DBH_max×H_max)) / 
     *            (2H + opt_s × DBH × e^(-opt_s × DBH / (H_max - 1.37)))
     */
    static double calcOptimalDbhIncrement(double dbh, double height,
                                           double dbh_max, double h_max,
                                           double ga, double opt_s) {
        if (dbh <= 0.0 || height <= 0.0 || dbh_max <= 0.0 || h_max <= 1.37) {
            return 0.0;
        }
        
        // 分子
        double size_ratio = (dbh * height) / (dbh_max * h_max);
        double numerator = ga * dbh * (1.0 - size_ratio);
        
        // 分母
        double h_range = h_max - 1.37;
        double exp_term = std::exp(-opt_s * dbh / h_range);
        double denominator = 2.0 * height + opt_s * dbh * exp_term;
        
        if (denominator <= 0.0) return 0.0;
        
        return std::max(0.0, numerator / denominator);
    }
    
    /**
     * @brief 计算实际胸径增量
     * @param dbh_opt 最优胸径增量 (cm)
     * @param f_env 环境响应因子
     * @return 实际胸径增量 (cm)
     * 
     * 公式(107): ΔDBH = ΔDBH_opt × (f_drought × f_temperature × f_light)
     */
    static double calcActualDbhIncrement(double dbh_opt, double f_env) {
        return dbh_opt * f_env;
    }
    
    /**
     * @brief 根据新胸径计算新高度
     * @param new_dbh 新胸径 (cm)
     * @param h_max 最大高度 (m)
     * @param opt_s 生长参数
     * @return 新高度 (m)
     * 
     * 公式(106下半): H = 1.37 + (H_max - 1.37) × (1 - e^(-opt_s × DBH / (H_max - 1.37)))
     */
    static double calcHeightFromDbh(double new_dbh, double h_max, double opt_s) {
        if (h_max <= 1.37 || new_dbh <= 0.0) {
            return 1.37;
        }
        
        double h_range = h_max - 1.37;
        double exp_term = std::exp(-opt_s * new_dbh / h_range);
        
        return 1.37 + h_range * (1.0 - exp_term);
    }
    
    /**
     * @brief 计算叶面积
     * @param dbh 胸径 (cm)
     * @param kc2, ka1, ka2 异速生长参数
     * @return 叶面积 (m²)
     * 
     * 公式(108): LA = KC2 × KA1 × DBH^KA2
     */
    static double calcLeafArea(double dbh, double kc2, double ka1, double ka2) {
        if (dbh <= 0.0) return 0.0;
        return kc2 * ka1 * std::pow(dbh, ka2);
    }
    
    /**
     * @brief 计算冠幅半径
     * @param dbh 胸径 (cm)
     * @param r_a, r_b 冠幅参数
     * @return 冠幅半径 (m)
     */
    static double calcCrownRadius(double dbh, double r_a, double r_b) {
        if (dbh <= 0.0) return r_b;
        return r_a * std::log(dbh) + r_b;
    }
    
    /**
     * @brief 计算成树生长（完整流程）
     * @param tree 成树对象
     * @param available_light 可用光照 [0, 1]
     * @param drought_index 干旱指数 [0, 1]
     * @param gdd 生长度日
     * @param species 物种参数
     * @return 生长响应和各属性增量
     */
    static std::pair<GrowthResponse, double> calcGrowth(
        const AdultTree& tree,
        double available_light,
        double drought_index,
        double gdd,
        const SpeciesProfile& species
    );
};

// ============================================================================
// 生长引擎
// ============================================================================

/**
 * @class GrowthEngine
 * @brief 生长引擎
 * 
 * 整合草本、幼树、成树的生长计算
 */
class GrowthEngine {
public:
    GrowthEngine() = default;
    
    /**
     * @brief 设置物种管理器
     */
    void setSpeciesManager(const SpeciesManager* species_mgr) {
        species_mgr_ = species_mgr;
    }
    
    /**
     * @brief 处理草本层生长
     * @param herb_layer 草本层
     * @param ground_light 地表光照网格
     * @param drought_index 干旱指数
     * @param gdd 生长度日
     */
    void processHerbGrowth(HerbLayer& herb_layer,
                           const Grid<double>& ground_light,
                           double drought_index,
                           double gdd);
    
    /**
     * @brief 处理幼树生长
     * @param saplings 幼树列表
     * @param ground_light 地表光照网格
     * @param drought_index 干旱指数
     * @param gdd 生长度日
     * @return 每棵幼树的生长响应
     */
    std::vector<GrowthResponse> processSaplingGrowth(
        std::vector<Sapling>& saplings,
        const Grid<double>& ground_light,
        double drought_index,
        double gdd
    );
    
    /**
     * @brief 处理成树生长
     * @param trees 成树列表
     * @param drought_index 干旱指数
     * @param gdd 生长度日
     * @return 每棵树的生长响应
     */
    std::vector<GrowthResponse> processTreeGrowth(
        std::vector<AdultTree>& trees,
        double drought_index,
        double gdd
    );
    
    /**
     * @brief 处理幼树晋升
     * @param saplings 幼树列表（会被修改）
     * @param trees 成树列表（会被修改）
     * @param next_tree_id 下一个成树ID
     * @return 晋升的树木数量
     */
    int processRecruitment(std::vector<Sapling>& saplings,
                           std::vector<AdultTree>& trees,
                           TreeId& next_tree_id);

private:
    const SpeciesManager* species_mgr_ = nullptr;
};

/**
 * @struct GrowthStats
 * @brief 生长统计
 */
struct GrowthStats {
    double herb_biomass_change = 0.0;    ///< 草本生物量变化
    double mean_sapling_height_growth = 0.0;  ///< 平均幼树高度增长
    double mean_tree_dbh_growth = 0.0;   ///< 平均成树DBH增长
    int saplings_recruited = 0;          ///< 晋升的幼树数
};

} // namespace sfe

#endif // SFE_GROWTH_HPP
