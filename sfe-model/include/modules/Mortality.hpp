/**
 * @file Mortality.hpp
 * @brief 死亡率计算模块
 * 
 * 亚高山森林演替模型 (SFE-Model)
 * 
 * 实现成树和幼树的死亡判定：
 * 1. 内在死亡率 P_mAge (公式110)
 * 2. 生长胁迫死亡率 P_stress (公式112)
 * 3. 总死亡概率 P_tree_mort
 * 
 * 关键公式：
 * - P_mAge = 1 - 0.01^(1/maxAge)
 * - P_stress 基于连续3年 f_env < 0.1
 * - P_tree_mort = min(P_stress + P_mAge, 1)
 */

#ifndef SFE_MORTALITY_HPP
#define SFE_MORTALITY_HPP

#include <vector>
#include <cmath>
#include <algorithm>
#include "core/GlobalConfig.hpp"
#include "core/Types.hpp"
#include "core/RandomGenerator.hpp"
#include "io/SpeciesProfile.hpp"

namespace sfe {

// 前向声明
class AdultTree;
class Sapling;

/**
 * @struct GrowthResponse
 * @brief 生长响应因子
 */
struct GrowthResponse {
    TreeId entity_id = 0;         ///< 实体ID
    double f_drought = 1.0;       ///< 干旱响应
    double f_temperature = 1.0;   ///< 温度响应
    double f_light = 1.0;         ///< 光照响应
    double f_env = 1.0;           ///< 综合环境响应
    double growth_rate = 0.0;     ///< 生长速率
    
    /**
     * @brief 计算综合环境响应
     * f_env = f_drought × f_temperature × f_light
     */
    double getEnvironmentResponse() const {
        return f_drought * f_temperature * f_light;
    }
    
    void updateFEnv() {
        f_env = getEnvironmentResponse();
    }
};

/**
 * @class MortalityCalculator
 * @brief 死亡率计算器
 */
class MortalityCalculator {
public:
    // 常量
    static constexpr double STRESS_THRESHOLD = 0.1;        ///< 胁迫阈值 (f_env < 0.1)
    static constexpr int STRESS_YEARS_THRESHOLD = 3;       ///< 连续胁迫年数阈值
    static constexpr double BASE_SURVIVAL_PROB = 0.01;     ///< 基础存活概率 (1%)
    
    /**
     * @brief 计算内在死亡率（基于最大年龄）
     * @param max_age 物种最大年龄 (年)
     * @return 年内在死亡概率 [0, 1]
     * 
     * 公式(110): P_mAge = 1 - 0.01^(1/maxAge)
     */
    static double calcIntrinsicMortality(int max_age) {
        if (max_age <= 0) return 1.0;
        return 1.0 - std::pow(BASE_SURVIVAL_PROB, 1.0 / max_age);
    }
    
    /**
     * @brief 计算生长胁迫死亡率
     * @param f_env 环境响应因子
     * @param stress_years 连续胁迫年数
     * @param max_age 物种最大年龄
     * @return 胁迫死亡概率 [0, 1]
     * 
     * 公式(112):
     * P_stress = 0                           if f_env >= 0.1
     *          = P_mAge                      if f_env < 0.1 and N < 3
     *          = (1-0.01^0.1) × (1-f_env/0.1) if f_env < 0.1 and N >= 3
     */
    static double calcStressMortality(double f_env, int stress_years, int max_age) {
        if (f_env >= STRESS_THRESHOLD) {
            return 0.0;
        }
        
        if (stress_years < STRESS_YEARS_THRESHOLD) {
            return calcIntrinsicMortality(max_age);
        }
        
        // 严重胁迫死亡率
        double severe_base = 1.0 - std::pow(BASE_SURVIVAL_PROB, 0.1);
        double stress_factor = 1.0 - f_env / STRESS_THRESHOLD;
        return severe_base * stress_factor;
    }
    
    /**
     * @brief 计算成树总死亡概率
     * @param f_env 环境响应因子
     * @param stress_years 连续胁迫年数
     * @param species 物种参数
     * @return 总死亡概率 [0, 1]
     * 
     * P_tree_mort = min(P_stress + P_mAge, 1)
     */
    static double calcTreeMortality(double f_env, int stress_years, 
                                     const SpeciesProfile& species) {
        double p_age = calcIntrinsicMortality(species.max_age_yr);
        double p_stress = calcStressMortality(f_env, stress_years, species.max_age_yr);
        return std::min(p_stress + p_age, 1.0);
    }
    
    /**
     * @brief 计算幼树死亡概率
     * @param f_env 环境响应因子
     * @param stress_years 连续胁迫年数
     * @param species 物种参数
     * @return 死亡概率 [0, 1]
     * 
     * 幼树仅使用胁迫死亡率
     */
    static double calcSaplingMortality(double f_env, int stress_years,
                                        const SpeciesProfile& species) {
        return calcStressMortality(f_env, stress_years, species.max_age_yr);
    }
    
    /**
     * @brief 判定是否死亡
     * @param mortality_prob 死亡概率
     * @param rng 随机数生成器
     * @return 是否死亡
     */
    static bool shouldDie(double mortality_prob, RandomGenerator& rng) {
        return rng.getUniform01() < mortality_prob;
    }
};

/**
 * @class StressTracker
 * @brief 胁迫追踪器
 * 
 * 用于追踪个体的连续胁迫年数
 */
class StressTracker {
public:
    /**
     * @brief 更新胁迫状态
     * @param f_env 当年环境响应因子
     * @param current_stress_years 当前连续胁迫年数
     * @return 更新后的胁迫年数
     */
    static int updateStressYears(double f_env, int current_stress_years) {
        if (f_env < MortalityCalculator::STRESS_THRESHOLD) {
            return current_stress_years + 1;
        }
        return 0;  // 重置
    }
};

/**
 * @class MortalityEngine
 * @brief 死亡判定引擎
 * 
 * 处理整个样方的死亡判定
 */
class MortalityEngine {
public:
    MortalityEngine() = default;
    
    /**
     * @brief 设置物种管理器
     */
    void setSpeciesManager(const SpeciesManager* species_mgr) {
        species_mgr_ = species_mgr;
    }
    
    /**
     * @brief 设置随机数生成器
     */
    void setRandomGenerator(RandomGenerator* rng) {
        rng_ = rng;
    }
    
    /**
     * @brief 处理成树死亡判定
     * @param trees 成树列表（会被修改：标记死亡）
     * @param growth_responses 每棵树的生长响应
     * @return 死亡树木数量
     */
    int processTreeMortality(std::vector<AdultTree>& trees,
                              const std::vector<GrowthResponse>& growth_responses);
    
    /**
     * @brief 处理幼树死亡判定
     * @param saplings 幼树列表（会被修改：标记死亡）
     * @param growth_responses 每棵幼树的生长响应
     * @return 死亡幼树数量
     */
    int processSaplingMortality(std::vector<Sapling>& saplings,
                                 const std::vector<GrowthResponse>& growth_responses);
    
    /**
     * @brief 移除死亡的成树
     * @param trees 成树列表
     * @return 移除的数量
     */
    static int removeDeadTrees(std::vector<AdultTree>& trees);
    
    /**
     * @brief 移除死亡的幼树
     * @param saplings 幼树列表
     * @return 移除的数量
     */
    static int removeDeadSaplings(std::vector<Sapling>& saplings);

private:
    const SpeciesManager* species_mgr_ = nullptr;
    RandomGenerator* rng_ = nullptr;
};

/**
 * @struct MortalityStats
 * @brief 死亡统计
 */
struct MortalityStats {
    int trees_died_intrinsic = 0;    ///< 内在死亡
    int trees_died_stress = 0;       ///< 胁迫死亡
    int saplings_died = 0;           ///< 幼树死亡
    
    int totalTreeDeaths() const {
        return trees_died_intrinsic + trees_died_stress;
    }
    
    int totalDeaths() const {
        return totalTreeDeaths() + saplings_died;
    }
};

} // namespace sfe

#endif // SFE_MORTALITY_HPP
