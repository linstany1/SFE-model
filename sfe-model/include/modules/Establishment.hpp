/**
 * @file Establishment.hpp
 * @brief 幼树定居模块
 * 
 * 亚高山森林演替模型 (SFE-Model)
 * 
 * 实现环境过滤器和定居判定：
 * 1. 非生长季温度检查 f_T_NGSes
 * 2. 生长度日检查 f_GDDes  
 * 3. 干旱检查 f_Dres
 * 4. 光照检查 f_ALes
 * 
 * 最终定居概率：P_est = P_seed × f_T_NGSes × f_GDDes × f_Dres × f_ALes
 */

#ifndef SFE_ESTABLISHMENT_HPP
#define SFE_ESTABLISHMENT_HPP

#include <vector>
#include <cmath>
#include "core/Grid.hpp"
#include "core/GlobalConfig.hpp"
#include "core/Types.hpp"
#include "core/RandomGenerator.hpp"
#include "io/SpeciesProfile.hpp"
#include "entities/Sapling.hpp"

namespace sfe {

// 前向声明
class Plot;
class PlotSeedBank;

/**
 * @struct EstablishmentFilters
 * @brief 环境过滤器结果
 */
struct EstablishmentFilters {
    double f_temp_ngs = 1.0;    ///< 非生长季温度过滤 [0, 1]
    double f_gdd = 1.0;         ///< 生长度日过滤 [0, 1]
    double f_drought = 1.0;    ///< 干旱过滤 [0, 1]
    double f_light = 1.0;       ///< 光照过滤 [0, 1]
    
    /**
     * @brief 计算综合环境概率
     */
    double getEnvironmentProbability() const {
        return f_temp_ngs * f_gdd * f_drought * f_light;
    }
};

/**
 * @struct CellEnvironment
 * @brief 网格环境条件
 */
struct CellEnvironment {
    double light_at_ground = 1.0;   ///< 地表光照 [0, 1]
    double drought_index = 0.0;     ///< 干旱指数 [0, 1]
    double gdd = 0.0;               ///< 生长度日
    double ngs_temp = 0.0;          ///< 非生长季温度 (°C)
};

/**
 * @class EstablishmentFilter
 * @brief 环境过滤器计算器
 */
class EstablishmentFilter {
public:
    /**
     * @brief 检查非生长季温度
     * @param ngs_temp 非生长季平均温度 (°C)
     * @param species 物种参数
     * @return 过滤值 [0, 1]
     * 
     * f_T_NGSes = 1 if ngsTmin <= T_ngs < ngsTmax, else 0
     */
    static double checkTemperature(double ngs_temp, const SpeciesProfile& species) {
        if (ngs_temp >= species.ngs_tmin_c && ngs_temp < species.ngs_tmax_c) {
            return 1.0;
        }
        return 0.0;
    }
    
    /**
     * @brief 检查生长度日
     * @param gdd 生长度日
     * @param species 物种参数
     * @return 过滤值 [0, 1]
     * 
     * f_GDDes = 1 if GDD >= DD_min, else 0
     */
    static double checkGDD(double gdd, const SpeciesProfile& species) {
        if (gdd >= species.dd_min) {
            return 1.0;
        }
        return 0.0;
    }
    
    /**
     * @brief 检查干旱条件
     * @param drought_index 干旱指数 [0, 1]
     * @param species 物种参数
     * @return 过滤值 [0, 1]
     * 
     * f_Dres = 1 if DI <= DrTol, else 0
     */
    static double checkDrought(double drought_index, const SpeciesProfile& species) {
        if (drought_index <= species.dr_tol) {
            return 1.0;
        }
        return 0.0;
    }
    
    /**
     * @brief 检查光照条件
     * @param light 可用光照 [0, 1]
     * @param species 物种参数
     * @return 过滤值 [0, 1]
     * 
     * f_ALes = 1 if L_min <= AL <= L_max, else 0
     */
    static double checkLight(double light, const SpeciesProfile& species) {
        if (light >= species.l_min && light <= species.l_max) {
            return 1.0;
        }
        return 0.0;
    }
    
    /**
     * @brief 计算所有环境过滤器
     * @param env 网格环境
     * @param species 物种参数
     * @return 过滤器结果
     */
    static EstablishmentFilters evaluate(const CellEnvironment& env, 
                                          const SpeciesProfile& species) {
        EstablishmentFilters filters;
        filters.f_temp_ngs = checkTemperature(env.ngs_temp, species);
        filters.f_gdd = checkGDD(env.gdd, species);
        filters.f_drought = checkDrought(env.drought_index, species);
        filters.f_light = checkLight(env.light_at_ground, species);
        return filters;
    }
};

/**
 * @struct NewSeedling
 * @brief 新定居的幼苗信息
 */
struct NewSeedling {
    SpeciesId species_id = 0;
    int cell_x = 0;
    int cell_y = 0;
    double x = 0.0;     ///< 精确x坐标
    double y = 0.0;     ///< 精确y坐标
};

/**
 * @class EstablishmentEngine
 * @brief 幼苗定居引擎
 * 
 * 处理整个样方的幼苗定居过程
 */
class EstablishmentEngine {
public:
    static constexpr int PLOT_SIZE = GlobalConfig::PLOT_SIZE_M;
    
    EstablishmentEngine() = default;
    
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
     * @brief 处理样方的幼苗定居
     * @param seed_bank 种子库
     * @param ground_light 地表光照网格
     * @param drought_index 样方干旱指数
     * @param gdd 样方生长度日
     * @param ngs_temp 非生长季温度
     * @param existing_saplings 现有幼苗（用于避免重叠）
     * @return 新定居的幼苗列表
     */
    std::vector<NewSeedling> processEstablishment(
        const PlotSeedBank& seed_bank,
        const Grid<double>& ground_light,
        double drought_index,
        double gdd,
        double ngs_temp,
        const std::vector<Sapling>& existing_saplings
    );
    
    /**
     * @brief 计算单个网格的定居概率
     * @param seed_probability 种子限制概率
     * @param env 环境条件
     * @param species 物种参数
     * @return 定居概率 [0, 1]
     */
    static double calcEstablishmentProbability(
        double seed_probability,
        const CellEnvironment& env,
        const SpeciesProfile& species
    );
    
    /**
     * @brief 为新幼苗生成随机坐标
     * @param cell_x 网格x索引
     * @param cell_y 网格y索引
     * @param occupied_positions 已占用位置列表
     * @return 随机坐标 (x, y)
     */
    std::pair<double, double> generateRandomPosition(
        int cell_x, int cell_y,
        const std::vector<std::pair<double, double>>& occupied_positions
    );

private:
    const SpeciesManager* species_mgr_ = nullptr;
    RandomGenerator* rng_ = nullptr;
    
    /**
     * @brief 检查位置是否已被占用
     */
    bool isPositionOccupied(double x, double y,
                            const std::vector<std::pair<double, double>>& occupied,
                            double min_distance = 0.1) const;
};

/**
 * @class SaplingRecruitment
 * @brief 幼树晋升处理器
 * 
 * 处理幼树晋升为成树的逻辑
 */
class SaplingRecruitment {
public:
    /**
     * @brief 检查幼树是否可以晋升
     * @param sapling 幼树对象
     * @return 是否可以晋升
     */
    static bool canRecruit(const Sapling& sapling) {
        return sapling.canRecruit();  // height >= 1.37m
    }
    
    /**
     * @brief 计算晋升时的初始胸径
     * @param height 幼树高度 (m)
     * @param species 物种参数
     * @return 初始胸径 (cm)
     * 
     * 使用反解高度-胸径关系
     */
    static double calcInitialDbh(double height, const SpeciesProfile& species);
};

} // namespace sfe

#endif // SFE_ESTABLISHMENT_HPP
