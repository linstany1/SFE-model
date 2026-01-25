/**
 * @file HerbLayer.hpp
 * @brief 草本层生物量管理类
 * 
 * 亚高山森林演替模型 (SFE-Model)
 * 
 * 草本层特性：
 * - 30×30 网格，每个 1m² 单元独立维护生物量
 * - 固定高度 1.0m
 * - 连续体模型（不作为个体存在）
 * - 最小生物量阈值 0.02 kg/m²
 * 
 * 关键公式（文档第6节）：
 * - LA_herb = 5.84 × B_herb [公式103]
 * - folW_herb = 0.5 × B_herb [公式104]
 * - 逻辑斯蒂生长 [公式101-102]
 */

#ifndef SFE_HERB_LAYER_HPP
#define SFE_HERB_LAYER_HPP

#include <algorithm>
#include <cmath>
#include "core/Grid.hpp"
#include "core/GlobalConfig.hpp"
#include "core/MathUtils.hpp"

namespace sfe {

/**
 * @class HerbLayer
 * @brief 草本层管理类
 * 
 * 管理样方内 30×30 个 1m² 网格单元的草本生物量。
 * 草本被视为覆盖地表的连续介质，快速占据烧后空地。
 */
class HerbLayer {
public:
    // 常量定义
    static constexpr int GRID_SIZE = GlobalConfig::PLOT_SIZE_M;  ///< 网格维度 (30)
    static constexpr double HERB_HEIGHT_M = GlobalConfig::HERB_HEIGHT_M;  ///< 固定高度 (1.0m)
    static constexpr double MIN_BIOMASS_KG = GlobalConfig::HERB_MIN_BIOMASS_KG;  ///< 最小生物量 (0.02 kg/m²)
    static constexpr double LA_COEFFICIENT = GlobalConfig::HERB_LA_COEFFICIENT;  ///< 叶面积系数 (5.84)
    static constexpr double FOLIAGE_FRACTION = GlobalConfig::HERB_FOLIAGE_FRACTION;  ///< 叶生物量比例 (0.5)
    
    // ========================================================================
    // 构造函数
    // ========================================================================
    
    /**
     * @brief 默认构造函数
     * 
     * 创建 30×30 网格，初始化为最小生物量
     */
    HerbLayer() 
        : biomass_(GRID_SIZE, GRID_SIZE, MIN_BIOMASS_KG) {
    }
    
    /**
     * @brief 使用指定初始生物量构造
     * @param initial_biomass 初始生物量 (kg/m²)
     */
    explicit HerbLayer(double initial_biomass)
        : biomass_(GRID_SIZE, GRID_SIZE, std::max(MIN_BIOMASS_KG, initial_biomass)) {
    }
    
    // ========================================================================
    // 重置与初始化
    // ========================================================================
    
    /**
     * @brief 重置为初始状态（最小生物量）
     * 
     * 用于火灾后重置
     */
    void reset() {
        biomass_.fill(MIN_BIOMASS_KG);
    }
    
    /**
     * @brief 火灾后重置
     */
    void resetAfterFire() {
        reset();
    }
    
    // ========================================================================
    // 生物量访问
    // ========================================================================
    
    /**
     * @brief 获取指定网格的生物量
     * @param x x坐标 [0, 29]
     * @param y y坐标 [0, 29]
     * @return 生物量 (kg/m²)
     */
    double getBiomass(int x, int y) const {
        return biomass_.get(x, y);
    }
    
    /**
     * @brief 设置指定网格的生物量
     * @param x x坐标
     * @param y y坐标
     * @param value 生物量 (kg/m²)，会被限制在 [MIN_BIOMASS, +∞)
     */
    void setBiomass(int x, int y, double value) {
        biomass_.set(x, y, std::max(MIN_BIOMASS_KG, value));
    }
    
    /**
     * @brief 设置所有网格的生物量为统一值
     * @param value 生物量 (kg/m²)
     */
    void setBiomass(double value) {
        double clamped = std::max(MIN_BIOMASS_KG, value);
        for (int y = 0; y < GRID_SIZE; ++y) {
            for (int x = 0; x < GRID_SIZE; ++x) {
                biomass_.set(x, y, clamped);
            }
        }
    }
    
    /**
     * @brief 获取生物量网格的只读引用
     */
    const Grid<double>& getBiomassGrid() const {
        return biomass_;
    }
    
    /**
     * @brief 获取生物量网格的可写引用
     */
    Grid<double>& getBiomassGrid() {
        return biomass_;
    }
    
    // ========================================================================
    // 叶面积计算
    // ========================================================================
    
    /**
     * @brief 计算指定网格的叶面积
     * @param x x坐标
     * @param y y坐标
     * @return 叶面积 (m²/m²，即LAI)
     * 
     * 文档公式(103): LA_herb = 5.84 × B_herb
     */
    double getLeafArea(int x, int y) const {
        return LA_COEFFICIENT * biomass_.get(x, y);
    }
    
    /**
     * @brief 计算指定网格的叶面积指数
     * 
     * 对于 1m² 网格，叶面积数值等于 LAI
     */
    double getLAI(int x, int y) const {
        return getLeafArea(x, y);
    }
    
    /**
     * @brief 计算整个样方的总叶面积
     * @return 总叶面积 (m²)
     */
    double getTotalLeafArea() const {
        return LA_COEFFICIENT * biomass_.sum();
    }
    
    /**
     * @brief 计算整个样方的平均 LAI
     */
    double getMeanLAI() const {
        return LA_COEFFICIENT * biomass_.mean();
    }
    
    // ========================================================================
    // 叶生物量计算
    // ========================================================================
    
    /**
     * @brief 计算指定网格的叶生物量
     * @param x x坐标
     * @param y y坐标
     * @return 叶生物量 (kg/m²)
     * 
     * 文档公式(104): folW_herb = 0.5 × B_herb
     */
    double getFoliageBiomass(int x, int y) const {
        return FOLIAGE_FRACTION * biomass_.get(x, y);
    }
    
    /**
     * @brief 计算整个样方的总叶生物量
     * @return 总叶生物量 (kg)
     */
    double getTotalFoliageBiomass() const {
        return FOLIAGE_FRACTION * biomass_.sum();
    }
    
    // ========================================================================
    // 统计信息
    // ========================================================================
    
    /**
     * @brief 获取总生物量
     * @return 样方总草本生物量 (kg)
     */
    double getTotalBiomass() const {
        return biomass_.sum();
    }
    
    /**
     * @brief 获取平均生物量
     * @return 平均生物量 (kg/m²)
     */
    double getMeanBiomass() const {
        return biomass_.mean();
    }
    
    /**
     * @brief 获取最大生物量
     */
    double getMaxBiomass() const {
        return biomass_.max();
    }
    
    /**
     * @brief 获取最小生物量
     */
    double getMinBiomass() const {
        return biomass_.min();
    }
    
    // ========================================================================
    // 生长模拟（基础框架，详细逻辑在模块中实现）
    // ========================================================================
    
    /**
     * @brief 对单个网格应用逻辑斯蒂生长
     * @param x x坐标
     * @param y y坐标
     * @param r_herb 生长速率参数
     * @param k_herb 环境容纳量 (kg/m²)
     * 
     * 文档公式(101-102):
     * B_{Y+1} = B_Y + r_herb × (1 - B/K_herb) × B
     * B_{Y+1} = max(B_{Y+1}, 0.02)
     */
    void applyLogisticGrowth(int x, int y, double r_herb, double k_herb) {
        double b = biomass_.get(x, y);
        
        if (k_herb > 0.0) {
            double growth = r_herb * (1.0 - b / k_herb) * b;
            b += growth;
        }
        
        // 强制最小生物量
        biomass_.set(x, y, std::max(MIN_BIOMASS_KG, b));
    }
    
    /**
     * @brief 对所有网格应用统一的生长参数
     * @param r_herb 生长速率参数
     * @param k_herb 环境容纳量 (kg/m²)
     */
    void applyGrowthAll(double r_herb, double k_herb) {
        for (int y = 0; y < GRID_SIZE; ++y) {
            for (int x = 0; x < GRID_SIZE; ++x) {
                applyLogisticGrowth(x, y, r_herb, k_herb);
            }
        }
    }
    
    /**
     * @brief 简化的生长接口
     * @param avg_light 平均光照 [0, 1]
     * @param gdd 生长度日
     * @param drought_index 干旱指数 [0, 1]
     */
    void grow(double avg_light, double gdd, double drought_index) {
        // 简化的环境响应计算
        double f_light = std::min(1.0, avg_light / 0.3);  // 草本饱和点低
        double f_temp = std::max(0.0, 1.0 - std::exp((100.0 - gdd) * 0.0027));
        double f_drought = std::sqrt(std::max(0.0, 1.0 - drought_index / 0.8));
        
        double f_env = f_light * f_temp * f_drought;
        double r_herb = 0.8 * f_env;  // 基础生长率 0.8 * 环境因子
        double k_herb = 0.5;  // 环境容纳量 0.5 kg/m²
        
        applyGrowthAll(r_herb, k_herb);
    }
    
    /**
     * @brief 强制所有网格满足最小生物量约束
     */
    void enforceMinBiomass() {
        biomass_.transform([](double val) {
            return std::max(MIN_BIOMASS_KG, val);
        });
    }
    
    // ========================================================================
    // 迭代器支持
    // ========================================================================
    
    /**
     * @brief 遍历所有网格
     * @param func 回调函数 func(x, y, biomass)
     */
    template <typename Func>
    void forEach(Func&& func) const {
        biomass_.forEach(std::forward<Func>(func));
    }
    
    /**
     * @brief 遍历并修改所有网格
     * @param func 回调函数 func(x, y, biomass&)
     */
    template <typename Func>
    void forEachMutable(Func&& func) {
        biomass_.forEach(std::forward<Func>(func));
        enforceMinBiomass();  // 确保约束
    }

private:
    Grid<double> biomass_;  ///< 生物量网格 (kg/m²)
};

} // namespace sfe

#endif // SFE_HERB_LAYER_HPP
