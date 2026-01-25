/**
 * @file Sapling.hpp
 * @brief 幼树（Sapling）实体类
 * 
 * 亚高山森林演替模型 (SFE-Model)
 * 
 * 幼树定义：高度 < 1.37m 的木本个体
 * 
 * 特性：
 * - 每株幼树作为独立个体模拟
 * - 根浅，易受干旱和庇荫影响
 * - 晋升条件：height >= 1.37m
 * 
 * 关键公式：
 * - 初始高度：5cm (0.05m)
 * - Bertalanffy 生长模型 [公式86]
 * - 生物量：sa × H^sb [公式90]
 * - 叶生物量：常绿40%，落叶20% [公式94-95]
 */

#ifndef SFE_SAPLING_HPP
#define SFE_SAPLING_HPP

#include <cmath>
#include <algorithm>
#include "core/Types.hpp"
#include "core/GlobalConfig.hpp"

namespace sfe {

/**
 * @class Sapling
 * @brief 幼树实体类
 * 
 * 存储幼树的状态变量，包括空间位置、生长状态等。
 * 复杂的生长逻辑在模块中实现。
 */
class Sapling {
public:
    // 常量
    static constexpr double HEIGHT_THRESHOLD_M = GlobalConfig::HEIGHT_THRESHOLD_M;  ///< 晋升阈值 (1.37m)
    static constexpr double INIT_HEIGHT_M = GlobalConfig::SEEDLING_INIT_HEIGHT_M;   ///< 初始高度 (0.05m)
    
    // ========================================================================
    // 构造函数
    // ========================================================================
    
    /**
     * @brief 默认构造函数
     */
    Sapling() = default;
    
    /**
     * @brief 完整构造函数
     * @param id 唯一标识符
     * @param species_id 物种ID
     * @param x x坐标 (m)
     * @param y y坐标 (m)
     * @param height 初始高度 (m)，默认 0.05m
     */
    Sapling(TreeId id, SpeciesId species_id, double x, double y, 
            double height = INIT_HEIGHT_M)
        : id_(id)
        , species_id_(species_id)
        , x_(x)
        , y_(y)
        , height_(std::max(INIT_HEIGHT_M, height))
        , age_(0)
        , stress_years_(0)
        , is_alive_(true)
    {
    }
    
    /**
     * @brief 工厂方法：创建新定居的幼苗
     * @param species_id 物种ID
     * @param x x坐标
     * @param y y坐标
     * @return 新幼苗对象
     */
    static Sapling createNewSeedling(SpeciesId species_id, double x, double y) {
        return Sapling(generateTreeId(), species_id, x, y, INIT_HEIGHT_M);
    }
    
    // ========================================================================
    // ID 和物种
    // ========================================================================
    
    TreeId getId() const { return id_; }
    void setId(TreeId id) { id_ = id; }
    
    SpeciesId getSpeciesId() const { return species_id_; }
    void setSpeciesId(SpeciesId id) { species_id_ = id; }
    
    // ========================================================================
    // 空间位置
    // ========================================================================
    
    double getX() const { return x_; }
    void setX(double x) { x_ = x; }
    
    double getY() const { return y_; }
    void setY(double y) { y_ = y; }
    
    /**
     * @brief 获取位置点
     */
    Point2D getPosition() const { return Point2D(x_, y_); }
    
    /**
     * @brief 设置位置
     */
    void setPosition(double x, double y) { x_ = x; y_ = y; }
    void setPosition(const Point2D& pos) { x_ = pos.x; y_ = pos.y; }
    
    /**
     * @brief 获取所在网格索引
     */
    CellIndex getCellIndex() const {
        return CellIndex(static_cast<int>(std::floor(x_)),
                         static_cast<int>(std::floor(y_)));
    }
    
    // ========================================================================
    // 生长状态
    // ========================================================================
    
    double getHeight() const { return height_; }
    void setHeight(double h) { height_ = std::max(0.0, h); }
    
    int getAge() const { return age_; }
    void setAge(int age) { age_ = age; }
    
    /**
     * @brief 增加年龄（每年调用一次）
     */
    void incrementAge() { ++age_; }
    
    // ========================================================================
    // 胁迫追踪（用于死亡率计算）
    // ========================================================================
    
    int getStressYears() const { return stress_years_; }
    void setStressYears(int years) { stress_years_ = years; }
    
    /**
     * @brief 增加胁迫年数
     */
    void incrementStressYears() { ++stress_years_; }
    
    /**
     * @brief 重置胁迫年数（当条件改善时）
     */
    void resetStressYears() { stress_years_ = 0; }
    
    // ========================================================================
    // 生存状态
    // ========================================================================
    
    bool isAlive() const { return is_alive_; }
    void setAlive(bool alive) { is_alive_ = alive; }
    
    /**
     * @brief 标记为死亡
     */
    void kill() { is_alive_ = false; }
    
    /**
     * @brief 执行高度生长
     * @param delta_height 高度增量 (m)
     */
    void grow(double delta_height) {
        height_ = std::max(INIT_HEIGHT_M, height_ + delta_height);
    }
    
    // ========================================================================
    // 晋升判定
    // ========================================================================
    
    /**
     * @brief 检查是否可以晋升为成树
     * @return true 如果 height >= 1.37m
     */
    bool canRecruit() const {
        return height_ >= HEIGHT_THRESHOLD_M;
    }
    
    /**
     * @brief 同上（语义别名）
     */
    bool isRecruited() const {
        return canRecruit();
    }
    
    // ========================================================================
    // 生物量计算（需要物种参数）
    // ========================================================================
    
    /**
     * @brief 计算地上生物量
     * @param sa 幼树生物量参数 sa
     * @param sb 幼树生物量参数 sb
     * @return 地上生物量 (kg)
     * 
     * 文档公式(90): Biomass_sapling = sa × H^sb
     */
    double calcBiomass(double sa, double sb) const {
        return sa * std::pow(height_, sb);
    }
    
    /**
     * @brief 获取缓存的生物量（简化版）
     * 使用通用参数估算幼树生物量
     */
    double getBiomass() const {
        // 使用通用经验公式：Biomass ≈ 0.1 * H^2.5
        return 0.1 * std::pow(height_, 2.5);
    }
    
    void setBiomass(double b) { biomass_ = b; }
    double getStoredBiomass() const { return biomass_; }
    
    /**
     * @brief 计算叶生物量
     * @param sa 幼树生物量参数 sa
     * @param sb 幼树生物量参数 sb
     * @param is_evergreen 是否常绿
     * @return 叶生物量 (kg)
     * 
     * 文档公式(94-95):
     * - 常绿: folW = 0.4 × Biomass
     * - 落叶: folW = 0.2 × Biomass
     */
    double calcFoliageBiomass(double sa, double sb, bool is_evergreen) const {
        double biomass = calcBiomass(sa, sb);
        double fraction = is_evergreen ? 
            GlobalConfig::SAPLING_FOLIAGE_EVERGREEN : 
            GlobalConfig::SAPLING_FOLIAGE_DECIDUOUS;
        return fraction * biomass;
    }
    
    // ========================================================================
    // 高度生长（Bertalanffy 模型框架）
    // ========================================================================
    
    /**
     * @brief 计算潜在高度（不含环境限制）
     * @param h_max 物种最大高度 (m)
     * @param gs 生长速率参数
     * @return 下一年的潜在高度 (m)
     * 
     * 文档公式(86):
     * H_{t+1} = H_max × [1 - (1 - (H_t/H_max)^(1/3)) × e^(-gs)]^3
     */
    double calcPotentialHeight(double h_max, double gs) const {
        if (h_max <= 0.0) return height_;
        
        double ratio = height_ / h_max;
        double inner = 1.0 - std::cbrt(ratio);  // (H/H_max)^(1/3)
        double growth_factor = 1.0 - inner * std::exp(-gs);
        
        return h_max * std::pow(growth_factor, 3.0);
    }
    
    /**
     * @brief 计算潜在高度增量
     * @param h_max 物种最大高度 (m)
     * @param gs 生长速率参数
     * @return 潜在高度增量 (m)
     */
    double calcPotentialHeightIncrement(double h_max, double gs) const {
        return calcPotentialHeight(h_max, gs) - height_;
    }

private:
    TreeId id_ = 0;                ///< 唯一标识符
    SpeciesId species_id_ = 0;    ///< 物种ID
    double x_ = 0.0;              ///< x坐标 (m)
    double y_ = 0.0;              ///< y坐标 (m)
    double height_ = INIT_HEIGHT_M;  ///< 高度 (m)
    double biomass_ = 0.0;        ///< 缓存的生物量 (kg)
    int age_ = 0;                 ///< 年龄 (年)
    int stress_years_ = 0;        ///< 连续胁迫年数
    bool is_alive_ = true;        ///< 生存状态
};

} // namespace sfe

#endif // SFE_SAPLING_HPP
