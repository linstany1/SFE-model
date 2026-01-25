/**
 * @file AdultTree.hpp
 * @brief 成树（Adult Tree）实体类
 * 
 * 亚高山森林演替模型 (SFE-Model)
 * 
 * 成树定义：高度 >= 1.37m 的木本个体
 * 
 * 特性：
 * - 每株成树作为独立个体模拟
 * - 具有深根和高冠，主导光截获和深层水利用
 * - 通过树冠占据空间并影响下层光照
 * 
 * 关键参数：
 * - DBH: 胸径 (cm)
 * - Height: 树高 (m)
 * - Crown geometry: 冠底高度、冠幅半径、冠形
 * 
 * 关键公式：
 * - 叶面积: LA = KC2 × KA1 × DBH^KA2 [公式108]
 * - 冠幅半径: R_max = R_a × ln(DBH) + R_b [文档2.1.2]
 * - 生物量: Biomass = B0 × DBH^B1 × H^B2 [公式115-116]
 */

#ifndef SFE_ADULT_TREE_HPP
#define SFE_ADULT_TREE_HPP

#include <cmath>
#include <algorithm>
#include "core/Types.hpp"
#include "core/GlobalConfig.hpp"
#include "core/MathUtils.hpp"
#include "entities/Sapling.hpp"

namespace sfe {

/**
 * @class AdultTree
 * @brief 成树实体类
 * 
 * 存储成树的完整状态变量，包括几何参数和生长状态。
 * 复杂的生长和竞争逻辑在模块中实现。
 */
class AdultTree {
public:
    // 常量
    static constexpr double HEIGHT_THRESHOLD_M = GlobalConfig::HEIGHT_THRESHOLD_M;
    static constexpr double MIN_CROWN_FRACTION = GlobalConfig::MIN_CROWN_FRACTION;
    
    // ========================================================================
    // 构造函数
    // ========================================================================
    
    /**
     * @brief 默认构造函数
     */
    AdultTree() = default;
    
    /**
     * @brief 完整构造函数
     * @param id 唯一标识符
     * @param species_id 物种ID
     * @param x x坐标 (m)
     * @param y y坐标 (m)
     * @param dbh 胸径 (cm)
     * @param height 树高 (m)
     * @param age 年龄 (年)
     */
    AdultTree(TreeId id, SpeciesId species_id, double x, double y,
              double dbh, double height, int age = 0)
        : id_(id)
        , species_id_(species_id)
        , x_(x)
        , y_(y)
        , dbh_cm_(dbh)
        , height_m_(height)
        , age_(age)
        , stress_years_(0)
        , is_alive_(true)
    {
        // 初始化冠底高度为树高的一定比例
        crown_base_m_ = 0.0;  // 初始无冠底，需要通过物种参数计算
    }
    
    /**
     * @brief 从幼树晋升创建成树
     * @param sapling 幼树对象
     * @param initial_dbh 初始胸径 (cm)
     * @return 新成树对象
     * 
     * 晋升时继承幼树的位置、年龄等信息
     */
    static AdultTree fromSapling(const Sapling& sapling, double initial_dbh) {
        AdultTree tree(
            sapling.getId(),
            sapling.getSpeciesId(),
            sapling.getX(),
            sapling.getY(),
            initial_dbh,
            sapling.getHeight(),
            sapling.getAge()
        );
        tree.setStressYears(sapling.getStressYears());
        return tree;
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
    
    Point2D getPosition() const { return Point2D(x_, y_); }
    
    void setPosition(double x, double y) { x_ = x; y_ = y; }
    void setPosition(const Point2D& pos) { x_ = pos.x; y_ = pos.y; }
    
    CellIndex getCellIndex() const {
        return CellIndex(static_cast<int>(std::floor(x_)),
                         static_cast<int>(std::floor(y_)));
    }
    
    // ========================================================================
    // 核心生长参数
    // ========================================================================
    
    /// 胸径 (cm)
    double getDbh() const { return dbh_cm_; }
    void setDbh(double dbh) { dbh_cm_ = std::max(0.0, dbh); }
    
    /// 树高 (m)
    double getHeight() const { return height_m_; }
    void setHeight(double h) { height_m_ = std::max(HEIGHT_THRESHOLD_M, h); }
    
    /// 年龄 (年)
    int getAge() const { return age_; }
    void setAge(int age) { age_ = age; }
    void incrementAge() { ++age_; }
    
    // ========================================================================
    // 冠层几何
    // ========================================================================
    
    /// 冠底高度 (m)
    double getCrownBase() const { return crown_base_m_; }
    void setCrownBase(double hb) { 
        // 约束: 0 <= crown_base <= 0.9 * height
        crown_base_m_ = math::clamp(hb, 0.0, 0.9 * height_m_);
    }
    
    /// 最大冠幅半径 (m)
    double getCrownRadius() const { return crown_radius_m_; }
    void setCrownRadius(double r) { crown_radius_m_ = std::max(0.0, r); }
    
    /// 冠形系数
    double getCrownShape() const { return crown_shape_; }
    void setCrownShape(double cs) { crown_shape_ = cs; }
    
    /// 冠层深度 (m)
    double getCrownDepth() const { 
        return height_m_ - crown_base_m_; 
    }
    
    /**
     * @brief 计算冠层投影面积 (m²)
     */
    double getCrownProjectionArea() const {
        return math::PI * crown_radius_m_ * crown_radius_m_;
    }
    
    /**
     * @brief 计算冠层体积 (m³)
     * 
     * 文档公式: V_crown = π × R_max² × (H - HB) / (2cs + 1)
     */
    double getCrownVolume() const {
        if (crown_shape_ < 0.0) return 0.0;
        double depth = getCrownDepth();
        return math::PI * crown_radius_m_ * crown_radius_m_ * depth / (2.0 * crown_shape_ + 1.0);
    }
    
    // ========================================================================
    // 叶面积与叶密度
    // ========================================================================
    
    /// 叶面积 (m²)
    double getLeafArea() const { return leaf_area_m2_; }
    void setLeafArea(double la) { leaf_area_m2_ = std::max(0.0, la); }
    
    /**
     * @brief 计算叶面积（使用异速生长参数）
     * @param ka1 参数 KA1
     * @param ka2 参数 KA2
     * @param kc2 参数 KC2
     * @return 叶面积 (m²)
     * 
     * 文档公式(108): LA = KC2 × KA1 × DBH^KA2
     */
    double calcLeafArea(double ka1, double ka2, double kc2) const {
        return kc2 * ka1 * std::pow(dbh_cm_, ka2);
    }
    
    /**
     * @brief 计算叶面积体密度 (m²/m³)
     * 
     * 文档公式(18): ρ = LA × (2cs + 1) / (π × R_max² × (H - HB))
     */
    double getLeafAreaDensity() const {
        double volume = getCrownVolume();
        if (volume <= 0.0) return 0.0;
        return leaf_area_m2_ / volume;
    }
    
    // ========================================================================
    // 生物量
    // ========================================================================
    
    /// 地上生物量 (kg)
    double getBiomass() const { return biomass_kg_; }
    void setBiomass(double b) { biomass_kg_ = std::max(0.0, b); }
    
    /// 叶生物量 (kg)
    double getFoliageBiomass() const { return foliage_biomass_kg_; }
    void setFoliageBiomass(double fb) { foliage_biomass_kg_ = std::max(0.0, fb); }
    
    /**
     * @brief 计算地上生物量
     * @param b0, b1, b2 生物量参数
     * @return 地上生物量 (kg)
     * 
     * 文档公式(115): Biomass = B0 × DBH^B1 × H^B2 (DBH >= 5cm)
     * 文档公式(116): Biomass = b0 × DBH^b1 × H^b2 (DBH < 5cm)
     */
    double calcBiomass(double b0, double b1, double b2) const {
        return b0 * std::pow(dbh_cm_, b1) * std::pow(height_m_, b2);
    }
    
    /**
     * @brief 计算叶生物量
     * @param dw 叶片干湿比
     * @param fw1, fw2 叶重参数
     * @return 叶生物量 (kg)
     * 
     * 文档公式(109): folW = dw × fw1 × DBH^fw2
     */
    double calcFoliageBiomass(double dw, double fw1, double fw2) const {
        return dw * fw1 * std::pow(dbh_cm_, fw2);
    }
    
    // ========================================================================
    // 冠幅计算
    // ========================================================================
    
    /**
     * @brief 计算最大冠幅半径
     * @param r_a 冠幅参数 R_a
     * @param r_b 冠幅参数 R_b
     * @return 最大冠幅半径 (m)
     * 
     * 文档公式: R_max = R_a × ln(DBH) + R_b
     */
    double calcMaxCrownRadius(double r_a, double r_b) const {
        if (dbh_cm_ <= 0.0) return r_b;
        return r_a * std::log(dbh_cm_) + r_b;
    }
    
    /**
     * @brief 计算指定高度处的冠层半径
     * @param z 高度 (m)
     * @return 冠层半径 (m)
     * 
     * 文档公式: R(z) = R_max × (1 - (z-HB)/(H-HB))^cs
     */
    double calcCrownRadiusAtHeight(double z) const {
        if (z < crown_base_m_ || z > height_m_) return 0.0;
        
        double depth = height_m_ - crown_base_m_;
        if (depth <= 0.0) return 0.0;
        
        double relative = (z - crown_base_m_) / depth;
        return crown_radius_m_ * std::pow(1.0 - relative, crown_shape_);
    }
    
    // ========================================================================
    // 胁迫追踪
    // ========================================================================
    
    int getStressYears() const { return stress_years_; }
    void setStressYears(int years) { stress_years_ = years; }
    void incrementStressYears() { ++stress_years_; }
    void resetStressYears() { stress_years_ = 0; }
    
    // ========================================================================
    // 生存状态
    // ========================================================================
    
    bool isAlive() const { return is_alive_; }
    void setAlive(bool alive) { is_alive_ = alive; }
    void kill() { is_alive_ = false; }
    
    // ========================================================================
    // 光照相关
    // ========================================================================
    
    /**
     * @brief 计算有效受光高度
     * @return 有效受光高度 (m)
     * 
     * 文档公式(28): Z_eff = H - 0.25 × (H - HB)
     * 假设光合作用集中在冠层上部75%区域
     */
    double calcEffectiveLightHeight() const {
        return height_m_ - 0.25 * (height_m_ - crown_base_m_);
    }
    
    /**
     * @brief 获取缓存的可用光照
     */
    double getAvailableLight() const { return available_light_; }
    void setAvailableLight(double al) { available_light_ = math::clamp01(al); }
    
    // ========================================================================
    // 生长计算框架
    // ========================================================================
    
    /**
     * @brief 计算最优胸径增量
     * @param dbh_max 最大胸径 (cm)
     * @param h_max 最大树高 (m)
     * @param ga 生长参数 ga
     * @param opt_s 异速生长参数 opt_s
     * @return 最优胸径增量 (cm)
     * 
     * 文档公式(105):
     * ΔDBH_opt = ga × DBH × (1 - DBH×H/(DBH_max×H_max)) / 
     *            (2H + opt_s × DBH × e^(-opt_s×DBH/(H_max-1.37)))
     */
    double calcOptimalDbhIncrement(double dbh_max, double h_max, 
                                    double ga, double opt_s) const {
        if (dbh_max <= 0.0 || h_max <= HEIGHT_THRESHOLD_M) return 0.0;
        
        double size_ratio = (dbh_cm_ * height_m_) / (dbh_max * h_max);
        double numerator = ga * dbh_cm_ * (1.0 - size_ratio);
        
        double exp_term = std::exp(-opt_s * dbh_cm_ / (h_max - HEIGHT_THRESHOLD_M));
        double denominator = 2.0 * height_m_ + opt_s * dbh_cm_ * exp_term;
        
        if (denominator <= 0.0) return 0.0;
        return numerator / denominator;
    }
    
    /**
     * @brief 计算成树高度（基于DBH的Mitscherlich模型）
     * @param h_max 最大树高 (m)
     * @param opt_s 异速生长参数
     * @param dbh 胸径 (cm)
     * @return 树高 (m)
     * 
     * 文档公式(106) 下半部分:
     * H = 1.37 + (H_max - 1.37) × (1 - e^(-opt_s × DBH / (H_max - 1.37)))
     */
    static double calcHeightFromDbh(double h_max, double opt_s, double dbh) {
        if (h_max <= HEIGHT_THRESHOLD_M) return HEIGHT_THRESHOLD_M;
        
        double h_range = h_max - HEIGHT_THRESHOLD_M;
        double exp_term = std::exp(-opt_s * dbh / h_range);
        return HEIGHT_THRESHOLD_M + h_range * (1.0 - exp_term);
    }
    
    // ========================================================================
    // 更新冠层几何（需要物种参数）
    // ========================================================================
    
    /**
     * @brief 更新所有冠层几何参数
     * @param r_a 冠幅参数 R_a
     * @param r_b 冠幅参数 R_b
     * @param cs 冠形系数
     * @param alpha_base 开阔地冠底比例
     * @param ka1, ka2, kc2 叶面积参数
     */
    void updateCrownGeometry(double r_a, double r_b, double cs,
                             double alpha_base,
                             double ka1, double ka2, double kc2) {
        crown_shape_ = cs;
        crown_radius_m_ = calcMaxCrownRadius(r_a, r_b);
        leaf_area_m2_ = calcLeafArea(ka1, ka2, kc2);
        
        // 开阔地条件下的默认冠底高度
        double hb_open = alpha_base * height_m_;
        crown_base_m_ = std::max(crown_base_m_, hb_open);
        
        // 确保冠底不超过树高的90%
        crown_base_m_ = std::min(crown_base_m_, 0.9 * height_m_);
    }

private:
    // 标识信息
    TreeId id_ = 0;
    SpeciesId species_id_ = 0;
    
    // 空间位置
    double x_ = 0.0;
    double y_ = 0.0;
    
    // 核心生长参数
    double dbh_cm_ = 0.0;           ///< 胸径 (cm)
    double height_m_ = HEIGHT_THRESHOLD_M;  ///< 树高 (m)
    int age_ = 0;                   ///< 年龄 (年)
    
    // 冠层几何
    double crown_base_m_ = 0.0;     ///< 冠底高度 (m)
    double crown_radius_m_ = 0.0;   ///< 最大冠幅半径 (m)
    double crown_shape_ = 1.0;      ///< 冠形系数 cs
    
    // 叶面积与生物量
    double leaf_area_m2_ = 0.0;     ///< 叶面积 (m²)
    double biomass_kg_ = 0.0;       ///< 地上生物量 (kg)
    double foliage_biomass_kg_ = 0.0; ///< 叶生物量 (kg)
    
    // 环境状态缓存
    double available_light_ = 1.0;  ///< 可用光照 [0, 1]
    
    // 胁迫追踪
    int stress_years_ = 0;          ///< 连续胁迫年数
    
    // 生存状态
    bool is_alive_ = true;
};

} // namespace sfe

#endif // SFE_ADULT_TREE_HPP
