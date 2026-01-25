/**
 * @file Plot.hpp
 * @brief 样方（Plot）核心容器类
 * 
 * 亚高山森林演替模型 (SFE-Model)
 * 
 * Plot 代表一个 30×30m 的模拟区域，是模型的核心空间单元。
 * 
 * 聚合的组件：
 * - SiteInfo: 样地基础信息（海拔、坡度等）
 * - Soil: 双层土壤水分模型
 * - HerbLayer: 草本生物量网格
 * - vector<AdultTree>: 成树列表
 * - vector<Sapling>: 幼树列表
 * - Grid<Cell>: 微环境网格（光照、种子等）
 */

#ifndef SFE_PLOT_HPP
#define SFE_PLOT_HPP

#include <vector>
#include <algorithm>
#include <functional>
#include <memory>

#include "core/Grid.hpp"
#include "core/GlobalConfig.hpp"
#include "core/Types.hpp"
#include "entities/SiteInfo.hpp"
#include "entities/Soil.hpp"
#include "entities/HerbLayer.hpp"
#include "entities/AdultTree.hpp"
#include "entities/Sapling.hpp"

namespace sfe {

/**
 * @struct CellState
 * @brief 单个 1×1m 网格单元的微环境状态
 */
struct CellState {
    double light_at_ground = 1.0;      ///< 地表光照强度 [0, 1] (LIF)
    double light_at_herb_top = 1.0;    ///< 草本顶部光照 [0, 1]
    double seed_density = 0.0;         ///< 种子密度 (seeds/m²)
    
    // 按物种的种子密度（可选，用于多物种种子追踪）
    std::vector<double> seed_by_species;
    
    /**
     * @brief 重置微环境状态
     */
    void reset() {
        light_at_ground = 1.0;
        light_at_herb_top = 1.0;
        seed_density = 0.0;
        seed_by_species.clear();
    }
};

/**
 * @class Plot
 * @brief 样方核心容器类
 * 
 * 管理 30×30m 区域内的所有生物和非生物组件。
 * 提供访问和管理成树、幼树、草本、土壤的接口。
 */
class Plot {
public:
    // 常量
    static constexpr int SIZE_M = GlobalConfig::PLOT_SIZE_M;  ///< 样方尺寸 (30m)
    static constexpr double AREA_M2 = GlobalConfig::PLOT_AREA_M2;  ///< 样方面积 (900m²)
    
    // ========================================================================
    // 构造函数
    // ========================================================================
    
    /**
     * @brief 默认构造函数
     */
    Plot() 
        : herb_layer_()
        , cells_(SIZE_M, SIZE_M)
    {
    }
    
    /**
     * @brief 使用样地信息构造
     * @param site_info 样地基础信息
     */
    explicit Plot(const SiteInfo& site_info)
        : site_info_(site_info)
        , soil_(site_info.getWhcMm())
        , herb_layer_()
        , cells_(SIZE_M, SIZE_M)
    {
    }
    
    /**
     * @brief 使用样地ID和WHC快速构造
     * @param plot_id 样方ID
     * @param whc_mm 土壤持水量 (mm)
     */
    Plot(PlotId plot_id, double whc_mm)
        : soil_(whc_mm)
        , herb_layer_()
        , cells_(SIZE_M, SIZE_M)
    {
        site_info_.plot_id = plot_id;
        site_info_.whc_mm = whc_mm;
    }
    
    // ========================================================================
    // 初始化与重置
    // ========================================================================
    
    /**
     * @brief 初始化样方（使用新的样地信息）
     * @param site_info 样地信息
     */
    void initialize(const SiteInfo& site_info) {
        site_info_ = site_info;
        soil_.initialize(site_info.getWhcMm());
        herb_layer_.reset();
        clearVegetation();
        resetCells();
    }
    
    /**
     * @brief 火灾后重置
     * 
     * 清除所有植被，重置土壤和草本到初始状态
     */
    void resetAfterFire() {
        // 清除所有树木
        trees_.clear();
        saplings_.clear();
        
        // 重置草本到最小生物量
        herb_layer_.resetAfterFire();
        
        // 重置土壤（保持满水状态）
        soil_.resetAfterFire();
        
        // 重置微环境
        resetCells();
        
        // 记录火灾状态
        burned_year_ = -1;  // 将由外部设置
        is_burned_ = true;
        is_in_fire_recovery_ = true;  // Phase 8: 设置恢复标记
    }
    
    /**
     * @brief Phase 8: 火灾处理方法
     * 
     * 执行火灾清理逻辑：
     * - 清空所有树木和幼苗
     * - 清空土壤种子库
     * - 重置草本为最小阈值 (0.02 kg/m²)
     * - 设置 is_in_fire_recovery = true
     * 
     * @param fire_year 火灾发生年份
     * @return 被移除的生物量 (kg)
     */
    double burn(Year fire_year) {
        double removed_biomass = 0.0;
        
        // 统计被移除的生物量
        for (const auto& tree : trees_) {
            removed_biomass += tree.getBiomass();
        }
        for (const auto& sapling : saplings_) {
            removed_biomass += sapling.getBiomass();
        }
        removed_biomass += herb_layer_.getTotalBiomass();
        
        // 记录死亡数量（用于日志）
        last_fire_trees_killed_ = static_cast<int>(trees_.size());
        last_fire_saplings_killed_ = static_cast<int>(saplings_.size());
        
        // 清空所有木本
        trees_.clear();
        saplings_.clear();
        
        // 清空种子库（重置微环境）
        resetCells();
        
        // 重置草本为最小阈值 (0.02 kg/m²)
        herb_layer_.setBiomass(0.02);
        
        // 关键状态翻转
        is_burned_ = true;
        burned_year_ = fire_year;
        is_in_fire_recovery_ = true;  // 这是启动外部种源传播的唯一触发条件
        
        return removed_biomass;
    }
    
    /**
     * @brief Phase 8: 检查是否处于火后恢复期
     * @return true 如果处于火后恢复状态
     */
    bool isInFireRecovery() const {
        return is_in_fire_recovery_;
    }
    
    /**
     * @brief Phase 8: 设置火后恢复状态
     */
    void setFireRecoveryState(bool state) {
        is_in_fire_recovery_ = state;
    }
    
    /**
     * @brief 获取上次火灾死亡的树木数量
     */
    int getLastFireTreesKilled() const { return last_fire_trees_killed_; }
    
    /**
     * @brief 获取上次火灾死亡的幼苗数量
     */
    int getLastFireSaplingsKilled() const { return last_fire_saplings_killed_; }
    
    /**
     * @brief 清除所有植被（保留土壤和草本结构）
     */
    void clearVegetation() {
        trees_.clear();
        saplings_.clear();
    }
    
    /**
     * @brief 重置所有网格单元状态
     */
    void resetCells() {
        cells_.forEach([](int, int, CellState& cell) {
            cell.reset();
        });
    }
    
    // ========================================================================
    // 样地信息访问
    // ========================================================================
    
    const SiteInfo& getSiteInfo() const { return site_info_; }
    SiteInfo& getSiteInfo() { return site_info_; }
    
    PlotId getPlotId() const { return site_info_.plot_id; }
    double getElevation() const { return site_info_.elevation_m; }
    double getLatitude() const { return site_info_.latitude; }
    
    // ========================================================================
    // 土壤访问
    // ========================================================================
    
    const Soil& getSoil() const { return soil_; }
    Soil& getSoil() { return soil_; }
    
    // ========================================================================
    // 草本层访问
    // ========================================================================
    
    const HerbLayer& getHerbLayer() const { return herb_layer_; }
    HerbLayer& getHerbLayer() { return herb_layer_; }
    
    // ========================================================================
    // 成树管理
    // ========================================================================
    
    const std::vector<AdultTree>& getTrees() const { return trees_; }
    std::vector<AdultTree>& getTrees() { return trees_; }
    
    /**
     * @brief 获取成树数量
     */
    size_t getTreeCount() const { return trees_.size(); }
    
    /**
     * @brief 获取存活成树数量
     */
    size_t getAliveTreeCount() const {
        return std::count_if(trees_.begin(), trees_.end(),
            [](const AdultTree& t) { return t.isAlive(); });
    }
    
    /**
     * @brief 添加成树
     * @param tree 成树对象
     */
    void addTree(const AdultTree& tree) {
        trees_.push_back(tree);
    }
    
    void addTree(AdultTree&& tree) {
        trees_.push_back(std::move(tree));
    }
    
    /**
     * @brief 移除死亡的成树
     * @return 移除的数量
     */
    size_t removeDeadTrees() {
        size_t before = trees_.size();
        trees_.erase(
            std::remove_if(trees_.begin(), trees_.end(),
                [](const AdultTree& t) { return !t.isAlive(); }),
            trees_.end()
        );
        return before - trees_.size();
    }
    
    /**
     * @brief 获取指定物种的成树
     * @param species_id 物种ID
     * @return 该物种的成树列表引用
     */
    std::vector<AdultTree*> getTreesBySpecies(SpeciesId species_id) {
        std::vector<AdultTree*> result;
        for (auto& tree : trees_) {
            if (tree.getSpeciesId() == species_id && tree.isAlive()) {
                result.push_back(&tree);
            }
        }
        return result;
    }
    
    // ========================================================================
    // 幼树管理
    // ========================================================================
    
    const std::vector<Sapling>& getSaplings() const { return saplings_; }
    std::vector<Sapling>& getSaplings() { return saplings_; }
    
    /**
     * @brief 获取幼树数量
     */
    size_t getSaplingCount() const { return saplings_.size(); }
    
    /**
     * @brief 获取存活幼树数量
     */
    size_t getAliveSaplingCount() const {
        return std::count_if(saplings_.begin(), saplings_.end(),
            [](const Sapling& s) { return s.isAlive(); });
    }
    
    /**
     * @brief 添加幼树
     */
    void addSapling(const Sapling& sapling) {
        saplings_.push_back(sapling);
    }
    
    void addSapling(Sapling&& sapling) {
        saplings_.push_back(std::move(sapling));
    }
    
    /**
     * @brief 移除死亡的幼树
     * @return 移除的数量
     */
    size_t removeDeadSaplings() {
        size_t before = saplings_.size();
        saplings_.erase(
            std::remove_if(saplings_.begin(), saplings_.end(),
                [](const Sapling& s) { return !s.isAlive(); }),
            saplings_.end()
        );
        return before - saplings_.size();
    }
    
    /**
     * @brief 处理幼树晋升（移动可晋升幼树到成树列表）
     * @param get_initial_dbh 函数：根据物种ID返回初始DBH
     * @return 晋升的数量
     */
    size_t processRecruitment(std::function<double(SpeciesId)> get_initial_dbh) {
        size_t recruited = 0;
        
        for (auto& sapling : saplings_) {
            if (sapling.isAlive() && sapling.canRecruit()) {
                double initial_dbh = get_initial_dbh(sapling.getSpeciesId());
                trees_.push_back(AdultTree::fromSapling(sapling, initial_dbh));
                sapling.kill();  // 标记为已处理
                ++recruited;
            }
        }
        
        // 移除已晋升的幼树
        removeDeadSaplings();
        
        return recruited;
    }
    
    // ========================================================================
    // 网格单元访问
    // ========================================================================
    
    const Grid<CellState>& getCells() const { return cells_; }
    Grid<CellState>& getCells() { return cells_; }
    
    const CellState& getCell(int x, int y) const { return cells_.get(x, y); }
    CellState& getCell(int x, int y) { return cells_.get(x, y); }
    
    // ========================================================================
    // 统计计算
    // ========================================================================
    
    /**
     * @brief 计算总叶面积指数（成树）
     */
    double calcTreeLAI() const {
        double total_la = 0.0;
        for (const auto& tree : trees_) {
            if (tree.isAlive()) {
                total_la += tree.getLeafArea();
            }
        }
        return total_la / AREA_M2;
    }
    
    /**
     * @brief 计算冠层盖度
     */
    double calcCanopyCover() const {
        double total_crown_area = 0.0;
        for (const auto& tree : trees_) {
            if (tree.isAlive()) {
                total_crown_area += tree.getCrownProjectionArea();
            }
        }
        // 简化处理：不考虑重叠
        return std::min(1.0, total_crown_area / AREA_M2);
    }
    
    /**
     * @brief 计算种源强度
     * @return 种源强度 [0, 1]
     * 
     * 文档公式: S_plot = min(LAI / 3.0, 1.0)
     */
    double calcSourceStrength() const {
        double lai = calcTreeLAI();
        return std::min(lai / GlobalConfig::LAI_SATURATION, 1.0);
    }
    
    /**
     * @brief 计算总地上生物量（成树）
     * @return 总生物量 (kg)
     */
    double calcTotalTreeBiomass() const {
        double total = 0.0;
        for (const auto& tree : trees_) {
            if (tree.isAlive()) {
                total += tree.getBiomass();
            }
        }
        return total;
    }
    
    /**
     * @brief 计算总叶生物量（所有功能群）
     */
    double calcTotalFoliageBiomass() const {
        double total = 0.0;
        
        // 成树
        for (const auto& tree : trees_) {
            if (tree.isAlive()) {
                total += tree.getFoliageBiomass();
            }
        }
        
        // 草本
        total += herb_layer_.getTotalFoliageBiomass();
        
        // 幼树（需要物种参数，此处简化忽略或由外部计算）
        
        return total;
    }
    
    // ========================================================================
    // 干扰状态
    // ========================================================================
    
    bool isBurned() const { return is_burned_; }
    void setBurned(bool burned) { is_burned_ = burned; }
    
    int getBurnedYear() const { return burned_year_; }
    void setBurnedYear(int year) { burned_year_ = year; }
    
    // ========================================================================
    // 遍历工具
    // ========================================================================
    
    /**
     * @brief 遍历所有存活成树
     */
    template <typename Func>
    void forEachAliveTree(Func&& func) {
        for (auto& tree : trees_) {
            if (tree.isAlive()) {
                func(tree);
            }
        }
    }
    
    template <typename Func>
    void forEachAliveTree(Func&& func) const {
        for (const auto& tree : trees_) {
            if (tree.isAlive()) {
                func(tree);
            }
        }
    }
    
    /**
     * @brief 遍历所有存活幼树
     */
    template <typename Func>
    void forEachAliveSapling(Func&& func) {
        for (auto& sapling : saplings_) {
            if (sapling.isAlive()) {
                func(sapling);
            }
        }
    }
    
    template <typename Func>
    void forEachAliveSapling(Func&& func) const {
        for (const auto& sapling : saplings_) {
            if (sapling.isAlive()) {
                func(sapling);
            }
        }
    }

private:
    // 样地基础信息
    SiteInfo site_info_;
    
    // 土壤水分
    Soil soil_;
    
    // 草本层
    HerbLayer herb_layer_;
    
    // 个体列表
    std::vector<AdultTree> trees_;      ///< 成树列表
    std::vector<Sapling> saplings_;     ///< 幼树列表
    
    // 微环境网格
    Grid<CellState> cells_;             ///< 30×30 网格单元
    
    // 干扰状态 (Phase 8 更新)
    bool is_burned_ = false;            ///< 是否被火烧过
    int burned_year_ = -1;              ///< 火烧年份 (-1 表示未火烧)
    bool is_in_fire_recovery_ = false;  ///< Phase 8: 是否处于火后恢复期（用于触发外部种源）
    int last_fire_trees_killed_ = 0;    ///< 上次火灾死亡树木数
    int last_fire_saplings_killed_ = 0; ///< 上次火灾死亡幼苗数
};

} // namespace sfe

#endif // SFE_PLOT_HPP
