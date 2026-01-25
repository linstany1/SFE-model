/**
 * @file LightEngine.hpp
 * @brief 光照计算引擎 (Phase 9 重构版)
 * 
 * 亚高山森林演替模型 (SFE-Model)
 * 
 * Phase 9 架构升级：
 * - 双层空间加速结构：
 *   1. SpatialHash (5×5m)：用于成树互遮阴计算 (AL_tree, AL_HB, AL_top)
 *   2. CanopyProjectionMap (1×1m)：用于近地层光照加速 (AL_herb, AL_sapling)
 * - 幼树分层光照计算：
 *   - H >= 1.0m：直接使用 calculateLightAtPoint(x, y, H)
 *   - H < 1.0m：受草本双重遮阴
 * - 动态冠底高度 (DCBH) 公式更新
 * 
 * 关键公式（符合文档 Ch 2）：
 * - 遮阴路径长度 L_path
 * - Beer-Lambert 光照衰减
 * - 有效受光高度 Z_eff = H - 0.25(H - HB)
 * - 冠底高度对数插值公式
 */

#ifndef SFE_LIGHT_ENGINE_HPP
#define SFE_LIGHT_ENGINE_HPP

#include <vector>
#include <cmath>
#include "core/Grid.hpp"
#include "core/GlobalConfig.hpp"
#include "core/MathUtils.hpp"
#include "modules/SpatialHash.hpp"

namespace sfe {

// 前向声明
class AdultTree;
class Sapling;
class Plot;

/**
 * @struct LightCalculationParams
 * @brief 光照计算参数
 */
struct LightCalculationParams {
    double extinction_ke = GlobalConfig::EXTINCTION_KE_DEFAULT;  ///< 消光系数
    double herb_height = GlobalConfig::HERB_HEIGHT_M;            ///< 草本层高度 (m)
    double seedling_height = GlobalConfig::SEEDLING_INIT_HEIGHT_M; ///< 幼苗初始高度 (m)
};

/**
 * @class LightEngine
 * @brief 光照计算引擎
 * 
 * Phase 9 架构：双层空间加速
 * - SpatialHash: 用于成树层光照计算
 * - CanopyProjectionMap: 用于近地层光照加速
 */
class LightEngine {
public:
    // 常量
    static constexpr int PLOT_SIZE = GlobalConfig::PLOT_SIZE_M;
    static constexpr double EFFECTIVE_CROWN_FRACTION = GlobalConfig::EFFECTIVE_LIGHT_FRACTION;
    
    /**
     * @brief 默认构造函数
     */
    LightEngine() 
        : ground_light_(PLOT_SIZE, PLOT_SIZE, 1.0)
        , herb_top_light_(PLOT_SIZE, PLOT_SIZE, 1.0)
        , canopy_map_(PLOT_SIZE, PLOT_SIZE)
        , trees_ptr_(nullptr)
    {
        // 初始化 canopy_map_ 为空向量
        for (int y = 0; y < PLOT_SIZE; ++y) {
            for (int x = 0; x < PLOT_SIZE; ++x) {
                canopy_map_.set(x, y, std::vector<int>());
            }
        }
    }
    
    /**
     * @brief 设置计算参数
     */
    void setParams(const LightCalculationParams& params) {
        params_ = params;
    }
    
    /**
     * @brief 重建空间索引（空间哈希）
     * @param trees 成树列表
     */
    void rebuildSpatialIndex(std::vector<AdultTree>& trees) {
        spatial_hash_.insertAll(trees);
        updateMaxCrownRadius(trees);
    }
    
    /**
     * @brief 构建冠层投影位图
     * @param trees 成树列表
     * 
     * 在每个模拟年初调用，将每棵成树的 ID 注册到其冠幅投影覆盖的所有 1×1m 网格中
     * 用于加速近地层光照计算（避免遍历所有成树）
     */
    void buildCanopyProjectionMap(const std::vector<AdultTree>& trees);
    
    /**
     * @brief 计算样方内所有光照
     * @param trees 成树列表
     * @param herb_lai_grid 草本LAI网格（可选）
     * 
     * 执行顺序：
     * 1. 重建空间哈希（成树互遮阴）
     * 2. 构建冠层投影位图（近地层加速）
     * 3. 计算成树光照 (AL_tree)
     * 4. 计算草本顶部光照 (AL_top_herb) - 使用 CanopyProjectionMap
     * 5. 计算地表光照（考虑草本遮挡）
     */
    void calculateAllLight(std::vector<AdultTree>& trees,
                           const Grid<double>* herb_lai_grid = nullptr);
    
    /**
     * @brief 计算单棵成树的可用光照
     * @param tree 目标树
     * @return 可用光照 AL_tree [0, 1]
     * 
     * 使用有效受光高度 Z_eff = H - 0.25(H - HB)
     */
    double calculateTreeLight(const AdultTree& tree) const;
    
    /**
     * @brief 计算指定点指定高度的可用光照
     * @param x x坐标 (m)
     * @param y y坐标 (m)
     * @param z 高度 (m)
     * @return 可用光照 [0, 1]
     * 
     * 使用 SpatialHash 加速邻域查询
     */
    double calculateLightAtPoint(double x, double y, double z) const;
    
    /**
     * @brief 计算幼树的可用光照 (Phase 9 新增)
     * @param sapling 幼树对象
     * @param herb_lai 所在网格的草本叶面积指数
     * @return 幼树可用光照 AL_sapling [0, 1]
     * 
     * 分层逻辑：
     * - H >= 1.0m：直接调用 calculateLightAtPoint(x, y, H)，仅受成树遮阴
     * - H < 1.0m：受成树+草本双重遮阴
     *   AL = AL_top_herb × exp(-k_e × LA_herb × (1.0 - H) / 1.0)
     */
    double calcSaplingLight(const Sapling& sapling, double herb_lai) const;
    
    /**
     * @brief 获取地表光照网格
     */
    const Grid<double>& getGroundLight() const { return ground_light_; }
    Grid<double>& getGroundLight() { return ground_light_; }
    
    /**
     * @brief 获取草本顶部光照网格
     */
    const Grid<double>& getHerbTopLight() const { return herb_top_light_; }
    
    /**
     * @brief 获取指定网格的地表光照
     */
    double getGroundLightAt(int x, int y) const {
        return ground_light_.get(x, y);
    }
    
    /**
     * @brief 获取指定网格的草本顶部光照
     */
    double getHerbTopLightAt(int x, int y) const {
        return herb_top_light_.get(x, y);
    }
    
    /**
     * @brief 获取空间哈希（只读）
     */
    const SpatialHash& getSpatialHash() const { return spatial_hash_; }
    
    /**
     * @brief 获取冠层投影位图（只读）
     */
    const Grid<std::vector<int>>& getCanopyMap() const { return canopy_map_; }
    
    // ========================================================================
    // 核心光照计算函数
    // ========================================================================
    
    /**
     * @brief 计算冠层进入高度
     * @param horizontal_dist 水平距离 d_ij (m)
     * @param crown_base 冠底高度 HB (m)
     * @param tree_height 树高 H (m)
     * @param crown_radius 最大冠幅半径 R_max (m)
     * @param crown_shape 冠形系数 cs
     * @return 进入高度 Z_in (m)
     * 
     * 公式: Z_in = HB + (H - HB) × [1 - (d/R_max)^(1/cs)]
     */
    static double calcCrownEntryHeight(double horizontal_dist,
                                        double crown_base, double tree_height,
                                        double crown_radius, double crown_shape);
    
    /**
     * @brief 计算遮阴路径长度
     * @param target_x 目标点x
     * @param target_y 目标点y
     * @param target_z 目标点高度z
     * @param shading_tree 遮阴树
     * @return 路径长度 L_path (m)
     * 
     * 公式: L_path = max(0, Z_in - max(z, HB))
     */
    static double calcShadingPathLength(double target_x, double target_y, double target_z,
                                         const AdultTree& shading_tree);
    
    /**
     * @brief 计算光照衰减因子
     * @param leaf_density 叶面积密度 ρ (m²/m³)
     * @param path_length 路径长度 L_path (m)
     * @param extinction_ke 消光系数
     * @return 衰减因子 (0, 1]
     * 
     * Beer-Lambert: exp(-ke × ρ × L_path)
     */
    static double calcLightAttenuation(double leaf_density, double path_length,
                                        double extinction_ke);

private:
    /**
     * @brief 计算单个网格点的草本顶部光照（使用 CanopyProjectionMap）
     * @param cell_x 网格x索引
     * @param cell_y 网格y索引
     * @return 草本顶部光照 AL_top_herb [0, 1]
     * 
     * Phase 9：使用 CanopyProjectionMap 加速，仅遍历覆盖此网格的树木
     */
    double calcGroundLightAtCell(int cell_x, int cell_y) const;
    
    /**
     * @brief 更新最大冠幅半径
     */
    void updateMaxCrownRadius(const std::vector<AdultTree>& trees);

private:
    LightCalculationParams params_;           ///< 计算参数
    SpatialHash spatial_hash_;                ///< 空间哈希（成树互遮阴）
    Grid<double> ground_light_;               ///< 地表光照网格 (H=0.05m)
    Grid<double> herb_top_light_;             ///< 草本顶部光照网格 (H=1.0m)
    Grid<std::vector<int>> canopy_map_;       ///< 冠层投影位图（Phase 9新增）
    const std::vector<AdultTree>* trees_ptr_; ///< 成树列表指针
    double max_crown_radius_ = 5.0;           ///< 当前最大冠幅半径
};

/**
 * @struct TreeLightInfo
 * @brief 树木光照信息
 */
struct TreeLightInfo {
    TreeId tree_id = 0;
    double available_light = 1.0;        ///< 可用光照 AL_tree [0, 1]
    double crown_base_light = 1.0;       ///< 冠底处光照 AL_HB
    double crown_top_light = 1.0;        ///< 冠顶处光照 AL_top
};

/**
 * @class CrownGeometryCalculator
 * @brief 冠层几何计算器
 */
class CrownGeometryCalculator {
public:
    /**
     * @brief 更新树木的冠层几何参数
     * @param tree 树木对象
     * @param h_max 物种最大树高
     * @param opt_s 异速生长参数
     * @param r_a 冠幅参数 R_a
     * @param r_b 冠幅参数 R_b
     * @param cs 冠形系数
     * @param alpha_base 开阔地冠底比例
     * @param ka1, ka2, kc2 叶面积参数
     */
    static void updateGeometry(AdultTree& tree,
                                double h_max, double opt_s,
                                double r_a, double r_b, double cs,
                                double alpha_base,
                                double ka1, double ka2, double kc2);
    
    /**
     * @brief 更新冠底高度（Phase 9 重构）
     * @param tree 树木对象
     * @param al_hb 冠底处光照 AL_HB
     * @param al_top 冠顶处光照 AL_top
     * @param l_min 物种光补偿点 L_min,s
     * @param alpha_base 开阔地冠底比例
     * 
     * 公式（文档 Ch 2.1.1）：
     * z* = HB + (H - HB) × [ln(L_min) - ln(AL_HB)] / [ln(AL_top) - ln(AL_HB)]
     * HB_new = min(0.9H, max(HB_old, HB_open, z*))
     * 
     * 修枝条件：AL_HB < L_min < AL_top
     */
    static void updateCrownBase(AdultTree& tree,
                                 double al_hb,
                                 double al_top,
                                 double l_min,
                                 double alpha_base);
};

} // namespace sfe

#endif // SFE_LIGHT_ENGINE_HPP
