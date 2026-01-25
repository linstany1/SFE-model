/**
 * @file SeedDispersal.hpp
 * @brief 种子散布模块
 * 
 * 亚高山森林演替模型 (SFE-Model)
 * 
 * 实现完整的种子散布过程：
 * 1. 双指数核函数 (公式1)
 * 2. 样方内部散布 - 基于个体位置
 * 3. 样方间外部散布 - 基于距离矩阵
 * 4. 种子限制概率 (公式289)
 * 
 * 关键公式：
 * K(d) = 1/(2πd) × [(1-k_LDD)/λ1² × e^(-d/λ1) + k_LDD/λ2² × e^(-d/λ2)]
 * P_seed = min(N_arrival / N_unlimited, 1.0)
 */

#ifndef SFE_SEED_DISPERSAL_HPP
#define SFE_SEED_DISPERSAL_HPP

#include <vector>
#include <map>
#include <cmath>
#include "core/Grid.hpp"
#include "core/GlobalConfig.hpp"
#include "core/Types.hpp"
#include "io/SpeciesProfile.hpp"
#include "io/SiteParser.hpp"
#include "io/DistanceMatrixParser.hpp"

namespace sfe {

// 前向声明
class AdultTree;
class Plot;

/**
 * @struct SeedRainResult
 * @brief 单个网格的种子雨结果
 */
struct SeedRainResult {
    SpeciesId species_id = 0;
    double seed_density = 0.0;      ///< 种子密度 (seeds/m²)
    double seed_probability = 0.0;  ///< 种子限制概率 [0, 1]
};

/**
 * @class SeedKernel
 * @brief 双指数核函数计算器
 */
class SeedKernel {
public:
    /**
     * @brief 计算双指数散布核函数
     * @param distance 距离 d (m)
     * @param lambda1 短距离特征尺度 λ1 (m)
     * @param lambda2 长距离特征尺度 λ2 (m)
     * @param k_ldd 长距离散布比例 [0, 1]
     * @return 核函数值 K(d) (1/m²)
     * 
     * 公式(1): K(d) = 1/(2πd) × [(1-k_LDD)/λ1² × e^(-d/λ1) + k_LDD/λ2² × e^(-d/λ2)]
     */
    static double evaluate(double distance, double lambda1, double lambda2, double k_ldd) {
        // 避免除以0
        if (distance < 0.1) distance = 0.1;
        
        double l1_sq = lambda1 * lambda1;
        double l2_sq = lambda2 * lambda2;
        
        // 短距离分量
        double short_dist = (1.0 - k_ldd) / l1_sq * std::exp(-distance / lambda1);
        
        // 长距离分量
        double long_dist = k_ldd / l2_sq * std::exp(-distance / lambda2);
        
        // 总核函数值
        return (short_dist + long_dist) / (2.0 * 3.1416 * distance);
    }
    
    /**
     * @brief 使用物种参数计算核函数
     */
    static double evaluate(double distance, const SpeciesProfile& species) {
        return evaluate(distance, species.lambda1_m, species.lambda2_m, species.k_ldd);
    }
};

/**
 * @class SeedDispersalCalculator
 * @brief 种子散布计算器
 * 
 * 计算样方内每个1m网格的种子到达密度
 */
class SeedDispersalCalculator {
public:
    // 常量
    static constexpr int PLOT_SIZE = GlobalConfig::PLOT_SIZE_M;
    static constexpr double SEED_UNLIMITED_THRESHOLD = GlobalConfig::SEED_UNLIMITED_THRESHOLD;
    static constexpr double PLOT_AREA = GlobalConfig::PLOT_AREA_M2;
    
    SeedDispersalCalculator() = default;
    
    /**
     * @brief 计算来自单棵树的种子贡献
     * @param tree 母树
     * @param species 物种参数
     * @param target_x 目标网格中心x
     * @param target_y 目标网格中心y
     * @return 种子密度贡献 (seeds/m²)
     */
    static double calcTreeContribution(const AdultTree& tree,
                                        const SpeciesProfile& species,
                                        double target_x, double target_y);
    
    /**
     * @brief 计算样方内部散布
     * @param trees 样方内成树列表
     * @param species 物种参数
     * @param[out] seed_grid 种子密度网格 (30×30)
     */
    static void calcInternalDispersal(const std::vector<AdultTree>& trees,
                                       const SpeciesProfile& species,
                                       Grid<double>& seed_grid);
    
    /**
     * @brief 计算来自单个外部样方的种子贡献
     * @param source_strength 种源强度 [0, 1]
     * @param composition_ratio 物种组成比例 [0, 1]
     * @param species 物种参数
     * @param distance 样方间距离 (m)
     * @return 均匀分布到每个网格的种子密度 (seeds/m²)
     */
    static double calcExternalContribution(double source_strength,
                                            double composition_ratio,
                                            const SpeciesProfile& species,
                                            double distance);
    
    /**
     * @brief 计算种子限制概率
     * @param seed_density 种子密度 (seeds/m²)
     * @return 种子限制概率 [0, 1]
     * 
     * 公式: P_seed = min(N_arrival / N_unlimited, 1.0)
     */
    static double calcSeedProbability(double seed_density) {
        return std::min(seed_density / SEED_UNLIMITED_THRESHOLD, 1.0);
    }
};

/**
 * @class PlotSeedBank
 * @brief 样方种子库管理器
 * 
 * 管理每个物种在每个1m网格的种子密度
 */
class PlotSeedBank {
public:
    static constexpr int PLOT_SIZE = GlobalConfig::PLOT_SIZE_M;
    static constexpr int MAX_SEED_AGE = 2;  ///< 种子最长存活年数
    
    PlotSeedBank() = default;
    
    /**
     * @brief 初始化种子库（指定物种数量）
     */
    void initialize(int num_species) {
        num_species_ = num_species;
        // 每个物种一个网格
        seed_density_.resize(num_species);
        for (auto& grid : seed_density_) {
            grid = Grid<double>(PLOT_SIZE, PLOT_SIZE, 0.0);
        }
    }
    
    /**
     * @brief 清空种子库
     */
    void clear() {
        for (auto& grid : seed_density_) {
            grid.fill(0.0);
        }
    }
    
    /**
     * @brief 添加种子（累加到现有密度）
     */
    void addSeeds(SpeciesId species_id, int x, int y, double density) {
        if (species_id >= 0 && species_id < num_species_) {
            double current = seed_density_[species_id].get(x, y);
            seed_density_[species_id].set(x, y, current + density);
        }
    }
    
    /**
     * @brief 设置种子密度
     */
    void setSeedDensity(SpeciesId species_id, int x, int y, double density) {
        if (species_id >= 0 && species_id < num_species_) {
            seed_density_[species_id].set(x, y, density);
        }
    }
    
    /**
     * @brief 获取种子密度
     */
    double getSeedDensity(SpeciesId species_id, int x, int y) const {
        if (species_id >= 0 && species_id < num_species_) {
            return seed_density_[species_id].get(x, y);
        }
        return 0.0;
    }
    
    /**
     * @brief 获取种子限制概率
     */
    double getSeedProbability(SpeciesId species_id, int x, int y) const {
        return SeedDispersalCalculator::calcSeedProbability(
            getSeedDensity(species_id, x, y)
        );
    }
    
    /**
     * @brief 获取指定物种的种子网格
     */
    const Grid<double>& getSeedGrid(SpeciesId species_id) const {
        return seed_density_.at(species_id);
    }
    
    Grid<double>& getSeedGrid(SpeciesId species_id) {
        return seed_density_.at(species_id);
    }
    
    /**
     * @brief 年末种子库老化（移除超龄种子）
     * 简化实现：假设当年未发芽的种子全部失效
     */
    void ageSeeds() {
        // 简化处理：每年重新计算种子雨，不累积
        clear();
    }
    
    /**
     * @brief 获取物种数量
     */
    int getNumSpecies() const { return num_species_; }

private:
    int num_species_ = 0;
    std::vector<Grid<double>> seed_density_;  ///< 每个物种的种子密度网格
};

/**
 * @class SeedDispersalEngine
 * @brief 种子散布引擎
 * 
 * 整合内部和外部散布，计算完整的种子雨
 */
class SeedDispersalEngine {
public:
    SeedDispersalEngine() = default;
    
    /**
     * @brief 设置物种管理器
     */
    void setSpeciesManager(const SpeciesManager* species_mgr) {
        species_mgr_ = species_mgr;
    }
    
    /**
     * @brief 设置距离矩阵
     */
    void setDistanceMatrix(const DistanceMatrix* dist_mat) {
        dist_matrix_ = dist_mat;
    }
    
    /**
     * @brief 设置外部种源组成
     */
    void setExternalSources(const std::map<PlotId, std::vector<SeedSourceInfo>>* sources) {
        external_sources_ = sources;
    }
    
    /**
     * @brief 计算目标样方的完整种子雨
     * @param target_plot_id 目标样方ID
     * @param target_trees 目标样方内的成树
     * @param source_plots_info 其他模拟样方的信息（用于样方间散布）
     * @param[out] seed_bank 更新后的种子库
     */
    void calculateSeedRain(
        PlotId target_plot_id,
        const std::vector<AdultTree>& target_trees,
        const std::map<PlotId, std::pair<double, std::vector<const AdultTree*>>>& source_plots_info,
        PlotSeedBank& seed_bank
    );
    
    /**
     * @brief 计算样方的种源强度
     * @param trees 样方内成树
     * @return 种源强度 [0, 1]
     */
    static double calcSourceStrength(const std::vector<AdultTree>& trees);

private:
    const SpeciesManager* species_mgr_ = nullptr;
    const DistanceMatrix* dist_matrix_ = nullptr;
    const std::map<PlotId, std::vector<SeedSourceInfo>>* external_sources_ = nullptr;
    
    /**
     * @brief 计算内部散布（样方内成树产种）
     */
    void calcInternalSeedRain(const std::vector<AdultTree>& trees,
                               PlotSeedBank& seed_bank);
    
    /**
     * @brief 计算来自其他模拟样方的种子输入
     */
    void calcInterPlotSeedRain(
        PlotId target_id,
        const std::map<PlotId, std::pair<double, std::vector<const AdultTree*>>>& source_plots,
        PlotSeedBank& seed_bank
    );
    
    /**
     * @brief 计算来自外部种源样方的种子输入
     */
    void calcExternalSeedRain(PlotId target_id, PlotSeedBank& seed_bank);
};

} // namespace sfe

#endif // SFE_SEED_DISPERSAL_HPP
