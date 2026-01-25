/**
 * @file SeedDispersal.cpp
 * @brief 种子散布模块实现
 * 
 * 亚高山森林演替模型 (SFE-Model)
 */

#include "modules/SeedDispersal.hpp"
#include "entities/AdultTree.hpp"
#include <algorithm>
#include <cmath>

namespace sfe {

// ============================================================================
// SeedDispersalCalculator 实现
// ============================================================================

double SeedDispersalCalculator::calcTreeContribution(
    const AdultTree& tree,
    const SpeciesProfile& species,
    double target_x, double target_y
) {
    if (!tree.isAlive()) return 0.0;
    
    // 检查树木是否达到性成熟
    if (tree.getAge() < species.maturity_age_yr) return 0.0;
    
    // 计算距离
    double dx = target_x - tree.getX();
    double dy = target_y - tree.getY();
    double distance = std::sqrt(dx * dx + dy * dy);
    
    // 计算核函数值
    double kernel = SeedKernel::evaluate(distance, species);
    
    // 种子产量 = F_m2 (繁殖力) × LA (叶面积)
    // F_m2 是每平方米叶面积的种子产量
    double seed_production = species.fecundity_f * tree.getLeafArea();
    
    // 到达目标点的种子密度
    return kernel * seed_production;
}

void SeedDispersalCalculator::calcInternalDispersal(
    const std::vector<AdultTree>& trees,
    const SpeciesProfile& species,
    Grid<double>& seed_grid
) {
    // 遍历每个1m网格
    for (int y = 0; y < PLOT_SIZE; ++y) {
        for (int x = 0; x < PLOT_SIZE; ++x) {
            // 网格中心坐标
            double cx = x + 0.5;
            double cy = y + 0.5;
            
            double total_density = 0.0;
            
            // 遍历所有成树
            for (const auto& tree : trees) {
                if (tree.getSpeciesId() != species.species_id) continue;
                
                total_density += calcTreeContribution(tree, species, cx, cy);
            }
            
            seed_grid.set(x, y, seed_grid.get(x, y) + total_density);
        }
    }
}

double SeedDispersalCalculator::calcExternalContribution(
    double source_strength,
    double composition_ratio,
    const SpeciesProfile& species,
    double distance
) {
    if (distance <= 0.0 || source_strength <= 0.0 || composition_ratio <= 0.0) {
        return 0.0;
    }
    
    // 核函数值
    double kernel = SeedKernel::evaluate(distance, species);
    
    // 外部种子密度 = S_external × Ratio × F_m2 × K(D) × Area_plot
    // 然后均匀分布到所有网格
    double total_seeds = source_strength * composition_ratio * 
                         species.fecundity_f * kernel * PLOT_AREA;
    
    // 每个网格的密度 (均匀分布)
    return total_seeds / PLOT_AREA;
}

// ============================================================================
// SeedDispersalEngine 实现
// ============================================================================

double SeedDispersalEngine::calcSourceStrength(const std::vector<AdultTree>& trees) {
    // S_plot = min(Σ LA_i / (Area_plot × 3.0), 1.0)
    double total_la = 0.0;
    for (const auto& tree : trees) {
        if (tree.isAlive()) {
            total_la += tree.getLeafArea();
        }
    }
    
    double lai = total_la / GlobalConfig::PLOT_AREA_M2;
    return std::min(lai / 3.0, 1.0);
}

void SeedDispersalEngine::calculateSeedRain(
    PlotId target_plot_id,
    const std::vector<AdultTree>& target_trees,
    const std::map<PlotId, std::pair<double, std::vector<const AdultTree*>>>& source_plots_info,
    PlotSeedBank& seed_bank
) {
    // 清空种子库（每年重新计算）
    seed_bank.ageSeeds();
    
    // 1. 内部散布（样方内成树产种）
    calcInternalSeedRain(target_trees, seed_bank);
    
    // 2. 来自其他模拟样方的种子
    calcInterPlotSeedRain(target_plot_id, source_plots_info, seed_bank);
    
    // 3. 来自外部种源样方的种子
    calcExternalSeedRain(target_plot_id, seed_bank);
}

void SeedDispersalEngine::calcInternalSeedRain(
    const std::vector<AdultTree>& trees,
    PlotSeedBank& seed_bank
) {
    if (!species_mgr_) return;
    
    // 对每个树种计算内部散布
    for (size_t sp = 1; sp < species_mgr_->getSpeciesCount(); ++sp) {  // 跳过草本(0)
        const auto& species = species_mgr_->getSpecies(static_cast<SpeciesId>(sp));
        
        auto& seed_grid = seed_bank.getSeedGrid(static_cast<SpeciesId>(sp));
        SeedDispersalCalculator::calcInternalDispersal(trees, species, seed_grid);
    }
}

void SeedDispersalEngine::calcInterPlotSeedRain(
    PlotId target_id,
    const std::map<PlotId, std::pair<double, std::vector<const AdultTree*>>>& source_plots,
    PlotSeedBank& seed_bank
) {
    if (!species_mgr_ || !dist_matrix_) return;
    
    constexpr int PLOT_SIZE = GlobalConfig::PLOT_SIZE_M;
    
    // 遍历所有其他模拟样方
    for (const auto& [source_id, info] : source_plots) {
        if (source_id == target_id) continue;  // 跳过自己
        
        double source_strength = info.first;
        const auto& source_trees = info.second;
        
        if (source_strength <= 0.0) continue;
        
        // 获取样方间距离
        double distance = dist_matrix_->getDistance(source_id, target_id);
        if (distance <= 0.0) continue;
        
        // 计算每个物种的贡献
        for (size_t sp = 1; sp < species_mgr_->getSpeciesCount(); ++sp) {
            const auto& species = species_mgr_->getSpecies(static_cast<SpeciesId>(sp));
            
            // 计算该物种在源样方的叶面积比例
            double species_la = 0.0;
            double total_la = 0.0;
            for (const auto* tree : source_trees) {
                if (tree && tree->isAlive()) {
                    total_la += tree->getLeafArea();
                    if (tree->getSpeciesId() == static_cast<SpeciesId>(sp)) {
                        species_la += tree->getLeafArea();
                    }
                }
            }
            
            double composition_ratio = (total_la > 0.0) ? (species_la / total_la) : 0.0;
            
            // 计算外部贡献
            double contribution = SeedDispersalCalculator::calcExternalContribution(
                source_strength, composition_ratio, species, distance
            );
            
            // 均匀添加到所有网格
            if (contribution > 0.0) {
                for (int y = 0; y < PLOT_SIZE; ++y) {
                    for (int x = 0; x < PLOT_SIZE; ++x) {
                        seed_bank.addSeeds(static_cast<SpeciesId>(sp), x, y, contribution);
                    }
                }
            }
        }
    }
}

void SeedDispersalEngine::calcExternalSeedRain(PlotId target_id, PlotSeedBank& seed_bank) {
    if (!species_mgr_ || !dist_matrix_ || !external_sources_) return;
    
    constexpr int PLOT_SIZE = GlobalConfig::PLOT_SIZE_M;
    
    // 遍历所有外部种源样方
    for (const auto& [source_id, compositions] : *external_sources_) {
        // 获取样方间距离
        double distance = dist_matrix_->getDistance(source_id, target_id);
        if (distance <= 0.0) continue;
        
        // 外部种源强度假设为饱和状态
        double source_strength = 1.0;
        
        // 对每个物种组成
        for (const auto& comp : compositions) {
            if (comp.species_id <= 0 || 
                comp.species_id >= static_cast<int>(species_mgr_->getSpeciesCount())) {
                continue;
            }
            
            const auto& species = species_mgr_->getSpecies(comp.species_id);
            
            // 计算外部贡献
            double contribution = SeedDispersalCalculator::calcExternalContribution(
                source_strength, comp.composition_ratio, species, distance
            );
            
            // 均匀添加到所有网格
            if (contribution > 0.0) {
                for (int y = 0; y < PLOT_SIZE; ++y) {
                    for (int x = 0; x < PLOT_SIZE; ++x) {
                        seed_bank.addSeeds(comp.species_id, x, y, contribution);
                    }
                }
            }
        }
    }
}

} // namespace sfe
