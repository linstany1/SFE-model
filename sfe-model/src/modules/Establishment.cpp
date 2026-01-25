/**
 * @file Establishment.cpp
 * @brief 幼树定居模块实现
 * 
 * 亚高山森林演替模型 (SFE-Model)
 */

#include "modules/Establishment.hpp"
#include "modules/SeedDispersal.hpp"
#include <algorithm>
#include <cmath>
#include <set>

namespace sfe {

// ============================================================================
// EstablishmentEngine 实现
// ============================================================================

double EstablishmentEngine::calcEstablishmentProbability(
    double seed_probability,
    const CellEnvironment& env,
    const SpeciesProfile& species
) {
    // P_est = P_seed × (f_T_NGSes × f_GDDes × f_Dres × f_ALes)
    auto filters = EstablishmentFilter::evaluate(env, species);
    return seed_probability * filters.getEnvironmentProbability();
}

std::vector<NewSeedling> EstablishmentEngine::processEstablishment(
    const PlotSeedBank& seed_bank,
    const Grid<double>& ground_light,
    double drought_index,
    double gdd,
    double ngs_temp,
    const std::vector<Sapling>& existing_saplings
) {
    std::vector<NewSeedling> new_seedlings;
    
    if (!species_mgr_ || !rng_) return new_seedlings;
    
    // 构建现有幼苗位置索引（按网格）
    std::map<std::pair<int, int>, std::vector<std::pair<double, double>>> cell_occupancy;
    for (const auto& sapling : existing_saplings) {
        if (sapling.isAlive()) {
            int cx = static_cast<int>(sapling.getX());
            int cy = static_cast<int>(sapling.getY());
            cell_occupancy[{cx, cy}].emplace_back(sapling.getX(), sapling.getY());
        }
    }
    
    // 遍历每个1m网格
    for (int y = 0; y < PLOT_SIZE; ++y) {
        for (int x = 0; x < PLOT_SIZE; ++x) {
            // 获取网格环境
            CellEnvironment env;
            env.light_at_ground = ground_light.get(x, y);
            env.drought_index = drought_index;
            env.gdd = gdd;
            env.ngs_temp = ngs_temp;
            
            // 本网格已定居的物种
            std::set<SpeciesId> established_species;
            
            // 本网格已占用的位置
            auto& occupied = cell_occupancy[{x, y}];
            
            // 对每个物种检查定居
            for (int sp = 1; sp < seed_bank.getNumSpecies(); ++sp) {
                // 每个网格每个物种最多定居一株
                if (established_species.count(static_cast<SpeciesId>(sp)) > 0) {
                    continue;
                }
                
                const auto& species = species_mgr_->getSpecies(static_cast<SpeciesId>(sp));
                
                // 获取种子限制概率
                double seed_prob = seed_bank.getSeedProbability(static_cast<SpeciesId>(sp), x, y);
                if (seed_prob <= 0.0) continue;
                
                // 计算定居概率
                double est_prob = calcEstablishmentProbability(seed_prob, env, species);
                if (est_prob <= 0.0) continue;
                
                // 随机判定
                double rand_val = rng_->getUniform01();
                if (rand_val < est_prob) {
                    // 定居成功！
                    NewSeedling seedling;
                    seedling.species_id = static_cast<SpeciesId>(sp);
                    seedling.cell_x = x;
                    seedling.cell_y = y;
                    
                    // 生成随机坐标
                    auto pos = generateRandomPosition(x, y, occupied);
                    seedling.x = pos.first;
                    seedling.y = pos.second;
                    
                    new_seedlings.push_back(seedling);
                    established_species.insert(static_cast<SpeciesId>(sp));
                    occupied.push_back(pos);
                }
            }
        }
    }
    
    return new_seedlings;
}

std::pair<double, double> EstablishmentEngine::generateRandomPosition(
    int cell_x, int cell_y,
    const std::vector<std::pair<double, double>>& occupied_positions
) {
    // 尝试多次生成不重叠的位置
    constexpr int MAX_ATTEMPTS = 10;
    constexpr double MIN_DISTANCE = 0.05;  // 最小间距 5cm
    
    for (int attempt = 0; attempt < MAX_ATTEMPTS; ++attempt) {
        double x = cell_x + rng_->getUniform01();
        double y = cell_y + rng_->getUniform01();
        
        if (!isPositionOccupied(x, y, occupied_positions, MIN_DISTANCE)) {
            return {x, y};
        }
    }
    
    // 如果尝试失败，返回网格中心附近
    return {cell_x + 0.5, cell_y + 0.5};
}

bool EstablishmentEngine::isPositionOccupied(
    double x, double y,
    const std::vector<std::pair<double, double>>& occupied,
    double min_distance
) const {
    for (const auto& pos : occupied) {
        double dx = x - pos.first;
        double dy = y - pos.second;
        if (dx * dx + dy * dy < min_distance * min_distance) {
            return true;
        }
    }
    return false;
}

// ============================================================================
// SaplingRecruitment 实现
// ============================================================================

double SaplingRecruitment::calcInitialDbh(double height, const SpeciesProfile& species) {
    // 使用高度-胸径关系反解
    // H = 1.37 + (H_max - 1.37) × (1 - e^(-opt_s × DBH / (H_max - 1.37)))
    // 反解得：DBH = -(H_max - 1.37) / opt_s × ln(1 - (H - 1.37) / (H_max - 1.37))
    
    double h_max = species.h_max_m;
    double opt_s = species.opt_s;
    
    if (h_max <= 1.37 || opt_s <= 0.0) {
        return 1.0;  // 默认初始胸径
    }
    
    double h_range = h_max - 1.37;
    double height_ratio = (height - 1.37) / h_range;
    
    // 限制范围避免数值问题
    height_ratio = std::min(0.99, std::max(0.01, height_ratio));
    
    double dbh = -(h_range / opt_s) * std::log(1.0 - height_ratio);
    
    // 确保初始胸径在合理范围
    return std::max(0.5, std::min(dbh, 5.0));  // 0.5-5.0 cm
}

} // namespace sfe
