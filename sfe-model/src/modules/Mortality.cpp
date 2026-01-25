/**
 * @file Mortality.cpp
 * @brief 死亡率计算模块实现
 * 
 * 亚高山森林演替模型 (SFE-Model)
 */

#include "modules/Mortality.hpp"
#include "entities/AdultTree.hpp"
#include "entities/Sapling.hpp"
#include <algorithm>

namespace sfe {

// ============================================================================
// MortalityEngine 实现
// ============================================================================

int MortalityEngine::processTreeMortality(
    std::vector<AdultTree>& trees,
    const std::vector<GrowthResponse>& growth_responses
) {
    if (!species_mgr_ || !rng_) return 0;
    
    int death_count = 0;
    
    for (size_t i = 0; i < trees.size(); ++i) {
        auto& tree = trees[i];
        if (!tree.isAlive()) continue;
        
        // 获取生长响应
        double f_env = 1.0;
        if (i < growth_responses.size()) {
            f_env = growth_responses[i].getEnvironmentResponse();
        }
        
        // 更新胁迫年数（使用 StressTracker 的返回值）
        int new_stress_years = StressTracker::updateStressYears(
            f_env, tree.getStressYears()
        );
        tree.setStressYears(new_stress_years);
        
        // 获取物种参数
        const auto& species = species_mgr_->getSpecies(tree.getSpeciesId());
        
        // 计算死亡概率
        double mortality_prob = MortalityCalculator::calcTreeMortality(
            f_env, tree.getStressYears(), species
        );
        
        // 判定死亡
        if (MortalityCalculator::shouldDie(mortality_prob, *rng_)) {
            tree.setAlive(false);
            ++death_count;
        }
    }
    
    return death_count;
}

int MortalityEngine::processSaplingMortality(
    std::vector<Sapling>& saplings,
    const std::vector<GrowthResponse>& growth_responses
) {
    if (!species_mgr_ || !rng_) return 0;
    
    int death_count = 0;
    
    for (size_t i = 0; i < saplings.size(); ++i) {
        auto& sapling = saplings[i];
        if (!sapling.isAlive()) continue;
        
        // 获取生长响应
        double f_env = 1.0;
        if (i < growth_responses.size()) {
            f_env = growth_responses[i].getEnvironmentResponse();
        }
        
        // 更新胁迫年数
        if (f_env < MortalityCalculator::STRESS_THRESHOLD) {
            sapling.incrementStressYears();
        } else {
            sapling.resetStressYears();
        }
        
        // 获取物种参数
        const auto& species = species_mgr_->getSpecies(sapling.getSpeciesId());
        
        // 计算死亡概率（幼树仅使用胁迫死亡）
        double mortality_prob = MortalityCalculator::calcSaplingMortality(
            f_env, sapling.getStressYears(), species
        );
        
        // 判定死亡
        if (MortalityCalculator::shouldDie(mortality_prob, *rng_)) {
            sapling.setAlive(false);
            ++death_count;
        }
    }
    
    return death_count;
}

int MortalityEngine::removeDeadTrees(std::vector<AdultTree>& trees) {
    size_t before = trees.size();
    
    trees.erase(
        std::remove_if(trees.begin(), trees.end(),
            [](const AdultTree& t) { return !t.isAlive(); }),
        trees.end()
    );
    
    return static_cast<int>(before - trees.size());
}

int MortalityEngine::removeDeadSaplings(std::vector<Sapling>& saplings) {
    size_t before = saplings.size();
    
    saplings.erase(
        std::remove_if(saplings.begin(), saplings.end(),
            [](const Sapling& s) { return !s.isAlive(); }),
        saplings.end()
    );
    
    return static_cast<int>(before - saplings.size());
}

} // namespace sfe
