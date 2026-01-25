/**
 * @file Growth.cpp
 * @brief 生长模块实现
 * 
 * 亚高山森林演替模型 (SFE-Model)
 */

#include "modules/Growth.hpp"
#include "entities/AdultTree.hpp"
#include "entities/Sapling.hpp"
#include "entities/HerbLayer.hpp"
#include <algorithm>
#include <cmath>

namespace sfe {

// ============================================================================
// HerbGrowthCalculator 实现
// ============================================================================

double HerbGrowthCalculator::calcGrowth(
    double current_biomass,
    double available_light,
    double drought_index,
    double gdd,
    const HerbGrowthParams& params
) {
    // 计算环境响应
    double f_light = EnvironmentResponse::calcLightResponse(available_light, params.shade_tol);
    double f_temp = EnvironmentResponse::calcTemperatureResponseHerb(gdd, params.dd_min);
    double f_drought = EnvironmentResponse::calcDroughtResponse(drought_index, params.dr_tol);
    
    double f_env = f_light * f_temp * f_drought;
    
    // 计算有效参数
    double K = calcEffectiveK(params.k_max, f_env);
    double r = calcEffectiveR(params.r_max, f_env);
    
    // 防止 K 过小导致除法问题
    if (K < MIN_BIOMASS) {
        K = MIN_BIOMASS;
    }
    
    // Logistic 生长 (公式101)
    // B_{Y+1} = B_Y + r × (1 - B/K) × B
    double growth_term = r * (1.0 - current_biomass / K) * current_biomass;
    double new_biomass = current_biomass + growth_term;
    
    // 应用最小生物量约束 (公式102)
    return std::max(new_biomass, MIN_BIOMASS);
}

double HerbGrowthCalculator::calcGrowth(
    double current_biomass,
    double available_light,
    double drought_index,
    double gdd,
    const SpeciesProfile& herb_species
) {
    HerbGrowthParams params;
    params.k_max = herb_species.k_max_herb;
    params.r_max = herb_species.r_max_herb;
    params.dd_min = herb_species.dd_min;
    params.dr_tol = herb_species.dr_tol;
    params.shade_tol = herb_species.shade_tol;
    
    return calcGrowth(current_biomass, available_light, drought_index, gdd, params);
}

// ============================================================================
// SaplingGrowthCalculator 实现
// ============================================================================

std::pair<GrowthResponse, double> SaplingGrowthCalculator::calcGrowth(
    const Sapling& sapling,
    double available_light,
    double drought_index,
    double gdd,
    const SpeciesProfile& species
) {
    GrowthResponse response;
    
    // 计算环境响应
    response.f_light = EnvironmentResponse::calcLightResponse(
        available_light, species.shade_tol
    );
    response.f_temperature = EnvironmentResponse::calcTemperatureResponse(
        gdd, species.dd_min
    );
    
    // 幼树使用放大的干旱指数 (公式88)
    double sapling_di = calcSaplingDroughtIndex(drought_index);
    response.f_drought = EnvironmentResponse::calcDroughtResponse(
        sapling_di, species.dr_tol
    );
    
    double f_env = response.getEnvironmentResponse();
    
    // 计算潜在高度 (公式86/106)
    double current_height = sapling.getHeight();
    double potential_height = calcPotentialHeight(
        current_height, species.h_max_m, species.gs
    );
    
    // 潜在增量
    double delta_h_opt = potential_height - current_height;
    
    // 实际增量 (公式87)
    double delta_h = delta_h_opt * f_env;
    
    // 新高度
    double new_height = current_height + delta_h;
    new_height = std::max(new_height, INITIAL_HEIGHT);  // 确保不低于初始高度
    
    return {response, new_height};
}

// ============================================================================
// TreeGrowthCalculator 实现
// ============================================================================

std::pair<GrowthResponse, double> TreeGrowthCalculator::calcGrowth(
    const AdultTree& tree,
    double available_light,
    double drought_index,
    double gdd,
    const SpeciesProfile& species
) {
    GrowthResponse response;
    
    // 计算环境响应
    response.f_light = EnvironmentResponse::calcLightResponse(
        available_light, species.shade_tol
    );
    response.f_temperature = EnvironmentResponse::calcTemperatureResponse(
        gdd, species.dd_min
    );
    response.f_drought = EnvironmentResponse::calcDroughtResponse(
        drought_index, species.dr_tol
    );
    
    double f_env = response.getEnvironmentResponse();
    
    // 计算最优胸径增量 (公式105)
    double dbh_opt = calcOptimalDbhIncrement(
        tree.getDbh(), tree.getHeight(),
        species.dbh_max_cm, species.h_max_m,
        species.ga, species.opt_s
    );
    
    // 实际胸径增量 (公式107)
    double delta_dbh = calcActualDbhIncrement(dbh_opt, f_env);
    
    return {response, delta_dbh};
}

// ============================================================================
// GrowthEngine 实现
// ============================================================================

void GrowthEngine::processHerbGrowth(
    HerbLayer& herb_layer,
    const Grid<double>& ground_light,
    double drought_index,
    double gdd
) {
    if (!species_mgr_) return;
    
    // 获取草本物种参数 (species_id = 0)
    const auto& herb_species = species_mgr_->getHerbSpecies();
    
    // 使用常量网格大小
    constexpr int GRID_SIZE = HerbLayer::GRID_SIZE;
    
    for (int y = 0; y < GRID_SIZE; ++y) {
        for (int x = 0; x < GRID_SIZE; ++x) {
            double current_biomass = herb_layer.getBiomass(x, y);
            double light = ground_light.get(x, y);
            
            double new_biomass = HerbGrowthCalculator::calcGrowth(
                current_biomass, light, drought_index, gdd, herb_species
            );
            
            herb_layer.setBiomass(x, y, new_biomass);
        }
    }
}

std::vector<GrowthResponse> GrowthEngine::processSaplingGrowth(
    std::vector<Sapling>& saplings,
    const Grid<double>& ground_light,
    double drought_index,
    double gdd
) {
    std::vector<GrowthResponse> responses;
    responses.reserve(saplings.size());
    
    if (!species_mgr_) {
        return responses;
    }
    
    for (auto& sapling : saplings) {
        if (!sapling.isAlive()) {
            responses.push_back(GrowthResponse{});
            continue;
        }
        
        const auto& species = species_mgr_->getSpecies(sapling.getSpeciesId());
        
        // 获取幼树位置的光照
        int cell_x = static_cast<int>(sapling.getX());
        int cell_y = static_cast<int>(sapling.getY());
        cell_x = std::max(0, std::min(cell_x, ground_light.width() - 1));
        cell_y = std::max(0, std::min(cell_y, ground_light.height() - 1));
        double light = ground_light.get(cell_x, cell_y);
        
        // 计算生长
        auto [response, new_height] = SaplingGrowthCalculator::calcGrowth(
            sapling, light, drought_index, gdd, species
        );
        
        // 更新幼树
        sapling.setHeight(new_height);
        sapling.incrementAge();
        
        responses.push_back(response);
    }
    
    return responses;
}

std::vector<GrowthResponse> GrowthEngine::processTreeGrowth(
    std::vector<AdultTree>& trees,
    double drought_index,
    double gdd
) {
    std::vector<GrowthResponse> responses;
    responses.reserve(trees.size());
    
    if (!species_mgr_) {
        return responses;
    }
    
    for (auto& tree : trees) {
        if (!tree.isAlive()) {
            responses.push_back(GrowthResponse{});
            continue;
        }
        
        const auto& species = species_mgr_->getSpecies(tree.getSpeciesId());
        
        // 获取树木可用光照
        double available_light = tree.getAvailableLight();
        
        // 计算生长
        auto [response, delta_dbh] = TreeGrowthCalculator::calcGrowth(
            tree, available_light, drought_index, gdd, species
        );
        
        // 更新胸径
        double new_dbh = tree.getDbh() + delta_dbh;
        tree.setDbh(new_dbh);
        
        // 更新高度 (公式106)
        double new_height = TreeGrowthCalculator::calcHeightFromDbh(
            new_dbh, species.h_max_m, species.opt_s
        );
        tree.setHeight(new_height);
        
        // 更新叶面积 (公式108)
        double new_la = TreeGrowthCalculator::calcLeafArea(
            new_dbh, species.kc2, species.ka1, species.ka2
        );
        tree.setLeafArea(new_la);
        
        // 更新冠幅
        double new_cr = TreeGrowthCalculator::calcCrownRadius(
            new_dbh, species.r_a, species.r_b
        );
        tree.setCrownRadius(new_cr);
        
        // 增加年龄
        tree.incrementAge();
        
        responses.push_back(response);
    }
    
    return responses;
}

int GrowthEngine::processRecruitment(
    std::vector<Sapling>& saplings,
    std::vector<AdultTree>& trees,
    TreeId& next_tree_id
) {
    if (!species_mgr_) return 0;
    
    int recruited_count = 0;
    
    // 找出可以晋升的幼树
    std::vector<size_t> to_remove;
    
    for (size_t i = 0; i < saplings.size(); ++i) {
        auto& sapling = saplings[i];
        
        if (!sapling.isAlive()) continue;
        
        if (SaplingGrowthCalculator::canRecruit(sapling.getHeight())) {
            // 获取物种参数
            const auto& species = species_mgr_->getSpecies(sapling.getSpeciesId());
            
            // 计算初始胸径（反解高度-胸径关系）
            double init_dbh = 1.0;  // 默认初始胸径 1 cm
            if (species.h_max_m > 1.37 && species.opt_s > 0) {
                double h_range = species.h_max_m - 1.37;
                double height_ratio = (sapling.getHeight() - 1.37) / h_range;
                height_ratio = std::max(0.01, std::min(0.99, height_ratio));
                init_dbh = -h_range / species.opt_s * std::log(1.0 - height_ratio);
                init_dbh = std::max(0.5, std::min(init_dbh, 5.0));
            }
            
            // 创建新成树
            AdultTree new_tree(
                next_tree_id++,
                sapling.getSpeciesId(),
                sapling.getX(),
                sapling.getY(),
                init_dbh,
                sapling.getHeight(),
                sapling.getAge()
            );
            
            // 计算初始属性
            new_tree.setLeafArea(TreeGrowthCalculator::calcLeafArea(
                init_dbh, species.kc2, species.ka1, species.ka2
            ));
            new_tree.setCrownRadius(TreeGrowthCalculator::calcCrownRadius(
                init_dbh, species.r_a, species.r_b
            ));
            
            trees.push_back(std::move(new_tree));
            to_remove.push_back(i);
            ++recruited_count;
        }
    }
    
    // 从后往前移除已晋升的幼树
    std::sort(to_remove.rbegin(), to_remove.rend());
    for (size_t idx : to_remove) {
        saplings.erase(saplings.begin() + static_cast<std::ptrdiff_t>(idx));
    }
    
    return recruited_count;
}

} // namespace sfe
