/**
 * @file LightEngine.cpp
 * @brief 光照计算引擎实现 (Phase 9 重构版)
 * 
 * 亚高山森林演替模型 (SFE-Model)
 * 
 * Phase 9 实现要点：
 * - 双层空间加速：SpatialHash + CanopyProjectionMap
 * - 幼树分层光照计算
 * - 动态冠底高度公式更新
 */

#include "modules/LightEngine.hpp"
#include "entities/AdultTree.hpp"
#include "entities/Sapling.hpp"
#include <algorithm>
#include <cmath>

namespace sfe {

// ============================================================================
// LightEngine 实现
// ============================================================================

void LightEngine::updateMaxCrownRadius(const std::vector<AdultTree>& trees) {
    max_crown_radius_ = 5.0;  // 默认最小搜索半径
    for (const auto& tree : trees) {
        if (tree.isAlive()) {
            max_crown_radius_ = std::max(max_crown_radius_, tree.getCrownRadius());
        }
    }
}

void LightEngine::buildCanopyProjectionMap(const std::vector<AdultTree>& trees) {
    // 保存树木列表指针
    trees_ptr_ = &trees;
    
    // 清空所有网格的树木ID列表
    for (int y = 0; y < PLOT_SIZE; ++y) {
        for (int x = 0; x < PLOT_SIZE; ++x) {
            canopy_map_.get(x, y).clear();
        }
    }
    
    // 遍历所有成树，将其ID注册到覆盖的网格中
    for (size_t idx = 0; idx < trees.size(); ++idx) {
        const auto& tree = trees[idx];
        if (!tree.isAlive()) continue;
        
        double tx = tree.getX();
        double ty = tree.getY();
        double r_max = tree.getCrownRadius();
        
        // 计算冠层投影范围 [x - R_max, x + R_max] × [y - R_max, y + R_max]
        int x_min = std::max(0, static_cast<int>(std::floor(tx - r_max)));
        int x_max = std::min(PLOT_SIZE - 1, static_cast<int>(std::floor(tx + r_max)));
        int y_min = std::max(0, static_cast<int>(std::floor(ty - r_max)));
        int y_max = std::min(PLOT_SIZE - 1, static_cast<int>(std::floor(ty + r_max)));
        
        // 注册到覆盖的网格
        for (int cy = y_min; cy <= y_max; ++cy) {
            for (int cx = x_min; cx <= x_max; ++cx) {
                // 更精确的检查：网格中心是否在冠幅圆内
                double cell_cx = cx + 0.5;
                double cell_cy = cy + 0.5;
                double dist = std::sqrt((cell_cx - tx) * (cell_cx - tx) + 
                                        (cell_cy - ty) * (cell_cy - ty));
                if (dist <= r_max + 0.71) {  // 0.71 ≈ √0.5，考虑网格对角线
                    canopy_map_.get(cx, cy).push_back(static_cast<int>(idx));
                }
            }
        }
    }
}

void LightEngine::calculateAllLight(std::vector<AdultTree>& trees,
                                     const Grid<double>* herb_lai_grid) {
    // Step 1: 重建空间索引（用于成树互遮阴）
    rebuildSpatialIndex(trees);
    
    // Step 2: 构建冠层投影位图（用于近地层光照加速）
    buildCanopyProjectionMap(trees);
    
    // Step 3: 计算每棵成树的可用光照
    for (auto& tree : trees) {
        if (!tree.isAlive()) continue;
        
        double al = calculateTreeLight(tree);
        tree.setAvailableLight(al);
    }
    
    // Step 4: 计算草本顶部光照网格 (H = 1.0m)
    // 使用 CanopyProjectionMap 加速
    for (int y = 0; y < PLOT_SIZE; ++y) {
        for (int x = 0; x < PLOT_SIZE; ++x) {
            double light = calcGroundLightAtCell(x, y);
            herb_top_light_.set(x, y, light);
        }
    }
    
    // Step 5: 计算真正的地表光照（考虑草本遮挡，H = 0.05m）
    if (herb_lai_grid) {
        for (int y = 0; y < PLOT_SIZE; ++y) {
            for (int x = 0; x < PLOT_SIZE; ++x) {
                double herb_lai = herb_lai_grid->get(x, y);
                double herb_light = herb_top_light_.get(x, y);
                
                // 公式：AL_ground = AL_herb_top × exp(-k_e × LA_herb × (1.0 - 0.05) / 1.0)
                double height_ratio = (params_.herb_height - params_.seedling_height) / 
                                      params_.herb_height;
                double ground = herb_light * std::exp(-params_.extinction_ke * herb_lai * height_ratio);
                ground_light_.set(x, y, math::clamp01(ground));
            }
        }
    } else {
        // 无草本信息时，地表光照等于草本顶部光照
        for (int y = 0; y < PLOT_SIZE; ++y) {
            for (int x = 0; x < PLOT_SIZE; ++x) {
                ground_light_.set(x, y, herb_top_light_.get(x, y));
            }
        }
    }
}

double LightEngine::calculateTreeLight(const AdultTree& tree) const {
    if (!tree.isAlive()) return 0.0;
    
    // 计算有效受光高度 (公式: Z_eff = H - 0.25 × (H - HB))
    double z_eff = tree.calcEffectiveLightHeight();
    
    return calculateLightAtPoint(tree.getX(), tree.getY(), z_eff);
}

double LightEngine::calculateLightAtPoint(double x, double y, double z) const {
    // 使用 SpatialHash 获取候选遮挡树
    auto candidates = spatial_hash_.getShadingCandidates(x, y, max_crown_radius_);
    
    // 累积光学厚度
    double total_optical_depth = 0.0;
    
    for (const AdultTree* shading_tree : candidates) {
        if (!shading_tree || !shading_tree->isAlive()) continue;
        
        // 计算遮阴路径长度
        double path_length = calcShadingPathLength(x, y, z, *shading_tree);
        
        if (path_length > 0.0) {
            // 获取叶面积密度
            double leaf_density = shading_tree->getLeafAreaDensity();
            
            // 累加光学厚度
            total_optical_depth += params_.extinction_ke * leaf_density * path_length;
        }
    }
    
    // Beer-Lambert 定律
    double al = std::exp(-total_optical_depth);
    return math::clamp01(al);
}

double LightEngine::calcSaplingLight(const Sapling& sapling, double herb_lai) const {
    double sap_height = sapling.getHeight();
    double x = sapling.getX();
    double y = sapling.getY();
    
    // 获取网格索引
    int cx = static_cast<int>(std::floor(x));
    int cy = static_cast<int>(std::floor(y));
    cx = std::max(0, std::min(PLOT_SIZE - 1, cx));
    cy = std::max(0, std::min(PLOT_SIZE - 1, cy));
    
    if (sap_height >= params_.herb_height) {
        // 高幼树 (H >= 1.0m)：仅受成树遮阴
        // 直接计算在幼树高度处的光照
        return calculateLightAtPoint(x, y, sap_height);
    } else {
        // 矮幼树 (H < 1.0m)：受成树+草本双重遮阴
        // AL = AL_top_herb × exp(-k_e × LA_herb × (H_herb - H_sap) / H_herb)
        double al_herb_top = herb_top_light_.get(cx, cy);
        
        if (herb_lai <= 0.0) {
            return al_herb_top;
        }
        
        double height_ratio = (params_.herb_height - sap_height) / params_.herb_height;
        double attenuation = std::exp(-params_.extinction_ke * herb_lai * height_ratio);
        
        return math::clamp01(al_herb_top * attenuation);
    }
}

double LightEngine::calcGroundLightAtCell(int cell_x, int cell_y) const {
    // 网格中心坐标
    double x = cell_x + 0.5;
    double y = cell_y + 0.5;
    double z = params_.herb_height;  // 草本顶部高度 (1.0m)
    
    // Phase 9: 使用 CanopyProjectionMap 加速
    // 仅遍历覆盖此网格的树木
    if (!trees_ptr_ || trees_ptr_->empty()) {
        return 1.0;
    }
    
    const auto& tree_indices = canopy_map_.get(cell_x, cell_y);
    
    if (tree_indices.empty()) {
        return 1.0;  // 无遮挡树
    }
    
    // 累积光学厚度
    double total_optical_depth = 0.0;
    
    for (int idx : tree_indices) {
        if (idx < 0 || idx >= static_cast<int>(trees_ptr_->size())) continue;
        
        const AdultTree& shading_tree = (*trees_ptr_)[idx];
        if (!shading_tree.isAlive()) continue;
        
        // 计算遮阴路径长度
        double path_length = calcShadingPathLength(x, y, z, shading_tree);
        
        if (path_length > 0.0) {
            double leaf_density = shading_tree.getLeafAreaDensity();
            total_optical_depth += params_.extinction_ke * leaf_density * path_length;
        }
    }
    
    // Beer-Lambert 定律
    double al = std::exp(-total_optical_depth);
    return math::clamp01(al);
}

double LightEngine::calcCrownEntryHeight(double horizontal_dist,
                                          double crown_base, double tree_height,
                                          double crown_radius, double crown_shape) {
    // 如果水平距离超过冠幅半径，光线不穿过冠层
    if (horizontal_dist >= crown_radius || crown_radius <= 0.0) {
        return crown_base;  // 返回冠底高度表示不遮挡
    }
    
    // 如果水平距离为0（正上方），入射高度为树顶
    if (horizontal_dist <= 0.0) {
        return tree_height;
    }
    
    // 公式: Z_in = HB + (H - HB) × [1 - (d/R_max)^(1/cs)]
    double crown_depth = tree_height - crown_base;
    if (crown_depth <= 0.0) return crown_base;
    
    double ratio = horizontal_dist / crown_radius;
    double exponent = 1.0 / std::max(0.1, crown_shape);  // 避免除以0
    double factor = 1.0 - std::pow(ratio, exponent);
    
    return crown_base + crown_depth * factor;
}

double LightEngine::calcShadingPathLength(double target_x, double target_y, double target_z,
                                           const AdultTree& shading_tree) {
    // 计算水平距离
    double dx = target_x - shading_tree.getX();
    double dy = target_y - shading_tree.getY();
    double horizontal_dist = std::sqrt(dx * dx + dy * dy);
    
    double crown_radius = shading_tree.getCrownRadius();
    double crown_base = shading_tree.getCrownBase();
    double tree_height = shading_tree.getHeight();
    double crown_shape = shading_tree.getCrownShape();
    
    // 条件1: d > R_max，不遮挡
    if (horizontal_dist > crown_radius) {
        return 0.0;
    }
    
    // 计算光线进入冠层的高度
    double z_in = calcCrownEntryHeight(horizontal_dist, crown_base, tree_height,
                                        crown_radius, crown_shape);
    
    // 条件2: 目标点高于冠层入口，不遮挡
    if (target_z >= z_in) {
        return 0.0;
    }
    
    // 公式: L_path = max(0, Z_in - max(z, HB))
    double effective_bottom = std::max(target_z, crown_base);
    double path_length = z_in - effective_bottom;
    
    return std::max(0.0, path_length);
}

double LightEngine::calcLightAttenuation(double leaf_density, double path_length,
                                          double extinction_ke) {
    if (path_length <= 0.0 || leaf_density <= 0.0) {
        return 1.0;
    }
    return std::exp(-extinction_ke * leaf_density * path_length);
}

// ============================================================================
// CrownGeometryCalculator 实现
// ============================================================================

void CrownGeometryCalculator::updateGeometry(AdultTree& tree,
                                              double h_max, double opt_s,
                                              double r_a, double r_b, double cs,
                                              double alpha_base,
                                              double ka1, double ka2, double kc2) {
    // 更新树高（基于当前DBH）
    double new_height = AdultTree::calcHeightFromDbh(h_max, opt_s, tree.getDbh());
    tree.setHeight(new_height);
    
    // 更新冠幅半径 R_max = R_a × ln(DBH) + R_b
    double crown_radius = tree.calcMaxCrownRadius(r_a, r_b);
    tree.setCrownRadius(std::max(0.1, crown_radius));
    
    // 更新冠形系数
    tree.setCrownShape(cs);
    
    // 更新叶面积 LA = KC2 × KA1 × DBH^KA2
    double leaf_area = tree.calcLeafArea(ka1, ka2, kc2);
    tree.setLeafArea(leaf_area);
    
    // 更新开阔地条件下的默认冠底高度
    double hb_open = alpha_base * tree.getHeight();
    if (tree.getCrownBase() < hb_open) {
        tree.setCrownBase(hb_open);
    }
}

void CrownGeometryCalculator::updateCrownBase(AdultTree& tree,
                                               double al_hb,
                                               double al_top,
                                               double l_min,
                                               double alpha_base) {
    double height = tree.getHeight();
    double current_hb = tree.getCrownBase();
    
    // 开阔地基准冠底高度
    double hb_open = alpha_base * height;
    
    // 默认 z* = current_hb（不修枝）
    double z_star = current_hb;
    
    // 修枝判定条件（文档 Ch 2.1.1）：
    // 仅当 AL_HB < L_min < AL_top 时，才需要修枝
    if (al_hb < l_min && l_min < al_top) {
        // 检查数值稳定性
        double ln_diff = std::log(al_top) - std::log(al_hb);
        
        // 若 |ln(AL_top) - ln(AL_HB)| <= 1e-6，不求解 z*
        if (std::abs(ln_diff) > 1e-6 && al_hb > 0.0 && al_top > 0.0) {
            // 公式：z* = HB + (H - HB) × [ln(L_min) - ln(AL_HB)] / [ln(AL_top) - ln(AL_HB)]
            double crown_depth = height - current_hb;
            double ln_l_min = std::log(l_min);
            double ln_al_hb = std::log(al_hb);
            double ln_al_top = std::log(al_top);
            
            z_star = current_hb + crown_depth * (ln_l_min - ln_al_hb) / (ln_al_top - ln_al_hb);
        }
    }
    // 若 AL_HB >= L_min，光照充足，不修枝
    // 若 L_min >= AL_top，即使冠顶光照也不足，由生长/死亡模块处理
    
    // 高度更新规则（文档公式）：
    // HB_new = min(0.9H, max(HB_old, HB_open, z*))
    double new_hb = std::max({current_hb, hb_open, z_star});
    new_hb = std::min(new_hb, 0.9 * height);  // 强制保留至少 10% 冠层
    
    tree.setCrownBase(new_hb);
}

} // namespace sfe
