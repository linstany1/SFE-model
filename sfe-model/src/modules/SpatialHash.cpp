/**
 * @file SpatialHash.cpp
 * @brief 空间哈希网格实现
 * 
 * 亚高山森林演替模型 (SFE-Model)
 */

#include "modules/SpatialHash.hpp"
#include "entities/AdultTree.hpp"
#include <algorithm>
#include <cmath>

namespace sfe {

void SpatialHash::insert(AdultTree* tree) {
    if (!tree) return;
    
    int gx = coordToGridX(tree->getX());
    int gy = coordToGridY(tree->getY());
    
    auto& bucket = buckets_[bucketIndex(gx, gy)];
    
    // 检查是否已存在
    auto it = std::find(bucket.begin(), bucket.end(), tree);
    if (it == bucket.end()) {
        bucket.push_back(tree);
        ++tree_count_;
    }
}

void SpatialHash::remove(AdultTree* tree) {
    if (!tree) return;
    
    int gx = coordToGridX(tree->getX());
    int gy = coordToGridY(tree->getY());
    
    auto& bucket = buckets_[bucketIndex(gx, gy)];
    auto it = std::find(bucket.begin(), bucket.end(), tree);
    
    if (it != bucket.end()) {
        bucket.erase(it);
        --tree_count_;
    }
}

void SpatialHash::update(AdultTree* tree, double old_x, double old_y) {
    if (!tree) return;
    
    int old_gx = coordToGridX(old_x);
    int old_gy = coordToGridY(old_y);
    int new_gx = coordToGridX(tree->getX());
    int new_gy = coordToGridY(tree->getY());
    
    // 如果网格索引没变，不需要更新
    if (old_gx == new_gx && old_gy == new_gy) {
        return;
    }
    
    // 从旧桶移除
    auto& old_bucket = buckets_[bucketIndex(old_gx, old_gy)];
    auto it = std::find(old_bucket.begin(), old_bucket.end(), tree);
    if (it != old_bucket.end()) {
        old_bucket.erase(it);
    }
    
    // 添加到新桶
    buckets_[bucketIndex(new_gx, new_gy)].push_back(tree);
}

std::vector<AdultTree*> SpatialHash::getCandidatesInRadius(
    double x, double y, double radius
) const {
    std::vector<AdultTree*> candidates;
    
    // 计算需要搜索的网格范围
    int gx_min = coordToGridX(x - radius);
    int gx_max = coordToGridX(x + radius);
    int gy_min = coordToGridY(y - radius);
    int gy_max = coordToGridY(y + radius);
    
    // 遍历所有相关桶
    for (int gy = gy_min; gy <= gy_max; ++gy) {
        for (int gx = gx_min; gx <= gx_max; ++gx) {
            if (!isValidGrid(gx, gy)) continue;
            
            const auto& bucket = buckets_[bucketIndex(gx, gy)];
            for (AdultTree* tree : bucket) {
                if (!tree || !tree->isAlive()) continue;
                
                // 粗略距离检查（精确检查由调用者完成）
                double dx = tree->getX() - x;
                double dy = tree->getY() - y;
                double dist_sq = dx * dx + dy * dy;
                
                // 使用冠幅半径扩展搜索范围
                double search_radius = radius + tree->getCrownRadius();
                if (dist_sq <= search_radius * search_radius) {
                    candidates.push_back(tree);
                }
            }
        }
    }
    
    return candidates;
}

std::vector<AdultTree*> SpatialHash::getInRectangle(
    double x_min, double y_min,
    double x_max, double y_max
) const {
    std::vector<AdultTree*> result;
    
    int gx_min = coordToGridX(x_min);
    int gx_max = coordToGridX(x_max);
    int gy_min = coordToGridY(y_min);
    int gy_max = coordToGridY(y_max);
    
    for (int gy = gy_min; gy <= gy_max; ++gy) {
        for (int gx = gx_min; gx <= gx_max; ++gx) {
            if (!isValidGrid(gx, gy)) continue;
            
            const auto& bucket = buckets_[bucketIndex(gx, gy)];
            for (AdultTree* tree : bucket) {
                if (!tree || !tree->isAlive()) continue;
                
                double tx = tree->getX();
                double ty = tree->getY();
                
                if (tx >= x_min && tx <= x_max && ty >= y_min && ty <= y_max) {
                    result.push_back(tree);
                }
            }
        }
    }
    
    return result;
}

} // namespace sfe
