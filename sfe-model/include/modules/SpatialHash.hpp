/**
 * @file SpatialHash.hpp
 * @brief 空间哈希网格，用于快速空间查询
 * 
 * 亚高山森林演替模型 (SFE-Model)
 * 
 * 功能：
 * - 5×5m 粗网格加速成树邻域查询
 * - 支持快速注册/移除树木
 * - 支持范围查询和邻域查询
 * 
 * 用于光照计算优化：
 * - 快速获取可能遮挡目标点的候选树木
 * - 避免遍历所有树木的 O(N²) 复杂度
 */

#ifndef SFE_SPATIAL_HASH_HPP
#define SFE_SPATIAL_HASH_HPP

#include <vector>
#include <array>
#include <algorithm>
#include <cmath>
#include "core/GlobalConfig.hpp"
#include "core/Types.hpp"

namespace sfe {

// 前向声明
class AdultTree;

/**
 * @class SpatialHash
 * @brief 空间哈希网格
 * 
 * 将 30×30m 样方划分为 6×6 个 5×5m 的粗网格桶，
 * 每个桶存储位于该区域的树木指针列表。
 */
class SpatialHash {
public:
    // 常量
    static constexpr int PLOT_SIZE = GlobalConfig::PLOT_SIZE_M;           ///< 样方尺寸 (30m)
    static constexpr int CELL_SIZE = GlobalConfig::HASH_GRID_SIZE_M;      ///< 哈希格子尺寸 (5m)
    static constexpr int GRID_DIM = PLOT_SIZE / CELL_SIZE;                ///< 网格维度 (6)
    static constexpr int TOTAL_CELLS = GRID_DIM * GRID_DIM;               ///< 总格子数 (36)
    
    /**
     * @brief 默认构造函数
     */
    SpatialHash() {
        clear();
    }
    
    /**
     * @brief 清空所有桶
     */
    void clear() {
        for (auto& bucket : buckets_) {
            bucket.clear();
        }
        tree_count_ = 0;
    }
    
    /**
     * @brief 注册一棵树
     * @param tree 树木指针
     */
    void insert(AdultTree* tree);
    
    /**
     * @brief 批量注册树木
     * @param trees 树木列表
     */
    void insertAll(std::vector<AdultTree>& trees) {
        clear();
        for (auto& tree : trees) {
            insert(&tree);
        }
    }
    
    /**
     * @brief 移除一棵树
     * @param tree 树木指针
     */
    void remove(AdultTree* tree);
    
    /**
     * @brief 更新树木位置（移除后重新插入）
     * @param tree 树木指针
     * @param old_x 旧x坐标
     * @param old_y 旧y坐标
     */
    void update(AdultTree* tree, double old_x, double old_y);
    
    /**
     * @brief 获取指定桶中的所有树木
     * @param gx 网格x索引 [0, GRID_DIM)
     * @param gy 网格y索引 [0, GRID_DIM)
     * @return 树木指针列表
     */
    const std::vector<AdultTree*>& getBucket(int gx, int gy) const {
        return buckets_[bucketIndex(gx, gy)];
    }
    
    /**
     * @brief 获取指定点周围的候选树木（包含邻近桶）
     * @param x 查询点x坐标 (m)
     * @param y 查询点y坐标 (m)
     * @param radius 搜索半径 (m)
     * @return 候选树木列表
     */
    std::vector<AdultTree*> getCandidatesInRadius(double x, double y, double radius) const;
    
    /**
     * @brief 获取可能遮挡指定点的候选树木
     * @param x 目标点x坐标
     * @param y 目标点y坐标
     * @param max_crown_radius 最大冠幅半径（用于确定搜索范围）
     * @return 候选遮挡树木列表
     */
    std::vector<AdultTree*> getShadingCandidates(double x, double y, 
                                                  double max_crown_radius) const {
        return getCandidatesInRadius(x, y, max_crown_radius);
    }
    
    /**
     * @brief 获取指定矩形区域内的所有树木
     * @param x_min 最小x
     * @param y_min 最小y
     * @param x_max 最大x
     * @param y_max 最大y
     * @return 区域内树木列表
     */
    std::vector<AdultTree*> getInRectangle(double x_min, double y_min,
                                            double x_max, double y_max) const;
    
    /**
     * @brief 获取已注册的树木数量
     */
    size_t getTreeCount() const { return tree_count_; }
    
    /**
     * @brief 坐标转换为网格索引
     */
    static int coordToGridX(double x) {
        int gx = static_cast<int>(x / CELL_SIZE);
        return std::max(0, std::min(gx, GRID_DIM - 1));
    }
    
    static int coordToGridY(double y) {
        int gy = static_cast<int>(y / CELL_SIZE);
        return std::max(0, std::min(gy, GRID_DIM - 1));
    }
    
    /**
     * @brief 检查网格索引是否有效
     */
    static bool isValidGrid(int gx, int gy) {
        return gx >= 0 && gx < GRID_DIM && gy >= 0 && gy < GRID_DIM;
    }

private:
    std::array<std::vector<AdultTree*>, TOTAL_CELLS> buckets_;  ///< 哈希桶
    size_t tree_count_ = 0;
    
    /**
     * @brief 计算桶的一维索引
     */
    static int bucketIndex(int gx, int gy) {
        return gy * GRID_DIM + gx;
    }
};

} // namespace sfe

#endif // SFE_SPATIAL_HASH_HPP
