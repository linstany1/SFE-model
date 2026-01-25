/**
 * @file Types.hpp
 * @brief 项目通用类型定义
 * 
 * 亚高山森林演替模型 (SFE-Model)
 * 
 * 定义项目中使用的基本类型和类型别名
 */

#ifndef SFE_TYPES_HPP
#define SFE_TYPES_HPP

#include <cstdint>
#include <vector>
#include <array>
#include <memory>
#include <string>
#include <map>

namespace sfe {

// ============================================================================
// 基本类型别名
// ============================================================================

using TreeId = uint32_t;           ///< 树木唯一标识符
using SpeciesId = int;             ///< 物种标识符
using PlotId = int;                ///< 样方标识符
using Year = int;                  ///< 年份
using Month = int;                 ///< 月份 [1, 12]

// ============================================================================
// 坐标类型
// ============================================================================

/**
 * @struct Point2D
 * @brief 二维点坐标（浮点）
 */
struct Point2D {
    double x = 0.0;
    double y = 0.0;
    
    Point2D() = default;
    Point2D(double x_, double y_) : x(x_), y(y_) {}
    
    bool operator==(const Point2D& other) const {
        return x == other.x && y == other.y;
    }
    
    bool operator!=(const Point2D& other) const {
        return !(*this == other);
    }
};

/**
 * @struct CellIndex
 * @brief 网格单元索引（整数）
 */
struct CellIndex {
    int x = 0;
    int y = 0;
    
    CellIndex() = default;
    CellIndex(int x_, int y_) : x(x_), y(y_) {}
    
    bool operator==(const CellIndex& other) const {
        return x == other.x && y == other.y;
    }
    
    bool operator!=(const CellIndex& other) const {
        return !(*this == other);
    }
    
    /**
     * @brief 从浮点坐标转换
     */
    static CellIndex fromPoint(const Point2D& p) {
        return CellIndex(static_cast<int>(std::floor(p.x)), 
                         static_cast<int>(std::floor(p.y)));
    }
};

// ============================================================================
// 月度数据类型
// ============================================================================

/**
 * @brief 12个月的数据数组
 */
template <typename T>
using MonthlyArray = std::array<T, 12>;

/**
 * @brief 月度气候数据
 */
using MonthlyDouble = MonthlyArray<double>;

// ============================================================================
// 容器类型别名
// ============================================================================

template <typename T>
using Vec = std::vector<T>;

template <typename K, typename V>
using Map = std::map<K, V>;

template <typename T>
using UniquePtr = std::unique_ptr<T>;

template <typename T>
using SharedPtr = std::shared_ptr<T>;

// ============================================================================
// ID 生成器
// ============================================================================

/**
 * @class IdGenerator
 * @brief 线程安全的唯一ID生成器
 */
class IdGenerator {
public:
    static IdGenerator& getInstance() {
        static IdGenerator instance;
        return instance;
    }
    
    /**
     * @brief 生成新的树木ID
     */
    TreeId nextTreeId() {
        return ++current_tree_id_;
    }
    
    /**
     * @brief 重置ID计数器（通常在新模拟开始时调用）
     */
    void reset() {
        current_tree_id_ = 0;
    }
    
    /**
     * @brief 获取当前最大ID
     */
    TreeId currentMaxId() const {
        return current_tree_id_;
    }

private:
    IdGenerator() : current_tree_id_(0) {}
    TreeId current_tree_id_;
};

/**
 * @brief 便捷函数：生成新的树木ID
 */
inline TreeId generateTreeId() {
    return IdGenerator::getInstance().nextTreeId();
}

// ============================================================================
// 枚举类型
// ============================================================================

/**
 * @enum TreeStage
 * @brief 树木生命阶段
 */
enum class TreeStage {
    Seedling,   ///< 幼苗 (刚定居)
    Sapling,    ///< 幼树 (H < 1.37m)
    Adult       ///< 成树 (H >= 1.37m)
};

/**
 * @enum PlotType
 * @brief 样方类型
 */
enum class PlotType {
    Target,     ///< 目标样方（参与动态模拟）
    Source      ///< 种源样方（仅作为种子源）
};

/**
 * @enum DisturbanceType
 * @brief 干扰类型
 */
enum class DisturbanceType {
    None,       ///< 无干扰
    Fire        ///< 火灾
};

/**
 * @brief 将 TreeStage 转换为字符串
 */
inline const char* toString(TreeStage stage) {
    switch (stage) {
        case TreeStage::Seedling: return "Seedling";
        case TreeStage::Sapling:  return "Sapling";
        case TreeStage::Adult:    return "Adult";
        default: return "Unknown";
    }
}

/**
 * @brief 将 PlotType 转换为字符串
 */
inline const char* toString(PlotType type) {
    switch (type) {
        case PlotType::Target: return "Target";
        case PlotType::Source: return "Source";
        default: return "Unknown";
    }
}

/**
 * @brief 将 DisturbanceType 转换为字符串
 */
inline const char* toString(DisturbanceType type) {
    switch (type) {
        case DisturbanceType::None: return "None";
        case DisturbanceType::Fire: return "Fire";
        default: return "Unknown";
    }
}

} // namespace sfe

#endif // SFE_TYPES_HPP
