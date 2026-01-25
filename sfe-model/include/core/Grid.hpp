/**
 * @file Grid.hpp
 * @brief 通用二维网格模板类，用于空间数据管理
 * 
 * 亚高山森林演替模型 (SFE-Model)
 * 
 * 用途：
 * - LIF（地表光照强度）网格
 * - HerbBiomass（草本生物量）网格
 * - SeedDensity（种子密度）网格
 * - 任何需要空间显式表示的状态变量
 * 
 * 特性：
 * - 模板化支持任意数据类型
 * - 边界检查（可配置）
 * - 高效的连续内存布局
 * - 支持迭代器访问
 */

#ifndef SFE_GRID_HPP
#define SFE_GRID_HPP

#include <vector>
#include <stdexcept>
#include <algorithm>
#include <functional>
#include <cmath>
#include <cassert>

namespace sfe {

/**
 * @class Grid
 * @brief 二维网格模板类
 * @tparam T 存储的数据类型
 * 
 * 坐标系统：
 * - 原点 (0, 0) 在左下角
 * - x 轴向右增加
 * - y 轴向上增加
 * - 内部使用行优先存储 (row-major)
 */
template <typename T>
class Grid {
public:
    using value_type = T;
    using iterator = typename std::vector<T>::iterator;
    using const_iterator = typename std::vector<T>::const_iterator;
    
    // ========================================================================
    // 构造函数
    // ========================================================================
    
    /**
     * @brief 默认构造函数，创建空网格
     */
    Grid() : width_(0), height_(0), data_() {}
    
    /**
     * @brief 创建指定大小的网格
     * @param width 宽度（x方向）
     * @param height 高度（y方向）
     */
    Grid(int width, int height) 
        : width_(width), height_(height), data_(width * height) {
        if (width < 0 || height < 0) {
            throw std::invalid_argument("Grid dimensions must be non-negative");
        }
    }
    
    /**
     * @brief 创建指定大小的网格并用初始值填充
     * @param width 宽度（x方向）
     * @param height 高度（y方向）
     * @param initial_value 初始值
     */
    Grid(int width, int height, const T& initial_value) 
        : width_(width), height_(height), data_(width * height, initial_value) {
        if (width < 0 || height < 0) {
            throw std::invalid_argument("Grid dimensions must be non-negative");
        }
    }
    
    // ========================================================================
    // 基本属性
    // ========================================================================
    
    /**
     * @brief 获取网格宽度
     */
    int width() const noexcept { return width_; }
    
    /**
     * @brief 获取网格高度
     */
    int height() const noexcept { return height_; }
    
    /**
     * @brief 获取网格总单元数
     */
    int size() const noexcept { return width_ * height_; }
    
    /**
     * @brief 检查网格是否为空
     */
    bool empty() const noexcept { return data_.empty(); }
    
    // ========================================================================
    // 元素访问（带边界检查）
    // ========================================================================
    
    /**
     * @brief 获取指定位置的值（带边界检查）
     * @param x x坐标 [0, width)
     * @param y y坐标 [0, height)
     * @return 对应位置的值引用
     * @throws std::out_of_range 坐标越界时
     */
    T& get(int x, int y) {
        checkBounds(x, y);
        return data_[index(x, y)];
    }
    
    /**
     * @brief 获取指定位置的值（const版本，带边界检查）
     */
    const T& get(int x, int y) const {
        checkBounds(x, y);
        return data_[index(x, y)];
    }
    
    /**
     * @brief 设置指定位置的值（带边界检查）
     * @param x x坐标
     * @param y y坐标
     * @param value 要设置的值
     */
    void set(int x, int y, const T& value) {
        checkBounds(x, y);
        data_[index(x, y)] = value;
    }
    
    // ========================================================================
    // 快速访问（无边界检查，性能关键路径使用）
    // ========================================================================
    
    /**
     * @brief 快速获取值（无边界检查）
     * @warning 调用者必须确保坐标有效
     */
    T& getUnchecked(int x, int y) noexcept {
        assert(isValid(x, y));
        return data_[index(x, y)];
    }
    
    const T& getUnchecked(int x, int y) const noexcept {
        assert(isValid(x, y));
        return data_[index(x, y)];
    }
    
    /**
     * @brief 快速设置值（无边界检查）
     */
    void setUnchecked(int x, int y, const T& value) noexcept {
        assert(isValid(x, y));
        data_[index(x, y)] = value;
    }
    
    // ========================================================================
    // 运算符重载
    // ========================================================================
    
    /**
     * @brief 括号运算符访问（带边界检查）
     */
    T& operator()(int x, int y) {
        return get(x, y);
    }
    
    const T& operator()(int x, int y) const {
        return get(x, y);
    }
    
    // ========================================================================
    // 边界检查
    // ========================================================================
    
    /**
     * @brief 检查坐标是否在有效范围内
     * @param x x坐标
     * @param y y坐标
     * @return true 如果坐标有效
     */
    bool isValid(int x, int y) const noexcept {
        return x >= 0 && x < width_ && y >= 0 && y < height_;
    }
    
    /**
     * @brief 检查浮点坐标对应的网格是否有效
     * @param x 浮点x坐标
     * @param y 浮点y坐标
     */
    bool isValidFloat(double x, double y) const noexcept {
        return isValid(static_cast<int>(std::floor(x)), 
                       static_cast<int>(std::floor(y)));
    }
    
    /**
     * @brief 将坐标限制在有效范围内
     * @param x x坐标（会被修改）
     * @param y y坐标（会被修改）
     */
    void clamp(int& x, int& y) const noexcept {
        x = std::max(0, std::min(x, width_ - 1));
        y = std::max(0, std::min(y, height_ - 1));
    }
    
    // ========================================================================
    // 批量操作
    // ========================================================================
    
    /**
     * @brief 用指定值填充整个网格
     * @param value 填充值
     */
    void fill(const T& value) {
        std::fill(data_.begin(), data_.end(), value);
    }
    
    /**
     * @brief 将所有值设为零/默认值
     */
    void clear() {
        fill(T{});
    }
    
    /**
     * @brief 重置网格大小
     * @param new_width 新宽度
     * @param new_height 新高度
     * @param initial_value 初始值
     */
    void resize(int new_width, int new_height, const T& initial_value = T{}) {
        if (new_width < 0 || new_height < 0) {
            throw std::invalid_argument("Grid dimensions must be non-negative");
        }
        width_ = new_width;
        height_ = new_height;
        data_.assign(new_width * new_height, initial_value);
    }
    
    // ========================================================================
    // 数值操作（仅对算术类型有效）
    // ========================================================================
    
    /**
     * @brief 计算所有元素的和
     * @return 元素总和
     */
    template <typename U = T>
    typename std::enable_if<std::is_arithmetic<U>::value, U>::type
    sum() const {
        U total = U{};
        for (const auto& val : data_) {
            total += val;
        }
        return total;
    }
    
    /**
     * @brief 计算所有元素的平均值
     * @return 平均值
     */
    template <typename U = T>
    typename std::enable_if<std::is_arithmetic<U>::value, double>::type
    mean() const {
        if (data_.empty()) return 0.0;
        return static_cast<double>(sum()) / static_cast<double>(data_.size());
    }
    
    /**
     * @brief 查找最大值
     */
    template <typename U = T>
    typename std::enable_if<std::is_arithmetic<U>::value, U>::type
    max() const {
        if (data_.empty()) return U{};
        return *std::max_element(data_.begin(), data_.end());
    }
    
    /**
     * @brief 查找最小值
     */
    template <typename U = T>
    typename std::enable_if<std::is_arithmetic<U>::value, U>::type
    min() const {
        if (data_.empty()) return U{};
        return *std::min_element(data_.begin(), data_.end());
    }
    
    /**
     * @brief 将所有值乘以标量
     * @param scalar 乘数
     */
    template <typename U = T>
    typename std::enable_if<std::is_arithmetic<U>::value, void>::type
    multiply(U scalar) {
        for (auto& val : data_) {
            val *= scalar;
        }
    }
    
    /**
     * @brief 将所有值加上常数
     * @param offset 加数
     */
    template <typename U = T>
    typename std::enable_if<std::is_arithmetic<U>::value, void>::type
    add(U offset) {
        for (auto& val : data_) {
            val += offset;
        }
    }
    
    /**
     * @brief 将所有值限制在 [min_val, max_val] 范围内
     */
    template <typename U = T>
    typename std::enable_if<std::is_arithmetic<U>::value, void>::type
    clampValues(U min_val, U max_val) {
        for (auto& val : data_) {
            val = std::max(min_val, std::min(val, max_val));
        }
    }
    
    // ========================================================================
    // 遍历与函数式操作
    // ========================================================================
    
    /**
     * @brief 对每个单元格应用函数
     * @param func 函数 func(x, y, value&)
     */
    template <typename Func>
    void forEach(Func&& func) {
        for (int y = 0; y < height_; ++y) {
            for (int x = 0; x < width_; ++x) {
                func(x, y, data_[index(x, y)]);
            }
        }
    }
    
    /**
     * @brief 对每个单元格应用函数（const版本）
     * @param func 函数 func(x, y, const value&)
     */
    template <typename Func>
    void forEach(Func&& func) const {
        for (int y = 0; y < height_; ++y) {
            for (int x = 0; x < width_; ++x) {
                func(x, y, data_[index(x, y)]);
            }
        }
    }
    
    /**
     * @brief 应用变换函数到每个元素
     * @param func 变换函数 T -> T
     */
    template <typename Func>
    void transform(Func&& func) {
        for (auto& val : data_) {
            val = func(val);
        }
    }
    
    // ========================================================================
    // 迭代器支持
    // ========================================================================
    
    iterator begin() noexcept { return data_.begin(); }
    iterator end() noexcept { return data_.end(); }
    const_iterator begin() const noexcept { return data_.begin(); }
    const_iterator end() const noexcept { return data_.end(); }
    const_iterator cbegin() const noexcept { return data_.cbegin(); }
    const_iterator cend() const noexcept { return data_.cend(); }
    
    // ========================================================================
    // 原始数据访问
    // ========================================================================
    
    /**
     * @brief 获取底层数据指针
     */
    T* data() noexcept { return data_.data(); }
    const T* data() const noexcept { return data_.data(); }
    
    /**
     * @brief 获取底层 vector 引用
     */
    std::vector<T>& rawData() noexcept { return data_; }
    const std::vector<T>& rawData() const noexcept { return data_; }

private:
    int width_;              ///< 网格宽度
    int height_;             ///< 网格高度
    std::vector<T> data_;    ///< 数据存储（行优先）
    
    /**
     * @brief 计算一维数组索引（行优先）
     */
    int index(int x, int y) const noexcept {
        return y * width_ + x;
    }
    
    /**
     * @brief 边界检查，越界时抛出异常
     */
    void checkBounds(int x, int y) const {
        if (!isValid(x, y)) {
            throw std::out_of_range(
                "Grid index out of range: (" + std::to_string(x) + ", " + 
                std::to_string(y) + ") not in [0," + std::to_string(width_) + 
                ") x [0," + std::to_string(height_) + ")"
            );
        }
    }
};

// ============================================================================
// 类型别名（常用特化）
// ============================================================================

/// 浮点数网格（用于光照、生物量等）
using GridDouble = Grid<double>;

/// 整数网格（用于计数等）
using GridInt = Grid<int>;

/// 布尔网格（用于掩码等）
using GridBool = Grid<bool>;

} // namespace sfe

#endif // SFE_GRID_HPP
