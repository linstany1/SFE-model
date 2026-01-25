/**
 * @file MathUtils.hpp
 * @brief 数学工具函数和常量
 * 
 * 亚高山森林演替模型 (SFE-Model)
 * 
 * 包含模型计算所需的数学工具函数：
 * - 几何计算（距离、面积）
 * - 数值限制（clamp、约束）
 * - 物理公式相关常量
 */

#ifndef SFE_MATH_UTILS_HPP
#define SFE_MATH_UTILS_HPP

#include <cmath>
#include <algorithm>
#include <limits>

namespace sfe {
namespace math {

// ============================================================================
// 数学常量
// ============================================================================

constexpr double PI = 3.1416;
constexpr double TWO_PI = 2.0 * PI;
constexpr double EPSILON = 1e-10;  ///< 浮点比较容差

// ============================================================================
// 基本数学函数
// ============================================================================

/**
 * @brief 计算 x 的平方
 */
template <typename T>
inline constexpr T square(T x) noexcept {
    return x * x;
}

/**
 * @brief 计算 x 的立方
 */
template <typename T>
inline constexpr T cube(T x) noexcept {
    return x * x * x;
}

/**
 * @brief 将值限制在 [min, max] 范围内
 * @param value 输入值
 * @param min_val 最小值
 * @param max_val 最大值
 * @return 限制后的值
 */
template <typename T>
inline constexpr T clamp(T value, T min_val, T max_val) noexcept {
    return std::max(min_val, std::min(value, max_val));
}

/**
 * @brief 将值限制在 [0, 1] 范围内
 */
inline double clamp01(double value) noexcept {
    return clamp(value, 0.0, 1.0);
}

/**
 * @brief 检查两个浮点数是否近似相等
 */
inline bool nearlyEqual(double a, double b, double epsilon = EPSILON) noexcept {
    return std::abs(a - b) < epsilon;
}

/**
 * @brief 检查值是否在有效范围内（非 NaN，非 Inf）
 */
inline bool isFinite(double value) noexcept {
    return std::isfinite(value);
}

/**
 * @brief 安全除法，避免除以零
 * @param numerator 分子
 * @param denominator 分母
 * @param default_value 分母为零时的默认值
 * @return 商或默认值
 */
inline double safeDivide(double numerator, double denominator, 
                         double default_value = 0.0) noexcept {
    if (std::abs(denominator) < EPSILON) {
        return default_value;
    }
    return numerator / denominator;
}

// ============================================================================
// 几何计算
// ============================================================================

/**
 * @brief 计算两点之间的欧氏距离
 * @param x1, y1 第一个点
 * @param x2, y2 第二个点
 * @return 距离
 */
inline double distance(double x1, double y1, double x2, double y2) noexcept {
    double dx = x2 - x1;
    double dy = y2 - y1;
    return std::sqrt(dx * dx + dy * dy);
}

/**
 * @brief 计算两点之间的欧氏距离的平方（避免开方运算）
 */
inline double distanceSquared(double x1, double y1, 
                               double x2, double y2) noexcept {
    double dx = x2 - x1;
    double dy = y2 - y1;
    return dx * dx + dy * dy;
}

/**
 * @brief 计算圆的面积
 * @param radius 半径
 */
inline double circleArea(double radius) noexcept {
    return PI * radius * radius;
}

/**
 * @brief 计算圆周长
 */
inline double circleCircumference(double radius) noexcept {
    return TWO_PI * radius;
}

// ============================================================================
// 模型专用计算
// ============================================================================

/**
 * @brief 计算生长度日校正函数 (Bugmann 1994)
 * @param temp_c 月平均温度 (°C)
 * @return 校正因子
 * 
 * 文档公式 (75):
 * - T <= 5.5: 8.52 * 10^(0.165*T)
 * - 5.5 < T <= 15.5: 187.2 * 10^(0.0908*T)
 * - T > 15.5: -31.8 + 2.377 * T
 */
inline double gddCorrectionFactor(double temp_c) noexcept {
    if (temp_c <= 5.5) {
        return 8.52 * std::pow(10.0, 0.165 * temp_c);
    } else if (temp_c <= 15.5) {
        return 187.2 * std::pow(10.0, 0.0908 * temp_c);
    } else {
        return -31.8 + 2.377 * temp_c;
    }
}

/**
 * @brief 计算单月生长度日 (GDD)
 * @param temp_c 月平均温度 (°C)
 * @param t_base 基准温度，默认 5.5°C
 * @param days_per_month 每月天数，默认 30.5
 * @return 月生长度日
 * 
 * 文档公式 (74)
 */
inline double monthlyGDD(double temp_c, double t_base = 5.5, 
                         double days_per_month = 30.5) noexcept {
    double base_gdd = std::max(temp_c - t_base, 0.0);
    double corr = gddCorrectionFactor(temp_c);
    return base_gdd * days_per_month * corr;
}

/**
 * @brief 温度响应函数 (成树)
 * @param gdd 实际生长度日
 * @param dd_min 最低需求生长度日
 * @return 温度响应因子 [0, 1]
 * 
 * 文档公式 (78): f_temp = max(1 - exp((DD_min - DD) * 0.0013), 0)
 */
inline double temperatureResponseTree(double gdd, double dd_min) noexcept {
    return clamp01(1.0 - std::exp((dd_min - gdd) * 0.0013));
}

/**
 * @brief 温度响应函数 (草本)
 * @param gdd 实际生长度日
 * @param dd_min 最低需求生长度日
 * @return 温度响应因子 [0, 1]
 * 
 * 文档公式 (79): f_temp_herb = max(1 - exp((DD_min - DD) * 0.0027), 0)
 */
inline double temperatureResponseHerb(double gdd, double dd_min) noexcept {
    return clamp01(1.0 - std::exp((dd_min - gdd) * 0.0027));
}

/**
 * @brief 干旱响应函数
 * @param di 干旱指数
 * @param dr_tol 耐旱阈值
 * @return 干旱响应因子 [0, 1]
 * 
 * 文档公式 (73): f_drought = sqrt(max(1 - DI/DrTol, 0))
 */
inline double droughtResponse(double di, double dr_tol) noexcept {
    if (dr_tol <= 0.0) return 0.0;
    return std::sqrt(clamp01(1.0 - di / dr_tol));
}

/**
 * @brief 光照响应函数 - 耐阴型
 * @param al 可用光照 [0, 1]
 * @return 响应因子
 * 
 * 文档公式 (35): f_tol = 1.0106 * (1 - exp(-4.6 * (AL - 0.03)))
 */
inline double lightResponseTolerant(double al) noexcept {
    return 1.0106 * (1.0 - std::exp(-4.6 * (al - 0.03)));
}

/**
 * @brief 光照响应函数 - 不耐阴型
 * @param al 可用光照 [0, 1]
 * @return 响应因子
 * 
 * 文档公式 (36): f_intol = 1.60 * (1 - exp(-1.16 * (AL - 0.15)))
 */
inline double lightResponseIntolerant(double al) noexcept {
    return 1.60 * (1.0 - std::exp(-1.16 * (al - 0.15)));
}

/**
 * @brief 综合光照响应函数
 * @param al 可用光照 [0, 1]
 * @param shade_tolerance 耐阴等级 [1, 5]
 * @return 光响应因子 [0, 1]
 * 
 * 文档公式 (37): f_light = max(0, min(1, f_intol + (ST-1)/4 * (f_tol - f_intol)))
 */
inline double lightResponse(double al, double shade_tolerance) noexcept {
    double f_tol = lightResponseTolerant(al);
    double f_intol = lightResponseIntolerant(al);
    double st_factor = (shade_tolerance - 1.0) / 4.0;
    return clamp01(f_intol + st_factor * (f_tol - f_intol));
}

/**
 * @brief 反解耐阴型光照阈值
 * @param f_prune 保枝阈值 (默认 0.05)
 * @return 对应的光照值
 * 
 * 文档公式 (9-反解函数tol): AL_tol(f) = 0.03 - (1/4.6)*ln(1 - f/1.0106)
 */
inline double inverseLightTolerant(double f_prune) noexcept {
    double arg = 1.0 - f_prune / 1.0106;
    if (arg <= 0.0) return 1.0;  // 边界情况
    return 0.03 - (1.0 / 4.6) * std::log(arg);
}

/**
 * @brief 反解不耐阴型光照阈值
 * @param f_prune 保枝阈值 (默认 0.05)
 * @return 对应的光照值
 * 
 * 文档公式 (9-反解函数intol): AL_intol(f) = 0.15 - (1/1.16)*ln(1 - f/1.6)
 */
inline double inverseLightIntolerant(double f_prune) noexcept {
    double arg = 1.0 - f_prune / 1.6;
    if (arg <= 0.0) return 1.0;  // 边界情况
    return 0.15 - (1.0 / 1.16) * std::log(arg);
}

/**
 * @brief 计算枝条保留所需的最低光照阈值
 * @param shade_tolerance 耐阴等级 [1, 5]
 * @param f_prune 保枝阈值 (默认 0.05)
 * @return 光照阈值
 * 
 * 文档公式 (9)
 */
inline double branchLightThreshold(double shade_tolerance, 
                                    double f_prune = 0.05) noexcept {
    double al_intol = inverseLightIntolerant(f_prune);
    double al_tol = inverseLightTolerant(f_prune);
    double st_factor = (shade_tolerance - 1.0) / 4.0;
    return al_intol + st_factor * (al_tol - al_intol);
}

/**
 * @brief 种子散布核函数（双指数核）
 * @param distance 距离 (m)
 * @param lambda1 短距离特征距离
 * @param lambda2 长距离特征距离
 * @param k_ldd 长距离散布比例
 * @return 核函数值 (seeds/m²)
 * 
 * 文档公式 (1)
 */
inline double dispersalKernel(double distance, double lambda1, 
                               double lambda2, double k_ldd) noexcept {
    if (distance <= 0.0) {
        // 距离为0时返回一个大值，但需要特殊处理
        return 1.0;  // 简化处理
    }
    
    double term1 = (1.0 - k_ldd) / (lambda1 * lambda1) * 
                   std::exp(-distance / lambda1);
    double term2 = k_ldd / (lambda2 * lambda2) * 
                   std::exp(-distance / lambda2);
    
    return (term1 + term2) / (TWO_PI * distance);
}

/**
 * @brief 计算内在死亡率 (基于最大年龄)
 * @param max_age 最大年龄
 * @return 年死亡概率
 * 
 * 文档公式 (110): P_mAge = 1 - 0.01^(1/maxAge)
 * 假设1%的树木能活到最大年龄
 */
inline double intrinsicMortality(int max_age) noexcept {
    if (max_age <= 0) return 1.0;
    return 1.0 - std::pow(0.01, 1.0 / static_cast<double>(max_age));
}

} // namespace math
} // namespace sfe

#endif // SFE_MATH_UTILS_HPP
