/**
 * @file SiteInfo.hpp
 * @brief 样地基础信息结构体
 * 
 * 亚高山森林演替模型 (SFE-Model)
 * Phase 8: 支持多气候ID映射和区域分组
 * 
 * 对应 site.txt 文件中的字段
 */

#ifndef SFE_SITE_INFO_HPP
#define SFE_SITE_INFO_HPP

#include <string>
#include "core/Types.hpp"

namespace sfe {

/**
 * @struct SiteInfo
 * @brief 样地基础配置信息
 * 
 * Phase 8 更新字段映射 site.txt:
 * - plot_id: int (系统内部唯一索引，从0开始连续编号)
 * - plot_name: string (用户自定义名称，如 "sw1")
 * - climate_id: int (关联的气候序列 ID)
 * - region_group: string (区域分组标签)
 * - elevation: double (海拔, m)
 * - latitude: double (纬度)
 * - WHC: double (土壤持水量, mm)
 * - slope: double (坡度, 度)
 * - aspect: double (坡向, 度)
 */
struct SiteInfo {
    // === 基础标识 (Phase 8 新增) ===
    PlotId plot_id = 0;               ///< 系统内部唯一索引（从0开始连续编号）
    std::string plot_name = "";       ///< 用户自定义名称（如 "sw1", "sw2"）
    int climate_id = 0;               ///< 关联的气候序列 ID（关键：用于多气候支持）
    std::string region_group = "";    ///< 区域分组标签（如 "Gongga_South"）
    
    // === 地形属性 ===
    double elevation_m = 0.0;         ///< 海拔 (m)
    double latitude = 0.0;            ///< 纬度 (度)
    double slope_deg = 0.0;           ///< 坡度 (度)
    double aspect_deg = 0.0;          ///< 坡向 (度, 0=北, 90=东)
    
    // === 土壤属性 ===
    double whc_mm = 0.0;              ///< 土壤持水量 (mm) - 统一使用 mm
    
    // ========================================================================
    // 构造函数
    // ========================================================================
    
    SiteInfo() = default;
    
    /**
     * @brief 完整构造函数 (Phase 8)
     */
    SiteInfo(PlotId id, const std::string& name, int clim_id, 
             const std::string& region, double elevation, double lat,
             double whc, double slope, double aspect)
        : plot_id(id)
        , plot_name(name)
        , climate_id(clim_id)
        , region_group(region)
        , elevation_m(elevation)
        , latitude(lat)
        , slope_deg(slope)
        , aspect_deg(aspect)
        , whc_mm(whc)
    {
    }
    
    /**
     * @brief 兼容旧版构造函数
     */
    SiteInfo(PlotId id, double elevation, double lat, 
             double whc, double slope, double aspect)
        : plot_id(id)
        , climate_id(0)
        , elevation_m(elevation)
        , latitude(lat)
        , slope_deg(slope)
        , aspect_deg(aspect)
        , whc_mm(whc)
    {
    }
    
    // ========================================================================
    // 转换方法
    // ========================================================================
    
    /**
     * @brief 获取持水量（兼容旧代码，单位：mm）
     */
    double getWhcMm() const {
        return whc_mm;
    }
    
    /**
     * @brief 获取坡度（弧度）
     */
    double getSlopeRad() const {
        return slope_deg * 3.1416 / 180.0;
    }
    
    /**
     * @brief 获取坡向（弧度）
     */
    double getAspectRad() const {
        return aspect_deg * 3.14159265358979 / 180.0;
    }
    
    /**
     * @brief 获取完整描述字符串（用于日志）
     */
    std::string toString() const {
        return plot_name.empty() 
            ? "Plot_" + std::to_string(plot_id)
            : plot_name + " (ID=" + std::to_string(plot_id) + ")";
    }
    
    /**
     * @brief 检查是否有效
     */
    bool isValid() const {
        return whc_mm > 0 && elevation_m > 0;
    }
};

} // namespace sfe

#endif // SFE_SITE_INFO_HPP
