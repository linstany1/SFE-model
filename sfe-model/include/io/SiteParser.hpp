/**
 * @file SiteParser.hpp
 * @brief 站点和种源组成解析器
 * 
 * 亚高山森林演替模型 (SFE-Model)
 * 
 * Task 2.3: 解析 site.txt 和 seed_source_composition.txt
 * 
 * site.txt 格式：
 * plot_id, elevation, latitude, WHC, slope, aspect
 * 
 * seed_source_composition.txt 格式：
 * plot_id, species_id, composition_ratio
 */

#ifndef SFE_SITE_PARSER_HPP
#define SFE_SITE_PARSER_HPP

#include <string>
#include <vector>
#include <map>
#include <stdexcept>
#include "io/FileParser.hpp"
#include "entities/SiteInfo.hpp"
#include "core/Types.hpp"

namespace sfe {

/**
 * @struct PlotConfig
 * @brief 样方配置（对应 site.txt）
 * 
 * Phase 8 更新：新增 plot_name, climate_id, region_group
 */
struct PlotConfig {
    PlotId plot_id = 0;
    std::string plot_name = "";         ///< Phase 8: 用户自定义名称
    int climate_id = 0;                 ///< Phase 8: 气候序列 ID
    std::string region_group = "";      ///< Phase 8: 区域分组标签
    double elevation_m = 0.0;
    double latitude = 0.0;
    double whc_mm = 0.0;                ///< 持水量 (mm)
    double slope_deg = 0.0;
    double aspect_deg = 0.0;
    
    /**
     * @brief 转换为 SiteInfo
     */
    SiteInfo toSiteInfo() const {
        SiteInfo info;
        info.plot_id = plot_id;
        info.plot_name = plot_name;
        info.climate_id = climate_id;
        info.region_group = region_group;
        info.elevation_m = elevation_m;
        info.latitude = latitude;
        info.whc_mm = whc_mm;
        info.slope_deg = slope_deg;
        info.aspect_deg = aspect_deg;
        return info;
    }
};

/**
 * @struct SeedSourceInfo
 * @brief 种源组成信息（对应 seed_source_composition.txt）
 */
struct SeedSourceInfo {
    PlotId plot_id = 0;
    SpeciesId species_id = 0;
    double composition_ratio = 0.0;  ///< 组成比例 [0, 1]
    
    SeedSourceInfo() = default;
    SeedSourceInfo(PlotId pid, SpeciesId sid, double ratio)
        : plot_id(pid), species_id(sid), composition_ratio(ratio) {}
};

/**
 * @class SiteParser
 * @brief 站点文件解析器
 */
class SiteParser {
public:
    /**
     * @brief 解析 site.txt 文件
     * @param filepath 文件路径
     * @return 样方配置列表
     */
    static std::vector<PlotConfig> parse(const std::string& filepath) {
        auto header = FileParser::readHeader(filepath);
        auto table = FileParser::readTable(filepath, true);
        
        // 建立列名映射
        std::map<std::string, int> col_map;
        for (size_t i = 0; i < header.size(); ++i) {
            std::string col_name = header[i];
            std::transform(col_name.begin(), col_name.end(), col_name.begin(), ::tolower);
            col_map[col_name] = static_cast<int>(i);
        }
        
        std::vector<PlotConfig> configs;
        
        for (const auto& row : table) {
            if (row.empty()) continue;
            
            PlotConfig config;
            
            config.plot_id = getInt(row, col_map, "plot_id", 0);
            
            // Phase 8 新增字段
            config.plot_name = getString(row, col_map, "plot_name", "");
            config.climate_id = getInt(row, col_map, "climate_id", 0);
            config.region_group = getString(row, col_map, "region_group", "");
            
            config.elevation_m = getDouble(row, col_map, "elevation", 0.0);
            config.latitude = getDouble(row, col_map, "latitude", 0.0);
            
            // WHC 统一使用 mm
            double whc = getDouble(row, col_map, "whc", 0.0);
            // 假设文件中单位为 mm（如果值很小则可能是 cm）
            if (whc < 50.0 && whc > 0.0) {
                // 可能是 cm，转换为 mm
                whc *= 10.0;
            }
            config.whc_mm = whc;
            
            config.slope_deg = getDouble(row, col_map, "slope", 0.0);
            config.aspect_deg = getDouble(row, col_map, "aspect", 0.0);
            
            configs.push_back(config);
        }
        
        return configs;
    }
    
    /**
     * @brief 解析并返回 SiteInfo 列表
     */
    static std::vector<SiteInfo> parseAsSiteInfo(const std::string& filepath) {
        auto configs = parse(filepath);
        std::vector<SiteInfo> sites;
        sites.reserve(configs.size());
        
        for (const auto& config : configs) {
            sites.push_back(config.toSiteInfo());
        }
        
        return sites;
    }

private:
    static std::string getString(const std::vector<std::string>& row,
                                 const std::map<std::string, int>& col_map,
                                 const std::string& col_name,
                                 const std::string& default_val) {
        auto it = col_map.find(col_name);
        if (it == col_map.end()) return default_val;
        int idx = it->second;
        if (idx < 0 || idx >= static_cast<int>(row.size())) return default_val;
        return row[idx];
    }

    static double getDouble(const std::vector<std::string>& row,
                            const std::map<std::string, int>& col_map,
                            const std::string& col_name,
                            double default_val) {
        auto it = col_map.find(col_name);
        if (it == col_map.end()) return default_val;
        int idx = it->second;
        if (idx < 0 || idx >= static_cast<int>(row.size())) return default_val;
        return FileParser::toDouble(row[idx], default_val);
    }
    
    static int getInt(const std::vector<std::string>& row,
                      const std::map<std::string, int>& col_map,
                      const std::string& col_name,
                      int default_val) {
        auto it = col_map.find(col_name);
        if (it == col_map.end()) return default_val;
        int idx = it->second;
        if (idx < 0 || idx >= static_cast<int>(row.size())) return default_val;
        return FileParser::toInt(row[idx], default_val);
    }
};

/**
 * @class SeedSourceParser
 * @brief 种源组成文件解析器
 */
class SeedSourceParser {
public:
    /**
     * @brief 解析 seed_source_composition.txt 文件
     * @param filepath 文件路径
     * @return 种源信息列表
     */
    static std::vector<SeedSourceInfo> parse(const std::string& filepath) {
        auto header = FileParser::readHeader(filepath);
        auto table = FileParser::readTable(filepath, true);
        
        std::map<std::string, int> col_map;
        for (size_t i = 0; i < header.size(); ++i) {
            std::string col_name = header[i];
            std::transform(col_name.begin(), col_name.end(), col_name.begin(), ::tolower);
            col_map[col_name] = static_cast<int>(i);
        }
        
        std::vector<SeedSourceInfo> sources;
        
        for (const auto& row : table) {
            if (row.size() < 3) continue;
            
            SeedSourceInfo info;
            
            info.plot_id = getInt(row, col_map, "plot_id", 0);
            info.species_id = getInt(row, col_map, "species_id", 0);
            info.composition_ratio = getDouble(row, col_map, "composition_ratio", 0.0);
            
            // 尝试其他可能的列名
            if (col_map.find("composition_ratio") == col_map.end()) {
                info.composition_ratio = getDouble(row, col_map, "ratio", 0.0);
            }
            if (col_map.find("composition_ratio") == col_map.end() &&
                col_map.find("ratio") == col_map.end()) {
                info.composition_ratio = getDouble(row, col_map, "composition", 0.0);
            }
            
            sources.push_back(info);
        }
        
        return sources;
    }
    
    /**
     * @brief 解析并按样方ID分组
     * @param filepath 文件路径
     * @return map<plot_id, vector<SeedSourceInfo>>
     */
    static std::map<PlotId, std::vector<SeedSourceInfo>> parseGrouped(const std::string& filepath) {
        auto sources = parse(filepath);
        std::map<PlotId, std::vector<SeedSourceInfo>> grouped;
        
        for (const auto& src : sources) {
            grouped[src.plot_id].push_back(src);
        }
        
        return grouped;
    }

private:
    static double getDouble(const std::vector<std::string>& row,
                            const std::map<std::string, int>& col_map,
                            const std::string& col_name,
                            double default_val) {
        auto it = col_map.find(col_name);
        if (it == col_map.end()) return default_val;
        int idx = it->second;
        if (idx < 0 || idx >= static_cast<int>(row.size())) return default_val;
        return FileParser::toDouble(row[idx], default_val);
    }
    
    static int getInt(const std::vector<std::string>& row,
                      const std::map<std::string, int>& col_map,
                      const std::string& col_name,
                      int default_val) {
        auto it = col_map.find(col_name);
        if (it == col_map.end()) return default_val;
        int idx = it->second;
        if (idx < 0 || idx >= static_cast<int>(row.size())) return default_val;
        return FileParser::toInt(row[idx], default_val);
    }
};

} // namespace sfe

#endif // SFE_SITE_PARSER_HPP
