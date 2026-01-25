/**
 * @file SpeciesParser.hpp
 * @brief 物种参数解析器
 * 
 * 亚高山森林演替模型 (SFE-Model)
 * 
 * 读取 species_params.txt 文件并加载到 SpeciesManager
 * 
 * 文件格式：Tab/Space 分隔的 CSV，首行为表头
 */

#ifndef SFE_SPECIES_PARSER_HPP
#define SFE_SPECIES_PARSER_HPP

#include <string>
#include <vector>
#include <map>
#include <stdexcept>
#include <algorithm>
#include "io/FileParser.hpp"
#include "io/SpeciesProfile.hpp"

namespace sfe {

/**
 * @class SpeciesParser
 * @brief 物种参数文件解析器
 */
class SpeciesParser {
public:
    /**
     * @brief 解析物种参数文件
     * @param filepath 文件路径
     * @param species_manager 物种管理器
     * @return 加载的物种数
     */
    static int parse(const std::string& filepath, SpeciesManager& species_manager) {
        // 读取表头
        auto header = FileParser::readHeader(filepath);
        if (header.empty()) {
            throw std::runtime_error("Species file has no header: " + filepath);
        }
        
        // 建立列名到索引的映射
        std::map<std::string, int> col_map;
        for (size_t i = 0; i < header.size(); ++i) {
            std::string col_name = header[i];
            // 转换为小写以便不区分大小写匹配
            std::transform(col_name.begin(), col_name.end(), col_name.begin(), ::tolower);
            col_map[col_name] = static_cast<int>(i);
        }
        
        // 读取数据行
        auto table = FileParser::readTable(filepath, true);
        
        int species_count = 0;
        
        for (const auto& row : table) {
            if (row.empty()) continue;
            
            SpeciesProfile profile;
            
            // 解析每个字段（使用辅助函数）
            profile.species_id = getInt(row, col_map, "species_id", 0);
            profile.species_name = getString(row, col_map, "species_name", "");
            profile.is_evergreen = getBool(row, col_map, "is_evergreen", false);
            
            // 散布参数
            profile.lambda1_m = getDouble(row, col_map, "lambda1_m", 10.0);
            profile.lambda2_m = getDouble(row, col_map, "lambda2_m", 100.0);
            profile.k_ldd = getDouble(row, col_map, "k_ldd", 0.1);
            profile.fecundity_f = getDouble(row, col_map, "fecundity_f", 1.0);
            profile.maturity_age_yr = getInt(row, col_map, "maturity_age_yr", 20);
            profile.ms_seed_kg = getDouble(row, col_map, "ms_seed_kg", 0.001);
            
            // 形态参数
            profile.h_max_m = getDouble(row, col_map, "h_max_m", 25.0);
            profile.dbh_max_cm = getDouble(row, col_map, "dbh_max_cm", 50.0);
            profile.r_a = getDouble(row, col_map, "r_a", 0.5);
            profile.r_b = getDouble(row, col_map, "r_b", 0.5);
            profile.cs = getDouble(row, col_map, "cs", 1.5);
            profile.alpha_base = getDouble(row, col_map, "alpha_base", 0.3);
            
            // 异速生长参数 - 叶面积
            profile.ka1 = getDouble(row, col_map, "ka1", 0.1);
            profile.ka2 = getDouble(row, col_map, "ka2", 2.0);
            profile.kc2 = getDouble(row, col_map, "kc2", 1.0);
            
            // 异速生长参数 - 生长
            profile.opt_s = getDouble(row, col_map, "opt_s", 50.0);
            profile.ga = getDouble(row, col_map, "ga", 200.0);
            profile.gs = getDouble(row, col_map, "gs", 0.2);
            
            // 异速生长参数 - 叶生物量
            profile.dw = getDouble(row, col_map, "dw", 1.0);
            profile.fw1 = getDouble(row, col_map, "fw1", 0.1);
            profile.fw2 = getDouble(row, col_map, "fw2", 2.0);
            
            // 异速生长参数 - 幼树生物量
            profile.sa = getDouble(row, col_map, "sa", 0.1);
            profile.sb = getDouble(row, col_map, "sb", 2.5);
            
            // 异速生长参数 - 成树生物量
            profile.b0 = getDouble(row, col_map, "b0", 0.1);
            profile.b1 = getDouble(row, col_map, "b1", 2.0);
            profile.b2 = getDouble(row, col_map, "b2", 1.0);
            profile.b0_small = getDouble(row, col_map, "b0_small", 0.1);
            profile.b1_small = getDouble(row, col_map, "b1_small", 2.0);
            profile.b2_small = getDouble(row, col_map, "b2_small", 1.0);
            
            // 耐受性参数
            profile.shade_tol = getDouble(row, col_map, "shade_tol", 3.0);
            profile.dr_tol = getDouble(row, col_map, "dr_tol", 0.5);
            profile.dd_min = getDouble(row, col_map, "dd_min", 300.0);
            
            // 幼苗定居参数
            profile.l_min = getDouble(row, col_map, "l_min", 0.05);
            // 尝试多种可能的列名
            if (col_map.find("l_min") == col_map.end()) {
                profile.l_min = getDouble(row, col_map, "lmin", 0.05);
            }
            profile.l_max = getDouble(row, col_map, "l_max", 1.0);
            if (col_map.find("l_max") == col_map.end()) {
                profile.l_max = getDouble(row, col_map, "lmax", 1.0);
            }
            profile.ngs_tmin_c = getDouble(row, col_map, "ngs_tmin_c", -15.0);
            profile.ngs_tmax_c = getDouble(row, col_map, "ngs_tmax_c", 5.0);
            
            // 草本参数
            profile.k_max_herb = getDouble(row, col_map, "k_max_herb", 0.8);
            profile.r_max_herb = getDouble(row, col_map, "r_max_herb", 1.0);
            
            // 其他参数
            profile.max_age_yr = getInt(row, col_map, "max_age_yr", 200);
            profile.extinction_ke = getDouble(row, col_map, "extinction_ke", 0.5);
            
            species_manager.addSpecies(profile);
            ++species_count;
        }
        
        return species_count;
    }
    
    /**
     * @brief 便捷方法：直接从文件加载并返回 SpeciesManager
     */
    static SpeciesManager load(const std::string& filepath) {
        SpeciesManager manager;
        parse(filepath, manager);
        return manager;
    }

private:
    /**
     * @brief 从行数据中获取 double 值
     */
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
    
    /**
     * @brief 从行数据中获取 int 值
     */
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
    
    /**
     * @brief 从行数据中获取 bool 值
     */
    static bool getBool(const std::vector<std::string>& row,
                        const std::map<std::string, int>& col_map,
                        const std::string& col_name,
                        bool default_val) {
        auto it = col_map.find(col_name);
        if (it == col_map.end()) return default_val;
        
        int idx = it->second;
        if (idx < 0 || idx >= static_cast<int>(row.size())) return default_val;
        
        return FileParser::toBool(row[idx], default_val);
    }
    
    /**
     * @brief 从行数据中获取 string 值
     */
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
};

} // namespace sfe

#endif // SFE_SPECIES_PARSER_HPP
