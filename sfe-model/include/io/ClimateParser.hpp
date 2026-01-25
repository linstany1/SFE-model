/**
 * @file ClimateParser.hpp
 * @brief 气候数据解析器
 * 
 * 亚高山森林演替模型 (SFE-Model)
 * 
 * 读取 climate.txt 文件并加载到 ClimateManager
 * 
 * 文件格式：
 * Year + Temp_Jan...Temp_Dec + Prec_Jan...Prec_Dec + PET_Jan...PET_Dec
 * 共 37 列，每行代表一年
 */

#ifndef SFE_CLIMATE_PARSER_HPP
#define SFE_CLIMATE_PARSER_HPP

#include <string>
#include <vector>
#include <stdexcept>
#include "io/FileParser.hpp"
#include "modules/ClimateData.hpp"

namespace sfe {

/**
 * @class ClimateParser
 * @brief 气候数据文件解析器
 */
class ClimateParser {
public:
    /**
     * @brief 解析气候文件并加载到 ClimateManager
     * @param filepath 文件路径
     * @param climate_manager 气候管理器
     * @return 加载的记录数
     * 
     * Phase 8 更新：支持两种格式
     * 1. 新格式：8 列（climate_id, year, month, Temp, Temp_min, Temp_max, Prec, PET）
     * 2. 旧格式：7 列（year, month, Temp, Temp_min, Temp_max, Prec, PET）
     */
    static int parse(const std::string& filepath, ClimateManager& climate_manager) {
        auto table = FileParser::readTable(filepath, true);
        
        if (table.empty()) {
            throw std::runtime_error("Climate file is empty: " + filepath);
        }
        
        // 检测格式：8 列为新格式，7 列为旧格式
        bool has_climate_id = (table[0].size() >= 8);
        
        std::map<std::pair<int, Year>, YearlyClimate> year_data;
        int records = 0;
        
        for (const auto& row : table) {
            int climate_id = 0;
            Year year;
            int month;
            double temp, temp_min, temp_max, prec, pet;
            
            if (has_climate_id && row.size() >= 8) {
                // 新格式：climate_id, year, month, Temp, Temp_min, Temp_max, Prec, PET
                climate_id = FileParser::toInt(row[0]);
                year = FileParser::toInt(row[1]);
                month = FileParser::toInt(row[2]);
                temp = FileParser::toDouble(row[3]);
                temp_min = FileParser::toDouble(row[4]);
                temp_max = FileParser::toDouble(row[5]);
                prec = FileParser::toDouble(row[6]);
                pet = FileParser::toDouble(row[7]);
            } else if (row.size() >= 7) {
                // 旧格式：year, month, Temp, Temp_min, Temp_max, Prec, PET
                climate_id = 0;  // 默认 climate_id
                year = FileParser::toInt(row[0]);
                month = FileParser::toInt(row[1]);
                temp = FileParser::toDouble(row[2]);
                temp_min = FileParser::toDouble(row[3]);
                temp_max = FileParser::toDouble(row[4]);
                prec = FileParser::toDouble(row[5]);
                pet = FileParser::toDouble(row[6]);
            } else {
                continue;  // 跳过格式不正确的行
            }
            
            if (month < 1 || month > 12) continue;
            
            auto key = std::make_pair(climate_id, year);
            if (year_data.find(key) == year_data.end()) {
                year_data[key].year = year;
            }
            
            auto& monthly = year_data[key].months[month - 1];
            monthly.temp = temp;
            monthly.temp_min = temp_min;
            monthly.temp_max = temp_max;
            monthly.prec = prec;
            monthly.pet = pet;
            
            records++;
        }
        
        // 添加到管理器（支持多气候序列）
        for (auto& [key, yearly] : year_data) {
            climate_manager.addYearlyClimate(key.first, yearly);
        }
        
        return records;
    }
    
    /**
     * @brief 便捷方法：直接从文件加载并返回 ClimateManager
     */
    static ClimateManager load(const std::string& filepath) {
        ClimateManager manager;
        parse(filepath, manager);
        return manager;
    }

private:
    /**
     * @brief 解析旧格式（每行一个月）
     * 格式：year, month, temp, prec, pet
     */
    static void parseOldFormat(const std::vector<std::vector<std::string>>& table,
                                ClimateManager& climate_manager) {
        std::map<Year, YearlyClimate> year_data;
        
        for (const auto& row : table) {
            if (row.size() < 5) continue;
            
            Year year = FileParser::toInt(row[0]);
            int month = FileParser::toInt(row[1]);  // 1-12
            
            if (month < 1 || month > 12) continue;
            
            if (year_data.find(year) == year_data.end()) {
                year_data[year].year = year;
            }
            
            auto& monthly = year_data[year].months[month - 1];
            monthly.temp = FileParser::toDouble(row[2]);
            monthly.prec = FileParser::toDouble(row[3]);
            monthly.pet = FileParser::toDouble(row[4]);
        }
        
        for (auto& [year, yearly] : year_data) {
            climate_manager.addYearlyClimate(yearly);
        }
    }
};

} // namespace sfe

#endif // SFE_CLIMATE_PARSER_HPP
