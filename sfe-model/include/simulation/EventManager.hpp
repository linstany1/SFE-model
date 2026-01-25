/**
 * @file EventManager.hpp
 * @brief 干扰事件管理器
 * 
 * 亚高山森林演替模型 (SFE-Model)
 * Phase 8: 实现基于事件的火灾干扰机制
 * 
 * 读取 events.txt 并提供查询接口
 */

#ifndef SFE_EVENT_MANAGER_HPP
#define SFE_EVENT_MANAGER_HPP

#include <string>
#include <vector>
#include <map>
#include <set>
#include <fstream>
#include <sstream>
#include <algorithm>
#include "core/Types.hpp"

namespace sfe {

/**
 * @struct DisturbanceEvent
 * @brief 干扰事件结构体
 * 
 * 对应 events.txt 中的一行
 */
struct DisturbanceEvent {
    Year year = 0;                    ///< 干扰发生年份
    int plot_id = -1;                 ///< 受干扰的目标样方 ID（-1 表示所有样方）
    std::string disturbance_type = "Fire";  ///< 干扰类型（如 "Fire", "Insect"）
    double intensity = 1.0;           ///< 干扰强度 (0.0-1.0)
    
    DisturbanceEvent() = default;
    
    DisturbanceEvent(Year y, int pid, const std::string& type, double intens)
        : year(y), plot_id(pid), disturbance_type(type), intensity(intens)
    {}
    
    /**
     * @brief 是否为全局事件（应用于所有样方）
     */
    bool isGlobal() const {
        return plot_id == -1;
    }
    
    /**
     * @brief 是否是火灾事件
     */
    bool isFire() const {
        return disturbance_type == "Fire" || disturbance_type == "fire";
    }
};

/**
 * @class EventManager
 * @brief 干扰事件管理器
 * 
 * 管理所有样方的干扰事件队列，支持查询"某年某样方"是否发生干扰
 */
class EventManager {
public:
    EventManager() = default;
    
    // ========================================================================
    // 数据加载
    // ========================================================================
    
    /**
     * @brief 从文件加载事件
     * @param filepath events.txt 文件路径
     * @return 加载的事件数量
     */
    int loadFromFile(const std::string& filepath) {
        std::ifstream file(filepath);
        if (!file.is_open()) {
            // 文件不存在是允许的（无干扰模拟）
            return 0;
        }
        
        events_.clear();
        events_by_year_.clear();
        
        std::string line;
        bool first_line = true;
        int count = 0;
        
        while (std::getline(file, line)) {
            // 跳过空行和注释
            if (line.empty() || line[0] == '#') continue;
            
            // 跳过表头
            if (first_line) {
                first_line = false;
                // 检查是否是数据行（首字符是数字）
                if (!std::isdigit(line[0]) && line[0] != '-') continue;
            }
            
            std::istringstream iss(line);
            DisturbanceEvent event;
            
            // 解析：year, plot_id, disturbance_type, intensity
            std::string token;
            std::vector<std::string> tokens;
            while (iss >> token) {
                tokens.push_back(token);
            }
            
            if (tokens.size() >= 4) {
                try {
                    event.year = std::stoi(tokens[0]);
                    event.plot_id = std::stoi(tokens[1]);
                    event.disturbance_type = tokens[2];
                    event.intensity = std::stod(tokens[3]);
                    
                    addEvent(event);
                    count++;
                } catch (...) {
                    // 解析错误，跳过此行
                    continue;
                }
            }
        }
        
        return count;
    }
    
    /**
     * @brief 手动添加事件
     */
    void addEvent(const DisturbanceEvent& event) {
        events_.push_back(event);
        events_by_year_[event.year].push_back(events_.size() - 1);
    }
    
    /**
     * @brief 清空所有事件
     */
    void clear() {
        events_.clear();
        events_by_year_.clear();
    }
    
    // ========================================================================
    // 查询接口
    // ========================================================================
    
    /**
     * @brief 检查某年某样方是否发生火灾
     * @param year 年份
     * @param plot_id 样方ID
     * @return 如果发生火灾返回 true
     */
    bool hasFireEvent(Year year, PlotId plot_id) const {
        auto it = events_by_year_.find(year);
        if (it == events_by_year_.end()) return false;
        
        for (size_t idx : it->second) {
            const auto& event = events_[idx];
            if (event.isFire()) {
                // 检查是否匹配：全局事件或特定样方
                if (event.isGlobal() || event.plot_id == static_cast<int>(plot_id)) {
                    return true;
                }
            }
        }
        return false;
    }
    
    /**
     * @brief 获取某年某样方的火灾事件（如果有）
     * @param year 年份
     * @param plot_id 样方ID
     * @return 指向事件的指针，若无则返回 nullptr
     */
    const DisturbanceEvent* getFireEvent(Year year, PlotId plot_id) const {
        auto it = events_by_year_.find(year);
        if (it == events_by_year_.end()) return nullptr;
        
        for (size_t idx : it->second) {
            const auto& event = events_[idx];
            if (event.isFire()) {
                if (event.isGlobal() || event.plot_id == static_cast<int>(plot_id)) {
                    return &event;
                }
            }
        }
        return nullptr;
    }
    
    /**
     * @brief 获取某年发生干扰的所有样方 ID
     * @param year 年份
     * @return 样方 ID 集合（-1 表示全局）
     */
    std::set<int> getAffectedPlots(Year year) const {
        std::set<int> result;
        auto it = events_by_year_.find(year);
        if (it == events_by_year_.end()) return result;
        
        for (size_t idx : it->second) {
            result.insert(events_[idx].plot_id);
        }
        return result;
    }
    
    /**
     * @brief 检查某年是否有任何干扰事件
     */
    bool hasAnyEvent(Year year) const {
        return events_by_year_.find(year) != events_by_year_.end();
    }
    
    /**
     * @brief 获取所有事件
     */
    const std::vector<DisturbanceEvent>& getAllEvents() const {
        return events_;
    }
    
    /**
     * @brief 获取事件总数
     */
    size_t getEventCount() const {
        return events_.size();
    }
    
    /**
     * @brief 获取有事件发生的年份列表
     */
    std::vector<Year> getEventYears() const {
        std::vector<Year> years;
        for (const auto& [year, _] : events_by_year_) {
            years.push_back(year);
        }
        std::sort(years.begin(), years.end());
        return years;
    }

private:
    std::vector<DisturbanceEvent> events_;                  ///< 所有事件列表
    std::map<Year, std::vector<size_t>> events_by_year_;   ///< 按年份索引的事件
};

/**
 * @class EventParser
 * @brief 事件文件解析器（静态工具类）
 */
class EventParser {
public:
    /**
     * @brief 解析 events.txt 并返回 EventManager
     */
    static EventManager load(const std::string& filepath) {
        EventManager manager;
        manager.loadFromFile(filepath);
        return manager;
    }
};

} // namespace sfe

#endif // SFE_EVENT_MANAGER_HPP
