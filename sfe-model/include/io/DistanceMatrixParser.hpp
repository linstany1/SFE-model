/**
 * @file DistanceMatrixParser.hpp
 * @brief 距离矩阵解析器
 * 
 * 亚高山森林演替模型 (SFE-Model)
 * 
 * Task 2.4: 解析 distance_matrix.txt
 * 
 * 关键点：这是非对称矩阵！
 * - 行 (Rows): 代表所有样方 (Sources)，包括模拟样方 (0-2) 和纯种源样方 (3-14)
 * - 列 (Columns): 代表模拟样方 (Target Plots)，例如 0, 1, 2
 * - D[i][j] 表示从 Source i 到 Target j 的距离
 * 
 * 文件格式：
 * 第一行可能是表头（目标样方ID）
 * 后续行：第一列可能是源样方ID，其余列是距离值
 */

#ifndef SFE_DISTANCE_MATRIX_PARSER_HPP
#define SFE_DISTANCE_MATRIX_PARSER_HPP

#include <string>
#include <vector>
#include <stdexcept>
#include "io/FileParser.hpp"
#include "core/Types.hpp"

namespace sfe {

/**
 * @class DistanceMatrix
 * @brief 距离矩阵类
 * 
 * 存储源样方到目标样方的距离
 */
class DistanceMatrix {
public:
    DistanceMatrix() = default;
    
    /**
     * @brief 使用指定维度初始化
     * @param num_sources 源样方数量（行数）
     * @param num_targets 目标样方数量（列数）
     */
    DistanceMatrix(int num_sources, int num_targets)
        : num_sources_(num_sources)
        , num_targets_(num_targets)
        , data_(num_sources, std::vector<double>(num_targets, 0.0))
    {
    }
    
    /**
     * @brief 设置距离值
     * @param source_id 源样方ID（行索引）
     * @param target_id 目标样方ID（列索引）
     * @param distance 距离 (m)
     */
    void setDistance(int source_id, int target_id, double distance) {
        if (source_id >= 0 && source_id < num_sources_ &&
            target_id >= 0 && target_id < num_targets_) {
            data_[source_id][target_id] = distance;
        }
    }
    
    /**
     * @brief 获取距离值
     * @param source_id 源样方ID
     * @param target_id 目标样方ID
     * @return 距离 (m)
     */
    double getDistance(int source_id, int target_id) const {
        if (source_id >= 0 && source_id < num_sources_ &&
            target_id >= 0 && target_id < num_targets_) {
            return data_[source_id][target_id];
        }
        return 0.0;
    }
    
    /**
     * @brief 获取从指定源到所有目标的距离
     */
    const std::vector<double>& getDistancesFromSource(int source_id) const {
        static const std::vector<double> empty;
        if (source_id >= 0 && source_id < num_sources_) {
            return data_[source_id];
        }
        return empty;
    }
    
    /**
     * @brief 获取源样方数量
     */
    int getNumSources() const { return num_sources_; }
    
    /**
     * @brief 获取目标样方数量
     */
    int getNumTargets() const { return num_targets_; }
    
    /**
     * @brief 检查是否有效
     */
    bool isValid() const {
        return num_sources_ > 0 && num_targets_ > 0;
    }
    
    /**
     * @brief 获取原始数据（用于调试）
     */
    const std::vector<std::vector<double>>& getData() const { return data_; }
    
    /**
     * @brief 调整大小
     */
    void resize(int num_sources, int num_targets) {
        num_sources_ = num_sources;
        num_targets_ = num_targets;
        data_.resize(num_sources);
        for (auto& row : data_) {
            row.resize(num_targets, 0.0);
        }
    }

private:
    int num_sources_ = 0;
    int num_targets_ = 0;
    std::vector<std::vector<double>> data_;
};

/**
 * @class DistanceMatrixParser
 * @brief 距离矩阵文件解析器
 */
class DistanceMatrixParser {
public:
    /**
     * @brief 解析距离矩阵文件
     * @param filepath 文件路径
     * @param has_header 是否有表头
     * @param has_row_id 第一列是否是行ID
     * @return DistanceMatrix 对象
     */
    static DistanceMatrix parse(const std::string& filepath, 
                                 bool has_header = true,
                                 bool has_row_id = true) {
        std::vector<std::string> lines = FileParser::readLines(filepath, has_header);
        
        if (lines.empty()) {
            throw std::runtime_error("Distance matrix file is empty: " + filepath);
        }
        
        // 解析第一行确定列数
        auto first_row = FileParser::split(lines[0]);
        int num_targets = static_cast<int>(first_row.size());
        if (has_row_id) {
            num_targets -= 1;  // 第一列是行ID
        }
        
        int num_sources = static_cast<int>(lines.size());
        
        DistanceMatrix matrix(num_sources, num_targets);
        
        for (int i = 0; i < num_sources; ++i) {
            auto tokens = FileParser::split(lines[i]);
            
            int col_offset = has_row_id ? 1 : 0;
            
            for (int j = 0; j < num_targets; ++j) {
                int col_idx = col_offset + j;
                if (col_idx < static_cast<int>(tokens.size())) {
                    double dist = FileParser::toDouble(tokens[col_idx], 0.0);
                    matrix.setDistance(i, j, dist);
                }
            }
        }
        
        return matrix;
    }
    
    /**
     * @brief 自动检测格式并解析
     * @param filepath 文件路径
     * @return DistanceMatrix 对象
     */
    static DistanceMatrix parseAuto(const std::string& filepath) {
        // 读取第一行（可能是表头）
        std::ifstream file(filepath);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open file: " + filepath);
        }
        
        std::string first_line;
        while (std::getline(file, first_line)) {
            first_line = FileParser::trim(first_line);
            if (!first_line.empty() && first_line[0] != '#') break;
        }
        file.close();
        
        // 检测是否有表头（第一个字段是否为非数字）
        auto tokens = FileParser::split(first_line);
        bool has_header = false;
        if (!tokens.empty()) {
            try {
                std::stod(tokens[0]);
            } catch (...) {
                has_header = true;  // 第一个字段不是数字，可能是表头
            }
        }
        
        // 检测是否有行ID
        // 假设如果第一列的值是递增的整数序列，则可能是行ID
        auto lines = FileParser::readLines(filepath, has_header);
        bool has_row_id = false;
        
        if (!lines.empty()) {
            auto first_tokens = FileParser::split(lines[0]);
            if (!first_tokens.empty()) {
                // 检查第一列是否看起来像ID
                try {
                    int first_val = std::stoi(first_tokens[0]);
                    // 如果第一个值是 0 或 1，可能是行ID
                    if (first_val == 0 || first_val == 1) {
                        has_row_id = true;
                    }
                } catch (...) {
                    has_row_id = false;
                }
            }
        }
        
        return parse(filepath, has_header, has_row_id);
    }
};

} // namespace sfe

#endif // SFE_DISTANCE_MATRIX_PARSER_HPP
