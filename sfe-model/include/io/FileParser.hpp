/**
 * @file FileParser.hpp
 * @brief 通用文件解析工具
 * 
 * 亚高山森林演替模型 (SFE-Model)
 * 
 * 提供读取和解析 txt 文件的通用工具函数
 */

#ifndef SFE_FILE_PARSER_HPP
#define SFE_FILE_PARSER_HPP

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <cctype>
#include <iostream>

namespace sfe {

/**
 * @class FileParser
 * @brief 通用文件解析工具类
 */
class FileParser {
public:
    /**
     * @brief 去除字符串两端的空白字符
     */
    static std::string trim(const std::string& str) {
        size_t first = str.find_first_not_of(" \t\r\n");
        if (first == std::string::npos) return "";
        size_t last = str.find_last_not_of(" \t\r\n");
        return str.substr(first, last - first + 1);
    }
    
    /**
     * @brief 按分隔符分割字符串
     * @param str 输入字符串
     * @param delimiters 分隔符集合
     * @return 分割后的字符串向量
     */
    static std::vector<std::string> split(const std::string& str, 
                                           const std::string& delimiters = " \t") {
        std::vector<std::string> tokens;
        size_t start = 0;
        size_t end = 0;
        
        while ((start = str.find_first_not_of(delimiters, end)) != std::string::npos) {
            end = str.find_first_of(delimiters, start);
            tokens.push_back(str.substr(start, end - start));
        }
        
        return tokens;
    }
    
    /**
     * @brief 读取文件所有行
     * @param filepath 文件路径
     * @param skip_header 是否跳过首行（表头）
     * @return 所有行的向量
     */
    static std::vector<std::string> readLines(const std::string& filepath, 
                                               bool skip_header = true) {
        std::vector<std::string> lines;
        std::ifstream file(filepath);
        
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open file: " + filepath);
        }
        
        std::string line;
        bool first_line = true;
        
        while (std::getline(file, line)) {
            // 跳过空行
            std::string trimmed = trim(line);
            if (trimmed.empty()) continue;
            
            // 跳过注释行
            if (trimmed[0] == '#' || trimmed[0] == '/') continue;
            
            // 跳过表头
            if (first_line && skip_header) {
                first_line = false;
                continue;
            }
            first_line = false;
            
            lines.push_back(trimmed);
        }
        
        return lines;
    }
    
    /**
     * @brief 读取文件为二维字符串表格
     * @param filepath 文件路径
     * @param skip_header 是否跳过首行
     * @param delimiters 分隔符
     * @return 二维字符串向量
     */
    static std::vector<std::vector<std::string>> readTable(
        const std::string& filepath,
        bool skip_header = true,
        const std::string& delimiters = " \t"
    ) {
        std::vector<std::vector<std::string>> table;
        auto lines = readLines(filepath, skip_header);
        
        for (const auto& line : lines) {
            auto tokens = split(line, delimiters);
            if (!tokens.empty()) {
                table.push_back(tokens);
            }
        }
        
        return table;
    }
    
    /**
     * @brief 读取表头（第一行）
     * @param filepath 文件路径
     * @return 表头列名向量
     */
    static std::vector<std::string> readHeader(const std::string& filepath,
                                                const std::string& delimiters = " \t") {
        std::ifstream file(filepath);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open file: " + filepath);
        }
        
        std::string line;
        while (std::getline(file, line)) {
            std::string trimmed = trim(line);
            if (trimmed.empty() || trimmed[0] == '#' || trimmed[0] == '/') {
                continue;
            }
            return split(trimmed, delimiters);
        }
        
        return {};
    }
    
    /**
     * @brief 安全转换为 double
     */
    static double toDouble(const std::string& str, double default_val = 0.0) {
        try {
            return std::stod(trim(str));
        } catch (...) {
            return default_val;
        }
    }
    
    /**
     * @brief 安全转换为 int
     */
    static int toInt(const std::string& str, int default_val = 0) {
        try {
            return std::stoi(trim(str));
        } catch (...) {
            return default_val;
        }
    }
    
    /**
     * @brief 安全转换为 bool (支持 0/1, true/false, yes/no)
     */
    static bool toBool(const std::string& str, bool default_val = false) {
        std::string s = trim(str);
        std::transform(s.begin(), s.end(), s.begin(), ::tolower);
        
        if (s == "1" || s == "true" || s == "yes") return true;
        if (s == "0" || s == "false" || s == "no") return false;
        return default_val;
    }
    
    /**
     * @brief 检查文件是否存在
     */
    static bool fileExists(const std::string& filepath) {
        std::ifstream file(filepath);
        return file.good();
    }
};

} // namespace sfe

#endif // SFE_FILE_PARSER_HPP
