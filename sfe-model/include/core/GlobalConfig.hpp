/**
 * @file GlobalConfig.hpp
 * @brief 全局配置结构体，存放模拟参数和系统配置
 * 
 * 亚高山森林演替模型 (SFE-Model)
 * 
 * 该结构体存储所有模拟参数，包括：
 * - 模拟年限和时间范围
 * - 输入/输出路径
 * - 随机种子（确保可复现性）
 * - 模型常量和阈值
 */

#ifndef SFE_GLOBAL_CONFIG_HPP
#define SFE_GLOBAL_CONFIG_HPP

#include <string>
#include <cstdint>
#include <vector>

namespace sfe {

/**
 * @struct GlobalConfig
 * @brief 存储模拟的全局配置参数
 */
struct GlobalConfig {
    // ========================================================================
    // 时间配置
    // ========================================================================
    int start_year = 1950;           ///< 模拟开始年份
    int end_year = 2100;             ///< 模拟结束年份
    int output_interval = 10;        ///< 输出间隔（年），用于详细个体输出
    
    // ========================================================================
    // 随机数配置（确保确定性可复现）
    // ========================================================================
    uint64_t random_seed = 42;       ///< 随机数种子
    
    // ========================================================================
    // 文件路径配置
    // ========================================================================
    std::string input_dir = "./data/";           ///< 输入数据目录
    std::string output_dir = "./output/";        ///< 输出数据目录
    
    // 输入文件名
    std::string site_file = "site.txt";                           ///< 样地定义文件
    std::string species_params_file = "species_params.txt";       ///< 物种参数文件
    std::string climate_file = "climate.txt";                     ///< 气候数据文件
    std::string distance_matrix_file = "distance_matrix.txt";     ///< 距离矩阵文件
    std::string seed_source_file = "seed_source_composition.txt"; ///< 种源组成文件
    
    // ========================================================================
    // 空间配置
    // ========================================================================
    static constexpr int PLOT_SIZE_M = 30;        ///< 样方大小（米）
    static constexpr int CELL_SIZE_M = 1;         ///< 网格单元大小（米）
    static constexpr int CELLS_PER_PLOT = PLOT_SIZE_M * PLOT_SIZE_M; ///< 每样方网格数
    static constexpr double PLOT_AREA_M2 = PLOT_SIZE_M * PLOT_SIZE_M; ///< 样方面积(m²)
    
    // 空间哈希网格配置（用于光照计算优化）
    static constexpr int HASH_GRID_SIZE_M = 5;    ///< 空间哈希网格大小（米）
    static constexpr int HASH_GRID_DIM = PLOT_SIZE_M / HASH_GRID_SIZE_M; ///< 哈希网格维度
    
    // ========================================================================
    // 物理/生态常量
    // ========================================================================
    
    // 高度阈值
    static constexpr double HEIGHT_THRESHOLD_M = 1.37;   ///< 幼树/成树划分高度阈值(m)
    static constexpr double HERB_HEIGHT_M = 1.0;         ///< 草本层固定高度(m)
    static constexpr double SEEDLING_INIT_HEIGHT_M = 0.05; ///< 幼苗初始高度(m)
    
    // 温度基准
    static constexpr double T_BASE_GDD = 5.5;     ///< 生长度日计算基准温度(°C)
    static constexpr double DAYS_PER_MONTH = 30.5; ///< 每月平均天数
    
    // 水分相关
    static constexpr double WHC_TOP_DEFAULT_CM = 5.6;   ///< 默认表层土壤持水量(cm)
    static constexpr double MAX_ET_RATE_MM = 12.0;      ///< 最大蒸散速率(mm/月)
    static constexpr double ALPHA_GAP = 0.06;           ///< 林窗截留系数
    static constexpr double ALPHA_FOREST = 0.19;        ///< 林内截留系数
    
    // 光照相关
    static constexpr double EXTINCTION_KE_DEFAULT = 0.5; ///< 默认消光系数
    static constexpr double LAI_SATURATION = 3.0;       ///< LAI饱和阈值（种源强度）
    // Phase 9: PRUNE_THRESHOLD 已移除，改用物种特异性参数 l_min
    
    // 种子相关
    static constexpr double SEED_UNLIMITED_THRESHOLD = 100.0; ///< 种子饱和阈值(seeds/m²)
    static constexpr int SEED_BANK_MAX_AGE = 2;         ///< 种子库最大存活年数
    
    // 草本相关
    static constexpr double HERB_MIN_BIOMASS_KG = 0.02; ///< 草本最小生物量(kg/m²)
    static constexpr double HERB_LA_COEFFICIENT = 5.84; ///< 草本叶面积系数
    static constexpr double HERB_FOLIAGE_FRACTION = 0.5; ///< 草本叶生物量比例
    
    // 幼树叶生物量比例
    static constexpr double SAPLING_FOLIAGE_EVERGREEN = 0.4;  ///< 常绿幼树叶比例
    static constexpr double SAPLING_FOLIAGE_DECIDUOUS = 0.2;  ///< 落叶幼树叶比例
    
    // 死亡率相关
    static constexpr double STRESS_THRESHOLD = 0.1;     ///< 胁迫死亡阈值 f_env
    static constexpr int STRESS_YEARS_THRESHOLD = 3;    ///< 胁迫死亡年数阈值
    static constexpr double LONGEVITY_SURVIVAL = 0.01;  ///< 存活到最大年龄的概率(1%)
    
    // 冠层几何
    static constexpr double MIN_CROWN_FRACTION = 0.1;   ///< 最小冠层比例（树高的10%）
    static constexpr double EFFECTIVE_LIGHT_FRACTION = 0.75; ///< 有效受光冠层比例
    
    // ========================================================================
    // 干扰配置
    // ========================================================================
    std::vector<int> fire_years;     ///< 火灾发生年份列表
    
    // ========================================================================
    // 输出开关
    // ========================================================================
    bool output_plot_summary = true;          ///< 输出样方汇总
    bool output_species_composition = true;   ///< 输出物种组成
    bool output_individual_trees = true;      ///< 输出个体树木详情
    bool output_disturbance_log = true;       ///< 输出干扰日志
    bool output_debug_seed = false;           ///< 输出种子调试信息
    
    // ========================================================================
    // 方法
    // ========================================================================
    
    /**
     * @brief 获取完整的输入文件路径
     * @param filename 文件名
     * @return 完整路径
     */
    std::string getInputPath(const std::string& filename) const {
        return input_dir + filename;
    }
    
    /**
     * @brief 获取完整的输出文件路径
     * @param filename 文件名
     * @return 完整路径
     */
    std::string getOutputPath(const std::string& filename) const {
        return output_dir + filename;
    }
    
    /**
     * @brief 计算模拟总年数
     * @return 模拟年数
     */
    int getSimulationYears() const {
        return end_year - start_year + 1;
    }
    
    /**
     * @brief 检查指定年份是否为火灾年
     * @param year 年份
     * @return 是否发生火灾
     */
    bool isFireYear(int year) const {
        for (int fy : fire_years) {
            if (fy == year) return true;
        }
        return false;
    }
    
    /**
     * @brief 检查是否应该输出详细个体数据
     * @param year_index 模拟年索引（从0开始）
     * @return 是否输出
     */
    bool shouldOutputIndividuals(int year_index) const {
        return output_individual_trees && (year_index % output_interval == 0);
    }
};

/**
 * @brief 全局配置实例（单例模式访问点）
 * 
 * 使用方法：
 * @code
 * auto& config = sfe::getGlobalConfig();
 * config.random_seed = 12345;
 * @endcode
 */
inline GlobalConfig& getGlobalConfig() {
    static GlobalConfig instance;
    return instance;
}

} // namespace sfe

#endif // SFE_GLOBAL_CONFIG_HPP
