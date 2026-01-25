/**
 * @file WaterBalance.hpp
 * @brief 水分平衡计算模块
 * 
 * 亚高山森林演替模型 (SFE-Model)
 * 
 * 实现完整的月度水分平衡计算：
 * 1. 计算样方总LAI和叶生物量
 * 2. 基于LAI修正PET
 * 3. 降雨截留计算
 * 4. 土壤入渗
 * 5. 蒸腾分配（优先序：草本 > 幼树 > 成树）
 * 6. 计算各功能群的干旱指数(DI)
 * 
 * 关键公式实现：
 * - 叶重分配比例 (50-52)
 * - 植被因子 F_vegetation (53)
 * - 截留系数 (54)
 * - 降水截留 (55-56)
 * - 林分供需 (58-59)
 * - 功能群水分竞争 (60-68)
 * - 干旱指数 (72-73)
 */

#ifndef SFE_WATER_BALANCE_HPP
#define SFE_WATER_BALANCE_HPP

#include <cmath>
#include <algorithm>
#include "core/GlobalConfig.hpp"
#include "core/MathUtils.hpp"
#include "modules/ClimateData.hpp"

namespace sfe {

// 前向声明
class Plot;
class Soil;

/**
 * @struct WaterBalanceInput
 * @brief 水分平衡计算的输入参数
 */
struct WaterBalanceInput {
    // 气候输入
    double temperature_c = 0.0;    ///< 月平均温度 (°C)
    double precipitation_mm = 0.0; ///< 月降水量 (mm)
    double pet_mm = 0.0;           ///< 月潜在蒸散发 (mm)
    
    // 植被状态
    double lai_tree = 0.0;         ///< 成树LAI
    double lai_herb = 0.0;         ///< 草本LAI
    double canopy_cover = 0.0;     ///< 冠层盖度 [0, 1]
    
    // 叶生物量 (kg)
    double foliage_tree_kg = 0.0;      ///< 成树叶生物量
    double foliage_sapling_kg = 0.0;   ///< 幼树叶生物量
    double foliage_herb_kg = 0.0;      ///< 草本叶生物量
    
    /**
     * @brief 计算样方总LAI
     */
    double getTotalLAI() const {
        return lai_tree + lai_herb;
    }
    
    /**
     * @brief 计算样方总叶生物量
     */
    double getTotalFoliage() const {
        return foliage_tree_kg + foliage_sapling_kg + foliage_herb_kg;
    }
};

/**
 * @struct WaterBalanceOutput
 * @brief 水分平衡计算的输出结果
 */
struct WaterBalanceOutput {
    // 水分通量 (mm)
    double interception_mm = 0.0;      ///< 冠层截留
    double throughfall_mm = 0.0;       ///< 穿透降水（到达土壤）
    
    // 林分总量
    double stand_demand_mm = 0.0;      ///< 林分总需求
    double stand_supply_mm = 0.0;      ///< 林分总供给
    double stand_et_mm = 0.0;          ///< 林分实际蒸腾
    double stand_aet_mm = 0.0;         ///< 林分实际蒸散发
    
    // 各功能群实际蒸腾 (mm)
    double et_herb_mm = 0.0;           ///< 草本蒸腾
    double et_sapling_mm = 0.0;        ///< 幼树蒸腾
    double et_tree_mm = 0.0;           ///< 成树蒸腾
    
    // 各功能群需求 (mm)
    double demand_herb_mm = 0.0;
    double demand_sapling_mm = 0.0;
    double demand_tree_mm = 0.0;
    
    // 各功能群供给 (mm)
    double supply_herb_mm = 0.0;
    double supply_sapling_mm = 0.0;
    double supply_tree_mm = 0.0;
    
    // 干旱指数 [0, 1]，越高越干旱
    double di_stand = 0.0;             ///< 林分干旱指数
    double di_herb = 0.0;              ///< 草本干旱指数
    double di_sapling = 0.0;           ///< 幼树干旱指数
    double di_tree = 0.0;              ///< 成树干旱指数
    
    // 土壤状态更新
    double soil_moisture_mm = 0.0;     ///< 更新后的土壤水分
    double runoff_mm = 0.0;            ///< 径流（溢出）
};

/**
 * @struct AnnualDroughtStats
 * @brief 年度干旱统计（仅温暖月份的平均）
 */
struct AnnualDroughtStats {
    double di_stand = 0.0;
    double di_herb = 0.0;
    double di_sapling = 0.0;
    double di_tree = 0.0;
    int warm_month_count = 0;          ///< 温暖月份计数（T >= 5.5°C）
    
    void reset() {
        di_stand = di_herb = di_sapling = di_tree = 0.0;
        warm_month_count = 0;
    }
    
    /**
     * @brief 累加月度干旱指数（仅温暖月份）
     */
    void accumulate(const WaterBalanceOutput& output, double temp_c) {
        if (temp_c >= GlobalConfig::T_BASE_GDD) {
            di_stand += output.di_stand;
            di_herb += output.di_herb;
            di_sapling += output.di_sapling;
            di_tree += output.di_tree;
            ++warm_month_count;
        }
    }
    
    /**
     * @brief 计算平均干旱指数
     */
    void finalize() {
        if (warm_month_count > 0) {
            double n = static_cast<double>(warm_month_count);
            di_stand /= n;
            di_herb /= n;
            di_sapling /= n;
            di_tree /= n;
        }
    }
};

/**
 * @class WaterBalanceCalculator
 * @brief 水分平衡计算器
 * 
 * 实现单月水分平衡的完整计算流程
 */
class WaterBalanceCalculator {
public:
    // 常量
    static constexpr double ALPHA_GAP = GlobalConfig::ALPHA_GAP;           ///< 林窗截留系数 (0.06)
    static constexpr double ALPHA_FOREST = GlobalConfig::ALPHA_FOREST;     ///< 林内截留系数 (0.19)
    static constexpr double MAX_ET_RATE = GlobalConfig::MAX_ET_RATE_MM;    ///< 最大蒸散率 (12 mm/月)
    
    /**
     * @brief 计算植被因子 F_vegetation
     * @param lai 样方LAI
     * @return 植被因子 [0.75, 1.0]
     * 
     * 公式(53): F_veg = 0.75 + LAI/32 (LAI <= 8), 1.0 (LAI > 8)
     */
    static double calcVegetationFactor(double lai) {
        if (lai > 8.0) return 1.0;
        return 0.75 + lai / 32.0;
    }
    
    /**
     * @brief 计算截留系数
     * @param canopy_cover 冠层盖度 [0, 1]
     * @return 截留系数
     * 
     * 公式(54): α_int = α_gap + (α_forest - α_gap) × CC
     */
    static double calcInterceptionCoefficient(double canopy_cover) {
        double cc = math::clamp01(canopy_cover);
        return ALPHA_GAP + (ALPHA_FOREST - ALPHA_GAP) * cc;
    }
    
    /**
     * @brief 计算降水截留量
     * @param precipitation 月降水量 (mm)
     * @param pet 月PET (mm)
     * @param canopy_cover 冠层盖度
     * @param lai 样方LAI
     * @return 截留量 (mm)
     * 
     * 公式(55): P_int = min(α_int × P, F_veg × PET)
     */
    static double calcInterception(double precipitation, double pet,
                                    double canopy_cover, double lai) {
        double alpha = calcInterceptionCoefficient(canopy_cover);
        double f_veg = calcVegetationFactor(lai);
        return std::min(alpha * precipitation, f_veg * pet);
    }
    
    /**
     * @brief 计算林分总水分需求
     * @param pet 月PET (mm)
     * @param interception 截留量 (mm)
     * @param lai 样方LAI
     * @return 需求量 (mm)
     * 
     * 公式(58): D_stand = F_veg × PET - P_int
     */
    static double calcStandDemand(double pet, double interception, double lai) {
        double f_veg = calcVegetationFactor(lai);
        return std::max(0.0, f_veg * pet - interception);
    }
    
    /**
     * @brief 计算林分总水分供给
     * @param soil_moisture 当前土壤水分 (mm)
     * @param whc_total 总持水量 (mm)
     * @return 供给量 (mm)
     * 
     * 公式(59): S_stand = min(c_ET, WHC_tot) × (soilmoist / WHC_tot)
     */
    static double calcStandSupply(double soil_moisture, double whc_total) {
        if (whc_total <= 0.0) return 0.0;
        double capacity = std::min(MAX_ET_RATE, whc_total);
        return capacity * (soil_moisture / whc_total);
    }
    
    /**
     * @brief 计算干旱指数
     * @param actual_et 实际蒸腾 (mm)
     * @param demand 需求 (mm)
     * @return 干旱指数 [0, 1]
     * 
     * 公式(72): DI = max(0, 1 - AT/D)
     */
    static double calcDroughtIndex(double actual_et, double demand) {
        if (demand <= 0.0) return 0.0;
        return math::clamp01(1.0 - actual_et / demand);
    }
    
    /**
     * @brief 执行完整的月度水分平衡计算
     * @param input 输入参数
     * @param soil 土壤对象（会被修改）
     * @return 输出结果
     */
    static WaterBalanceOutput calculate(const WaterBalanceInput& input, Soil& soil);
};

/**
 * @class PlotWaterBalance
 * @brief 样方级别的水分平衡管理器
 * 
 * 管理整个样方的年度水分平衡计算
 */
class PlotWaterBalance {
public:
    PlotWaterBalance() = default;
    
    /**
     * @brief 重置年度统计
     */
    void resetAnnualStats() {
        annual_stats_.reset();
        monthly_outputs_.clear();
    }
    
    /**
     * @brief 处理单月水分平衡
     * @param input 输入参数
     * @param soil 土壤对象
     * @return 月度输出
     */
    WaterBalanceOutput processMonth(const WaterBalanceInput& input, Soil& soil) {
        auto output = WaterBalanceCalculator::calculate(input, soil);
        
        // 累加到年度统计
        annual_stats_.accumulate(output, input.temperature_c);
        
        // 保存月度输出
        monthly_outputs_.push_back(output);
        
        return output;
    }
    
    /**
     * @brief 完成年度计算，获取平均干旱指数
     * @return 年度干旱统计
     */
    const AnnualDroughtStats& finalizeYear() {
        annual_stats_.finalize();
        return annual_stats_;
    }
    
    /**
     * @brief 获取当前年度统计
     */
    const AnnualDroughtStats& getAnnualStats() const {
        return annual_stats_;
    }
    
    /**
     * @brief 获取月度输出历史
     */
    const std::vector<WaterBalanceOutput>& getMonthlyOutputs() const {
        return monthly_outputs_;
    }

private:
    AnnualDroughtStats annual_stats_;
    std::vector<WaterBalanceOutput> monthly_outputs_;
};

} // namespace sfe

#endif // SFE_WATER_BALANCE_HPP
