/**
 * @file Soil.hpp
 * @brief 土壤水分双层桶模型
 * 
 * 亚高山森林演替模型 (SFE-Model)
 * 
 * 实现双层土壤水分模型：
 * - 上层土壤（表层）：草本和幼树优先利用
 * - 下层土壤（次表层）：成树主要利用
 * 
 * 关键公式（参考文档3.2节）：
 * - WHC_top = min(5.6, WHC_tot)  [公式48]
 * - WHC_sub = WHC_tot - WHC_top  [公式49]
 * - 水分更新采用自下而上填充机制
 */

#ifndef SFE_SOIL_HPP
#define SFE_SOIL_HPP

#include <algorithm>
#include <cmath>
#include "core/GlobalConfig.hpp"

namespace sfe {

/**
 * @class Soil
 * @brief 双层土壤水分桶模型
 * 
 * 水分分配优先级：
 * 1. 草本抢占表层水
 * 2. 幼树抢占表层剩余水
 * 3. 成树利用深层+表层剩余水
 * 
 * 填充规则：
 * - 自下而上填充（先填满下层，溢出到上层）
 */
class Soil {
public:
    // ========================================================================
    // 构造函数
    // ========================================================================
    
    /**
     * @brief 默认构造函数
     */
    Soil() = default;
    
    /**
     * @brief 使用总持水量初始化
     * @param whc_total_mm 土壤总持水量 (mm)
     * 
     * 自动计算上下层持水量，并假设土壤初始为满水状态
     */
    explicit Soil(double whc_total_mm) {
        initialize(whc_total_mm);
    }
    
    /**
     * @brief 初始化土壤参数
     * @param whc_total_mm 土壤总持水量 (mm)
     * 
     * 文档公式(48)(49):
     * WHC_top = min(5.6cm, WHC_tot) -> 转换为mm需要*10
     * WHC_sub = WHC_tot - WHC_top
     */
    void initialize(double whc_total_mm) {
        whc_total_mm_ = std::max(0.0, whc_total_mm);
        
        // 默认表层持水量 5.6 cm = 56 mm
        constexpr double WHC_TOP_DEFAULT_MM = GlobalConfig::WHC_TOP_DEFAULT_CM * 10.0;
        
        whc_top_mm_ = std::min(WHC_TOP_DEFAULT_MM, whc_total_mm_);
        whc_sub_mm_ = whc_total_mm_ - whc_top_mm_;
        
        // 初始化为满水状态
        reset();
    }
    
    /**
     * @brief 重置为满水状态（模拟开始时或年初）
     */
    void reset() {
        sw_top_mm_ = whc_top_mm_;
        sw_sub_mm_ = whc_sub_mm_;
        soilmoist_mm_ = whc_total_mm_;
    }
    
    /**
     * @brief 火灾后重置（与reset相同，但语义更清晰）
     */
    void resetAfterFire() {
        reset();
    }
    
    // ========================================================================
    // 水分更新核心逻辑
    // ========================================================================
    
    /**
     * @brief 更新月度水分平衡（桶模型核心）
     * @param precipitation_mm 到达土壤的有效降水 (mm)
     * @param actual_et_mm 实际蒸散发总量 (mm)
     * 
     * 文档公式(71):
     * soilmoist_{m+1} = max(min(soilmoist_m + P_soil - AT_stand, WHC_tot), 0)
     * 
     * 填充规则：自下而上
     * - 先填满下层，超出部分填入上层
     * - 超出总持水量的水分溢出（径流）
     */
    void updateMonthly(double precipitation_mm, double actual_et_mm) {
        // 计算新的总土壤水分
        double new_soilmoist = soilmoist_mm_ + precipitation_mm - actual_et_mm;
        
        // 约束在 [0, WHC_tot] 范围内
        new_soilmoist = std::max(0.0, std::min(new_soilmoist, whc_total_mm_));
        
        soilmoist_mm_ = new_soilmoist;
        
        // 自下而上分配水分
        redistributeWater();
    }
    
    /**
     * @brief 桶模型更新（别名，与文档术语一致）
     */
    void bucketUpdate(double precipitation_mm, double actual_et_mm) {
        updateMonthly(precipitation_mm, actual_et_mm);
    }
    
    /**
     * @brief 计算月初表层水分
     * @return 表层当前水分 (mm)
     * 
     * 文档公式(57):
     * Water_top = max(soilmoist - WHC_sub, 0)
     */
    double calcTopWater() const {
        return std::max(soilmoist_mm_ - whc_sub_mm_, 0.0);
    }
    
    /**
     * @brief 从表层消耗水分
     * @param amount_mm 消耗量 (mm)
     * @return 实际消耗量 (mm)
     */
    double consumeFromTop(double amount_mm) {
        double available = sw_top_mm_;
        double consumed = std::min(amount_mm, available);
        sw_top_mm_ -= consumed;
        soilmoist_mm_ -= consumed;
        return consumed;
    }
    
    /**
     * @brief 从下层消耗水分
     * @param amount_mm 消耗量 (mm)
     * @return 实际消耗量 (mm)
     */
    double consumeFromSub(double amount_mm) {
        double available = sw_sub_mm_;
        double consumed = std::min(amount_mm, available);
        sw_sub_mm_ -= consumed;
        soilmoist_mm_ -= consumed;
        return consumed;
    }
    
    /**
     * @brief 从任意层消耗水分（先表层后下层）
     * @param amount_mm 消耗量 (mm)
     * @return 实际消耗量 (mm)
     */
    double consume(double amount_mm) {
        double from_top = consumeFromTop(amount_mm);
        double remaining = amount_mm - from_top;
        double from_sub = consumeFromSub(remaining);
        return from_top + from_sub;
    }
    
    /**
     * @brief 添加水分（降水入渗）
     * @param amount_mm 入渗量 (mm)
     * @return 溢出量（径流）(mm)
     */
    double addWater(double amount_mm) {
        double new_total = soilmoist_mm_ + amount_mm;
        double overflow = 0.0;
        
        if (new_total > whc_total_mm_) {
            overflow = new_total - whc_total_mm_;
            new_total = whc_total_mm_;
        }
        
        soilmoist_mm_ = new_total;
        redistributeWater();
        
        return overflow;
    }
    
    // ========================================================================
    // 水分供给计算
    // ========================================================================
    
    /**
     * @brief 计算林分总水分供给能力
     * @return 供给能力 (mm)
     * 
     * 文档公式(59):
     * S_stand = min(c_ET, WHC_tot) * (soilmoist / WHC_tot)
     * c_ET = 12 mm/月
     */
    double calcStandSupply() const {
        if (whc_total_mm_ <= 0.0) return 0.0;
        
        constexpr double MAX_ET_RATE = GlobalConfig::MAX_ET_RATE_MM;
        double supply_capacity = std::min(MAX_ET_RATE, whc_total_mm_);
        return supply_capacity * (soilmoist_mm_ / whc_total_mm_);
    }
    
    /**
     * @brief 计算草本水分供给
     * @param p_soil_mm 到达土壤的降水 (mm)
     * @param stand_supply_mm 林分总供给 (mm)
     * @return 草本可用供给 (mm)
     * 
     * 文档公式(60):
     * S_herb = min(0.5 * P_soil + Water_top, S_stand)
     */
    double calcHerbSupply(double p_soil_mm, double stand_supply_mm) const {
        double top_water = calcTopWater();
        return std::min(0.5 * p_soil_mm + top_water, stand_supply_mm);
    }
    
    // ========================================================================
    // Getter 方法
    // ========================================================================
    
    /// 获取总持水量 (mm)
    double getWhcTotal() const { return whc_total_mm_; }
    
    /// 获取上层持水量 (mm)
    double getWhcTop() const { return whc_top_mm_; }
    
    /// 获取下层持水量 (mm)
    double getWhcSub() const { return whc_sub_mm_; }
    
    /// 获取当前总土壤水分 (mm)
    double getSoilMoist() const { return soilmoist_mm_; }
    
    /// 获取当前上层水分 (mm)
    double getSwTop() const { return sw_top_mm_; }
    
    /// 获取当前下层水分 (mm)
    double getSwSub() const { return sw_sub_mm_; }
    
    /// 获取相对土壤湿度 [0, 1]
    double getRelativeMoisture() const {
        if (whc_total_mm_ <= 0.0) return 0.0;
        return soilmoist_mm_ / whc_total_mm_;
    }
    
    /// 检查是否处于干旱状态
    bool isDry(double threshold = 0.3) const {
        return getRelativeMoisture() < threshold;
    }

private:
    // 持水量参数（常量，初始化后不变）
    double whc_total_mm_ = 0.0;    ///< 总持水量 (mm)
    double whc_top_mm_ = 0.0;      ///< 上层持水量 (mm)
    double whc_sub_mm_ = 0.0;      ///< 下层持水量 (mm)
    
    // 当前水分状态（动态变化）
    double soilmoist_mm_ = 0.0;    ///< 当前总土壤水分 (mm)
    double sw_top_mm_ = 0.0;       ///< 当前上层水分 (mm)
    double sw_sub_mm_ = 0.0;       ///< 当前下层水分 (mm)
    
    /**
     * @brief 重新分配水分到上下层（自下而上填充）
     * 
     * 溢出机制：
     * - 下层先填满
     * - 超出下层容量的部分进入上层
     */
    void redistributeWater() {
        if (soilmoist_mm_ <= whc_sub_mm_) {
            // 总水分不超过下层容量，全部在下层
            sw_sub_mm_ = soilmoist_mm_;
            sw_top_mm_ = 0.0;
        } else {
            // 下层填满，剩余进入上层
            sw_sub_mm_ = whc_sub_mm_;
            sw_top_mm_ = soilmoist_mm_ - whc_sub_mm_;
        }
        
        // 确保不超过各层容量
        sw_top_mm_ = std::min(sw_top_mm_, whc_top_mm_);
        sw_sub_mm_ = std::min(sw_sub_mm_, whc_sub_mm_);
    }
};

} // namespace sfe

#endif // SFE_SOIL_HPP
