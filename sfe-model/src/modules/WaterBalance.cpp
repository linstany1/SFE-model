/**
 * @file WaterBalance.cpp
 * @brief 水分平衡计算实现
 * 
 * 亚高山森林演替模型 (SFE-Model)
 */

#include "modules/WaterBalance.hpp"
#include "entities/Soil.hpp"
#include <algorithm>
#include <cmath>

namespace sfe {

WaterBalanceOutput WaterBalanceCalculator::calculate(
    const WaterBalanceInput& input, 
    Soil& soil
) {
    WaterBalanceOutput output;
    
    // ========================================================================
    // Step 1: 计算基本参数
    // ========================================================================
    
    double lai_total = input.getTotalLAI();
    double foliage_total = input.getTotalFoliage();
    
    // 计算叶重分配比例 (公式50-52)
    double f_herb = 0.0, f_sapling = 0.0, f_tree = 0.0;
    if (foliage_total > 0.0) {
        f_herb = input.foliage_herb_kg / foliage_total;
        f_sapling = input.foliage_sapling_kg / foliage_total;
        f_tree = input.foliage_tree_kg / foliage_total;
    } else {
        // 无植被时，假设全为草本需求
        f_herb = 1.0;
    }
    
    // ========================================================================
    // Step 2: 降水截留与入渗 (公式54-56)
    // ========================================================================
    
    output.interception_mm = calcInterception(
        input.precipitation_mm, 
        input.pet_mm,
        input.canopy_cover, 
        lai_total
    );
    
    output.throughfall_mm = input.precipitation_mm - output.interception_mm;
    double p_soil = output.throughfall_mm;  // 到达土壤的水
    
    // ========================================================================
    // Step 3: 计算林分总需求和供给 (公式58-59)
    // ========================================================================
    
    output.stand_demand_mm = calcStandDemand(
        input.pet_mm, 
        output.interception_mm, 
        lai_total
    );
    
    output.stand_supply_mm = calcStandSupply(
        soil.getSoilMoist(), 
        soil.getWhcTotal()
    );
    
    // ========================================================================
    // Step 4: 各功能群需求分配
    // ========================================================================
    
    output.demand_herb_mm = output.stand_demand_mm * f_herb;
    output.demand_sapling_mm = output.stand_demand_mm * f_sapling;
    output.demand_tree_mm = output.stand_demand_mm * f_tree;
    
    // ========================================================================
    // Step 5: 水分竞争分配（优先序：草本 -> 幼树 -> 成树）
    // ========================================================================
    
    // 计算月初表层水分 (公式57)
    double water_top = soil.calcTopWater();
    
    // 草本供给和实际蒸腾 (公式60-62)
    // S_herb = min(0.5 * P_soil + Water_top, S_stand)
    output.supply_herb_mm = std::min(
        0.5 * p_soil + water_top, 
        output.stand_supply_mm
    );
    output.et_herb_mm = std::min(output.supply_herb_mm, output.demand_herb_mm);
    
    // 幼树供给和实际蒸腾 (公式63-65)
    // S_sapling = S_herb - AT_herb
    output.supply_sapling_mm = output.supply_herb_mm - output.et_herb_mm;
    output.et_sapling_mm = std::min(output.supply_sapling_mm, output.demand_sapling_mm);
    
    // 成树供给和实际蒸腾 (公式66-68)
    // S_tree = S_stand - (AT_herb + AT_sapling)
    output.supply_tree_mm = output.stand_supply_mm - 
                            (output.et_herb_mm + output.et_sapling_mm);
    output.supply_tree_mm = std::max(0.0, output.supply_tree_mm);
    output.et_tree_mm = std::min(output.supply_tree_mm, output.demand_tree_mm);
    
    // ========================================================================
    // Step 6: 林分总蒸腾和蒸散发 (公式69-70)
    // ========================================================================
    
    output.stand_et_mm = std::min(output.stand_supply_mm, output.stand_demand_mm);
    output.stand_aet_mm = output.stand_et_mm + output.interception_mm;
    
    // ========================================================================
    // Step 7: 更新土壤水分 (公式71)
    // ========================================================================
    
    // soilmoist_{m+1} = max(min(soilmoist + P_soil - AT_stand, WHC_tot), 0)
    double new_moisture = soil.getSoilMoist() + p_soil - output.stand_et_mm;
    
    // 计算径流（超出持水量的部分）
    if (new_moisture > soil.getWhcTotal()) {
        output.runoff_mm = new_moisture - soil.getWhcTotal();
        new_moisture = soil.getWhcTotal();
    } else if (new_moisture < 0.0) {
        new_moisture = 0.0;
    }
    
    // 更新土壤状态
    soil.bucketUpdate(p_soil, output.stand_et_mm);
    output.soil_moisture_mm = soil.getSoilMoist();
    
    // ========================================================================
    // Step 8: 计算干旱指数 (公式72)
    // ========================================================================
    
    output.di_stand = calcDroughtIndex(output.stand_et_mm, output.stand_demand_mm);
    output.di_herb = calcDroughtIndex(output.et_herb_mm, output.demand_herb_mm);
    output.di_sapling = calcDroughtIndex(output.et_sapling_mm, output.demand_sapling_mm);
    output.di_tree = calcDroughtIndex(output.et_tree_mm, output.demand_tree_mm);
    
    return output;
}

} // namespace sfe
