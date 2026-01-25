/**
 * @file SpeciesProfile.hpp
 * @brief 物种参数结构体定义
 * 
 * 亚高山森林演替模型 (SFE-Model)
 * 
 * 对应 species_params.txt 文件中的所有列。
 * 包含 30+ 个物种参数，涵盖：
 * - 散布参数
 * - 形态参数
 * - 异速生长参数
 * - 耐受性参数
 * - 草本参数
 */

#ifndef SFE_SPECIES_PROFILE_HPP
#define SFE_SPECIES_PROFILE_HPP

#include <string>
#include <vector>
#include <map>
#include <cmath>
#include "core/Types.hpp"

namespace sfe {

/**
 * @struct SpeciesProfile
 * @brief 物种参数结构体
 * 
 * 字段严格对应 species_params.txt 中的列名
 */
struct SpeciesProfile {
    // ========================================================================
    // 基本标识
    // ========================================================================
    SpeciesId species_id = 0;          ///< 物种ID
    std::string species_name;          ///< 物种名称（可选）
    bool is_evergreen = false;         ///< 是否常绿 (0/1)
    
    // ========================================================================
    // 散布参数 (Dispersal Parameters)
    // ========================================================================
    double lambda1_m = 10.0;           ///< 短距离散布尺度参数 (m)
    double lambda2_m = 100.0;          ///< 长距离散布尺度参数 (m)
    double k_ldd = 0.1;                ///< 长距离散布比例 [0, 1]
    double fecundity_f = 1.0;          ///< 繁殖力系数
    int maturity_age_yr = 20;          ///< 性成熟年龄 (年)
    double ms_seed_kg = 0.001;         ///< 种子质量 (kg)
    
    // ========================================================================
    // 形态参数 (Morphological Parameters)
    // ========================================================================
    double h_max_m = 25.0;             ///< 最大树高 (m)
    double dbh_max_cm = 50.0;          ///< 最大胸径 (cm)
    double r_a = 0.5;                  ///< 冠幅参数 R_a
    double r_b = 0.5;                  ///< 冠幅参数 R_b
    double cs = 1.5;                   ///< 冠形系数
    double alpha_base = 0.3;           ///< 开阔地冠底高度比例
    
    // ========================================================================
    // 异速生长参数 - 叶面积 (Allometric: Leaf Area)
    // ========================================================================
    double ka1 = 0.1;                  ///< 叶面积参数 KA1
    double ka2 = 2.0;                  ///< 叶面积参数 KA2
    double kc2 = 1.0;                  ///< 叶面积参数 KC2
    
    // ========================================================================
    // 异速生长参数 - 生长 (Allometric: Growth)
    // ========================================================================
    double opt_s = 50.0;               ///< 成树异速生长参数
    double ga = 200.0;                 ///< 成树DBH生长参数
    double gs = 0.2;                   ///< 幼树高度生长参数
    
    // ========================================================================
    // 异速生长参数 - 叶生物量 (Allometric: Foliage Biomass)
    // ========================================================================
    double dw = 1.0;                   ///< 叶片干湿比
    double fw1 = 0.1;                  ///< 叶重参数 fw1
    double fw2 = 2.0;                  ///< 叶重参数 fw2
    
    // ========================================================================
    // 异速生长参数 - 幼树生物量 (Allometric: Sapling Biomass)
    // ========================================================================
    double sa = 0.1;                   ///< 幼树生物量参数 sa
    double sb = 2.5;                   ///< 幼树生物量参数 sb
    
    // ========================================================================
    // 异速生长参数 - 成树生物量 (Allometric: Adult Tree Biomass)
    // ========================================================================
    // DBH >= 5cm
    double b0 = 0.1;                   ///< 生物量参数 B0
    double b1 = 2.0;                   ///< 生物量参数 B1
    double b2 = 1.0;                   ///< 生物量参数 B2
    
    // DBH < 5cm (小树)
    double b0_small = 0.1;             ///< 小树生物量参数 b0
    double b1_small = 2.0;             ///< 小树生物量参数 b1
    double b2_small = 1.0;             ///< 小树生物量参数 b2
    
    // ========================================================================
    // 耐受性参数 (Tolerance Parameters)
    // ========================================================================
    double shade_tol = 3.0;            ///< 耐阴等级 [1, 5]
    double dr_tol = 0.5;               ///< 耐旱性 [0, 1]
    double dd_min = 300.0;             ///< 最小生长度日 (°C·day)
    
    // ========================================================================
    // 幼苗定居参数 (Seedling Establishment)
    // ========================================================================
    double l_min = 0.05;               ///< 定居最低光照需求 [0, 1]
    double l_max = 1.0;                ///< 定居最高光照需求 [0, 1]
    double ngs_tmin_c = -15.0;         ///< 非生长季最低温度限制 (°C)
    double ngs_tmax_c = 5.0;           ///< 非生长季最高温度限制 (°C)
    
    // ========================================================================
    // 草本参数 (仅 species_id = 0 有效)
    // ========================================================================
    double k_max_herb = 0.8;           ///< 草本最大容纳量 (kg/m²)
    double r_max_herb = 1.0;           ///< 草本最大生长速率
    
    // ========================================================================
    // 其他参数
    // ========================================================================
    int max_age_yr = 200;              ///< 最大年龄 (年)
    double extinction_ke = 0.5;        ///< 消光系数
    
    // ========================================================================
    // 辅助方法
    // ========================================================================
    
    /**
     * @brief 检查物种是否为草本（species_id == 0）
     */
    bool isHerb() const { return species_id == 0; }
    
    /**
     * @brief 检查树木是否达到性成熟
     */
    bool isMature(int age) const { return age >= maturity_age_yr; }
    
    /**
     * @brief 计算叶面积
     * @param dbh_cm 胸径 (cm)
     * @return 叶面积 (m²)
     */
    double calcLeafArea(double dbh_cm) const {
        return kc2 * ka1 * std::pow(dbh_cm, ka2);
    }
    
    /**
     * @brief 计算冠幅半径
     * @param dbh_cm 胸径 (cm)
     * @return 冠幅半径 (m)
     */
    double calcCrownRadius(double dbh_cm) const {
        if (dbh_cm <= 0.0) return r_b;
        return r_a * std::log(dbh_cm) + r_b;
    }
    
    /**
     * @brief 计算成树生物量
     * @param dbh_cm 胸径 (cm)
     * @param height_m 树高 (m)
     * @return 地上生物量 (kg)
     */
    double calcTreeBiomass(double dbh_cm, double height_m) const {
        if (dbh_cm >= 5.0) {
            return b0 * std::pow(dbh_cm, b1) * std::pow(height_m, b2);
        } else {
            return b0_small * std::pow(dbh_cm, b1_small) * std::pow(height_m, b2_small);
        }
    }
    
    /**
     * @brief 计算幼树生物量
     * @param height_m 高度 (m)
     * @return 地上生物量 (kg)
     */
    double calcSaplingBiomass(double height_m) const {
        return sa * std::pow(height_m, sb);
    }
    
    /**
     * @brief 计算成树叶生物量
     * @param dbh_cm 胸径 (cm)
     * @return 叶生物量 (kg)
     */
    double calcFoliageBiomass(double dbh_cm) const {
        return dw * fw1 * std::pow(dbh_cm, fw2);
    }
    
    /**
     * @brief 获取幼树叶生物量比例
     * @return 叶生物量占总生物量比例
     */
    double getSaplingFoliageFraction() const {
        return is_evergreen ? 0.4 : 0.2;
    }
};

/**
 * @class SpeciesManager
 * @brief 物种参数管理器
 * 
 * 存储和管理所有物种参数，提供按ID查询
 */
class SpeciesManager {
public:
    SpeciesManager() = default;
    
    /**
     * @brief 添加物种参数
     */
    void addSpecies(const SpeciesProfile& profile) {
        if (profile.species_id >= static_cast<int>(species_.size())) {
            species_.resize(profile.species_id + 1);
        }
        species_[profile.species_id] = profile;
        species_map_[profile.species_id] = &species_[profile.species_id];
    }
    
    /**
     * @brief 获取物种参数（按ID）
     */
    const SpeciesProfile& getSpecies(SpeciesId id) const {
        return species_.at(id);
    }
    
    /**
     * @brief 检查物种是否存在
     */
    bool hasSpecies(SpeciesId id) const {
        return id >= 0 && id < static_cast<int>(species_.size());
    }
    
    /**
     * @brief 获取物种数量
     */
    size_t getSpeciesCount() const { return species_.size(); }
    
    /**
     * @brief 获取所有物种列表
     */
    const std::vector<SpeciesProfile>& getAllSpecies() const { return species_; }
    
    /**
     * @brief 获取草本物种参数（species_id = 0）
     */
    const SpeciesProfile& getHerbSpecies() const { return species_.at(0); }
    
    /**
     * @brief 清空所有物种
     */
    void clear() { 
        species_.clear(); 
        species_map_.clear();
    }

private:
    std::vector<SpeciesProfile> species_;
    std::map<SpeciesId, SpeciesProfile*> species_map_;
};

} // namespace sfe

#endif // SFE_SPECIES_PROFILE_HPP
