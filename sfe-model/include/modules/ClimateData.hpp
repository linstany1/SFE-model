/**
 * @file ClimateData.hpp
 * @brief 气候数据结构与管理器
 * 
 * 亚高山森林演替模型 (SFE-Model)
 * Phase 8: 支持多气候序列 (std::map<int, ClimateSeries>)
 */

#ifndef SFE_CLIMATE_DATA_HPP
#define SFE_CLIMATE_DATA_HPP

#include <map>
#include <vector>
#include <array>
#include <algorithm>
#include <cmath>
#include <random>
#include "core/Types.hpp"

namespace sfe {

// ============================================================================
// 气候数据结构
// ============================================================================

/**
 * @struct MonthlyClimate
 * @brief 月度气候数据
 */
struct MonthlyClimate {
    double temp = 0.0;       ///< 月均温 (°C)
    double temp_min = 0.0;   ///< 月最低温 (°C)
    double temp_max = 0.0;   ///< 月最高温 (°C)
    double prec = 0.0;       ///< 月降水量 (mm)
    double pet = 0.0;        ///< 月潜在蒸散 (mm)
    
    bool isGrowingSeason(double t_base = 5.5) const {
        return temp >= t_base;
    }
    
    double calcGDD(double t_base = 5.5) const {
        return std::max(0.0, temp - t_base) * 30.0;
    }
};

/**
 * @struct YearlyClimate
 * @brief 年度气候数据（12个月）
 */
struct YearlyClimate {
    Year year = 0;
    std::array<MonthlyClimate, 12> months;
    
    const MonthlyClimate& getMonth(int month) const {
        return months[std::clamp(month - 1, 0, 11)];
    }
    
    MonthlyClimate& getMonth(int month) {
        return months[std::clamp(month - 1, 0, 11)];
    }
    
    double calcAnnualGDD(double t_base = 5.5) const {
        double gdd = 0.0;
        for (const auto& m : months) {
            gdd += m.calcGDD(t_base);
        }
        return gdd;
    }
    
    double calcGrowingSeasonGDD(double t_base = 5.5) const {
        double gdd = 0.0;
        for (const auto& m : months) {
            if (m.isGrowingSeason(t_base)) {
                gdd += m.calcGDD(t_base);
            }
        }
        return gdd;
    }
    
    double calcAnnualPrecip() const {
        double total = 0.0;
        for (const auto& m : months) {
            total += m.prec;
        }
        return total;
    }
    
    double calcMeanTemp() const {
        double sum = 0.0;
        for (const auto& m : months) {
            sum += m.temp;
        }
        return sum / 12.0;
    }
};

/**
 * @struct LocalClimate
 * @brief 经过降尺度的本地气候
 */
struct LocalClimate {
    double temp = 0.0;
    double prec = 0.0;
    double pet = 0.0;
    double gdd_month = 0.0;
    
    LocalClimate() = default;
    LocalClimate(double t, double p, double e, double t_base = 5.5)
        : temp(t), prec(p), pet(e)
        , gdd_month(std::max(0.0, t - t_base) * 30.0)
    {}
};

// ============================================================================
// 气候降尺度工具
// ============================================================================

class ClimateDownscaler {
public:
    static constexpr double LAPSE_RATE = 0.006;
    static constexpr double REFERENCE_ELEVATION = 3500.0;
    static constexpr double PREC_CORRECTION = 0.0003;
    
    static double downscaleTemperature(double temp, double elevation, 
                                        double ref_elevation = REFERENCE_ELEVATION) {
        return temp - LAPSE_RATE * (elevation - ref_elevation);
    }
    
    static double downscalePrecipitation(double prec, double elevation,
                                          double ref_elevation = REFERENCE_ELEVATION) {
        double factor = 1.0 + PREC_CORRECTION * (elevation - ref_elevation);
        return prec * std::max(0.5, std::min(1.5, factor));
    }
    
    static double downscalePET(double pet, double elevation,
                                double ref_elevation = REFERENCE_ELEVATION) {
        double factor = 1.0 - 0.0001 * (elevation - ref_elevation);
        return pet * std::max(0.7, std::min(1.3, factor));
    }
    
    static LocalClimate downscale(const MonthlyClimate& climate, double elevation,
                                   double ref_elevation = REFERENCE_ELEVATION) {
        double t = downscaleTemperature(climate.temp, elevation, ref_elevation);
        double p = downscalePrecipitation(climate.prec, elevation, ref_elevation);
        double e = downscalePET(climate.pet, elevation, ref_elevation);
        return LocalClimate(t, p, e);
    }
};

// ============================================================================
// 气候管理器 (Phase 8: 多气候序列支持)
// ============================================================================

/**
 * @class ClimateManager
 * @brief 气候数据管理器
 * 
 * Phase 8 更新：
 * - 支持多套气候序列 (std::map<int, ClimateSeries>)
 * - 通过 climate_id 查找对应的气候数据
 * - 提供随机重采样接口用于预热阶段
 */
class ClimateManager {
public:
    /// 气候序列类型定义：年份 -> 年度气候
    using ClimateSeries = std::map<Year, YearlyClimate>;
    
    ClimateManager() = default;
    
    // ========================================================================
    // 数据加载 - 支持多气候序列
    // ========================================================================
    
    /**
     * @brief 添加年度气候数据到指定气候序列
     */
    void addYearlyClimate(int climate_id, const YearlyClimate& year_climate) {
        multi_climate_data_[climate_id][year_climate.year] = year_climate;
        if (climate_id == default_climate_id_ || multi_climate_data_.size() == 1) {
            climate_data_[year_climate.year] = year_climate;
        }
    }
    
    /**
     * @brief 添加年度气候数据（兼容旧接口）
     */
    void addYearlyClimate(const YearlyClimate& year_climate) {
        addYearlyClimate(default_climate_id_, year_climate);
    }
    
    void setDefaultClimateId(int climate_id) {
        default_climate_id_ = climate_id;
        if (multi_climate_data_.count(climate_id) > 0) {
            climate_data_ = multi_climate_data_[climate_id];
        }
    }
    
    // ========================================================================
    // 查询接口
    // ========================================================================
    
    bool hasYear(int climate_id, Year year) const {
        auto it = multi_climate_data_.find(climate_id);
        if (it == multi_climate_data_.end()) return false;
        return it->second.find(year) != it->second.end();
    }
    
    bool hasYear(Year year) const {
        return climate_data_.find(year) != climate_data_.end();
    }
    
    bool hasClimateId(int climate_id) const {
        return multi_climate_data_.find(climate_id) != multi_climate_data_.end();
    }
    
    const YearlyClimate& getYearlyClimate(int climate_id, Year year) const {
        auto cit = multi_climate_data_.find(climate_id);
        if (cit != multi_climate_data_.end()) {
            auto yit = cit->second.find(year);
            if (yit != cit->second.end()) {
                return yit->second;
            }
            return findNearestYear(cit->second, year);
        }
        return getYearlyClimate(year);
    }
    
    const YearlyClimate& getYearlyClimate(Year year) const {
        auto it = climate_data_.find(year);
        if (it == climate_data_.end()) {
            return findNearestYear(climate_data_, year);
        }
        return it->second;
    }
    
    const MonthlyClimate& getMonthlyClimate(Year year, int month) const {
        return getYearlyClimate(year).getMonth(month);
    }
    
    const MonthlyClimate& getMonthlyClimate(int climate_id, Year year, int month) const {
        return getYearlyClimate(climate_id, year).getMonth(month);
    }
    
    // ========================================================================
    // 随机重采样 - 用于预热阶段
    // ========================================================================
    
    template<typename RNG>
    const YearlyClimate& getRandomYear(int climate_id, int window_years, RNG& rng) const {
        const auto& series = getClimateSeries(climate_id);
        if (series.empty()) {
            static YearlyClimate empty;
            return empty;
        }
        
        std::vector<Year> years;
        int count = 0;
        for (const auto& [year, _] : series) {
            if (count >= window_years) break;
            years.push_back(year);
            count++;
        }
        
        if (years.empty()) {
            return series.begin()->second;
        }
        
        std::uniform_int_distribution<size_t> dist(0, years.size() - 1);
        Year selected = years[dist(rng)];
        return series.at(selected);
    }
    
    std::vector<Year> getYearList(int climate_id) const {
        std::vector<Year> years;
        auto it = multi_climate_data_.find(climate_id);
        if (it != multi_climate_data_.end()) {
            for (const auto& [year, _] : it->second) {
                years.push_back(year);
            }
        }
        return years;
    }
    
    const ClimateSeries& getClimateSeries(int climate_id) const {
        auto it = multi_climate_data_.find(climate_id);
        if (it != multi_climate_data_.end()) {
            return it->second;
        }
        return climate_data_;
    }
    
    // ========================================================================
    // 派生计算
    // ========================================================================
    
    LocalClimate getLocalClimate(Year year, int month, double elevation) const {
        const auto& raw = getMonthlyClimate(year, month);
        return ClimateDownscaler::downscale(raw, elevation);
    }
    
    double calcGDD(Year year, double elevation, bool is_evergreen) const {
        double gdd = 0.0;
        int start_month = is_evergreen ? 1 : 4;
        int end_month = is_evergreen ? 12 : 10;
        
        for (int m = start_month; m <= end_month; ++m) {
            LocalClimate local = getLocalClimate(year, m, elevation);
            gdd += local.gdd_month;
        }
        return gdd;
    }
    
    double calcNonGrowingSeasonTemp(Year year, double elevation) const {
        double sum = 0.0;
        int count = 0;
        
        if (hasYear(year - 1)) {
            for (int m = 11; m <= 12; ++m) {
                sum += getLocalClimate(year - 1, m, elevation).temp;
                ++count;
            }
        } else {
            for (int m = 11; m <= 12; ++m) {
                sum += getLocalClimate(year, m, elevation).temp;
                ++count;
            }
        }
        
        for (int m = 1; m <= 3; ++m) {
            sum += getLocalClimate(year, m, elevation).temp;
            ++count;
        }
        
        return count > 0 ? sum / count : 0.0;
    }
    
    std::pair<double, double> getNGSTemperatureRange(Year year) const {
        double min_temp_min = 1000.0;
        double max_temp_max = -1000.0;
        
        Year prev_year = year - 1;
        bool use_prev = hasYear(prev_year);
        
        if (use_prev) {
            for (int m = 11; m <= 12; ++m) {
                const auto& mc = getMonthlyClimate(prev_year, m);
                min_temp_min = std::min(min_temp_min, mc.temp_min);
                max_temp_max = std::max(max_temp_max, mc.temp_max);
            }
        } else {
            if (hasYear(year)) {
                for (int m = 11; m <= 12; ++m) {
                    const auto& mc = getMonthlyClimate(year, m);
                    min_temp_min = std::min(min_temp_min, mc.temp_min);
                    max_temp_max = std::max(max_temp_max, mc.temp_max);
                }
            }
        }
        
        if (hasYear(year)) {
            for (int m = 1; m <= 3; ++m) {
                const auto& mc = getMonthlyClimate(year, m);
                min_temp_min = std::min(min_temp_min, mc.temp_min);
                max_temp_max = std::max(max_temp_max, mc.temp_max);
            }
        }
        
        return {min_temp_min, max_temp_max};
    }
    
    bool checkNGSTemperature(Year year, double ngs_tmin_c, double ngs_tmax_c) const {
        auto [min_temp, max_temp] = getNGSTemperatureRange(year);
        return (min_temp >= ngs_tmin_c) && (max_temp < ngs_tmax_c);
    }
    
    std::pair<Year, Year> getYearRange() const {
        if (climate_data_.empty()) return {0, 0};
        return {climate_data_.begin()->first, climate_data_.rbegin()->first};
    }
    
    std::pair<Year, Year> getYearRange(int climate_id) const {
        auto it = multi_climate_data_.find(climate_id);
        if (it == multi_climate_data_.end() || it->second.empty()) {
            return {0, 0};
        }
        return {it->second.begin()->first, it->second.rbegin()->first};
    }
    
    size_t getYearCount() const {
        return climate_data_.size();
    }
    
    size_t getClimateIdCount() const {
        return multi_climate_data_.size();
    }
    
    std::vector<int> getClimateIds() const {
        std::vector<int> ids;
        for (const auto& [id, _] : multi_climate_data_) {
            ids.push_back(id);
        }
        return ids;
    }
    
    void clear() {
        climate_data_.clear();
        multi_climate_data_.clear();
    }

private:
    int default_climate_id_ = 0;
    ClimateSeries climate_data_;
    std::map<int, ClimateSeries> multi_climate_data_;
    
    static const YearlyClimate& findNearestYear(const ClimateSeries& series, Year year) {
        if (series.empty()) {
            static YearlyClimate empty;
            return empty;
        }
        
        auto it = series.lower_bound(year);
        if (it == series.end()) {
            return series.rbegin()->second;
        }
        if (it == series.begin()) {
            return it->second;
        }
        
        auto prev = std::prev(it);
        if (year - prev->first <= it->first - year) {
            return prev->second;
        }
        return it->second;
    }
};

/**
 * @struct AnnualClimateStats
 * @brief 年度气候统计
 */
struct AnnualClimateStats {
    Year year = 0;
    double gdd_annual = 0.0;
    double gdd_growing_season = 0.0;
    double ngs_temp = 0.0;
    double mean_temp = 0.0;
    double total_prec = 0.0;
    
    double getGDD(bool is_evergreen) const {
        return is_evergreen ? gdd_annual : gdd_growing_season;
    }
};

} // namespace sfe

#endif // SFE_CLIMATE_DATA_HPP
