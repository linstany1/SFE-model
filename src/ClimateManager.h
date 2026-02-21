#ifndef CLIMATE_MANAGER_H
#define CLIMATE_MANAGER_H

#include <map>
#include <vector>
#include <string>
#include <array>
#include <random>
#include "Config.h"

// Monthly climate data structure
struct MonthlyClimate {
    int month;           // 1-12
    double temp_mean;    // Monthly mean temperature (Â°C)
    double temp_min;     // Monthly minimum temperature (Â°C)
    double temp_max;     // Monthly maximum temperature (Â°C)
    double precipitation;// Precipitation (mm)
    double pet;          // Potential evapotranspiration (mm)
    
    MonthlyClimate();
};

// Annual climate data (12 months)
struct YearClimate {
    int year;
    std::array<MonthlyClimate, 13> months;  // 1-based indexing (index 0 unused)
    
    YearClimate();
    
    const MonthlyClimate& operator[](int month) const;
    MonthlyClimate& operator[](int month);
    
    // GDD calculations with leap year support
    double calcGDD_Deciduous(int target_year) const;
    double calcGDD_Evergreen(int target_year) const;
    double calcGDD(bool is_evergreen, int target_year) const;
};

// Climate series for a single climate_id
class ClimateSeries {
private:
    int climate_id_;
    std::map<int, YearClimate> years_;
    std::vector<int> available_years_;
    
public:
    ClimateSeries();
    
    void setClimateId(int id);
    int getClimateId() const;
    
    void addMonthData(int year, const MonthlyClimate& data);
    const YearClimate& getYear(int year) const;
    bool hasYear(int year) const;
    const std::vector<int>& getAvailableYears() const;
    std::vector<int> getYearsInWindow(int start_year, int window_size) const;
    int getFirstYear() const;
    int getLastYear() const;
    int getRandomYear(std::mt19937& rng, int window_size) const;
};

// Main climate manager class
class ClimateManager {
private:
    std::map<int, ClimateSeries> series_;
    
    // Cache for previous year climate
    mutable std::map<std::pair<int, int>, YearClimate> prev_year_cache_;
    
public:
    bool loadFromFile(const std::string& filename);
    
    const ClimateSeries& getSeries(int climate_id) const;
    const YearClimate& getClimate(int climate_id, int year) const;
    const YearClimate& getRandomYearClimate(int climate_id, std::mt19937& rng, 
                                             int window_size = Config::CLIMATE_BASE_WINDOW) const;
    
    // Get previous year's climate (for winter temperature calculation)
    const YearClimate& getPreviousYearClimate(int climate_id, int year, std::mt19937& rng) const;
    
    bool hasClimateId(int climate_id) const;
    std::vector<int> getClimateIds() const;
    
    // Winter Temperature Calculation
    // T_winter = mean of monthly mean temps for Dec(n-1), Jan(n), Feb(n)
    static double calcWinterTemp(const YearClimate& prev_year, 
                                  const YearClimate& curr_year);
    
    // Clear cache (call at end of simulation year if needed)
    void clearCache() const;
};

#endif // CLIMATE_MANAGER_H