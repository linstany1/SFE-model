#include "ClimateManager.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <cmath>

MonthlyClimate::MonthlyClimate() : month(0), temp_mean(0), temp_min(0), temp_max(0),
                                   precipitation(0), pet(0) {}

YearClimate::YearClimate() : year(0) {}

const MonthlyClimate& YearClimate::operator[](int month) const {
    return months[month];
}

MonthlyClimate& YearClimate::operator[](int month) {
    return months[month];
}

// =====================================================================
// Sine-Curve GDD Calculation (per month)
// =====================================================================
// Must match the R parameterization code (calc_gdd_sine in
// 03_gdd_and_temp_fitting.R) that was used to fit species dd_min.
//
// The previous simple formula  GDD += days * max(0, Tmean - Tbase)
// underestimates GDD by ~23% at this site because it ignores sub-monthly
// temperature variation. When Tmean < Tbase but Tmax > Tbase, the sine
// method correctly captures the partial-day GDD contribution.
//
// Reference: Cesaraccio et al. 2001, Int J Biometeorol 45:161-169
// =====================================================================
static double calcMonthlyGDD_Sine(const MonthlyClimate& mc, int days_in_month,
                                   double Tbase) {
    double Tmax = mc.temp_max;
    double Tmin = mc.temp_min;
    double Tavg = (Tmax + Tmin) / 2.0;
    double A    = (Tmax - Tmin) / 2.0;    // daily temperature amplitude
    double Dm   = static_cast<double>(days_in_month);

    // Case 1: entire month above base temperature
    if (Tmin >= Tbase) {
        return (Tavg - Tbase) * Dm;
    }

    // Case 2: entire month below base temperature
    if (Tmax < Tbase) {
        return 0.0;
    }

    // Case 3: partial — sine-curve analytical integration
    if (A > 0.0) {
        double arg = (Tbase - Tavg) / A;
        arg = std::max(-1.0, std::min(1.0, arg));
        double theta = std::asin(arg);
        double gdd = (Dm / M_PI) *
            ((Tavg - Tbase) * (M_PI / 2.0 - theta) + A * std::cos(theta));
        return std::max(0.0, gdd);
    }

    return 0.0;
}

// GDD calculation with leap year support
double YearClimate::calcGDD_Deciduous(int target_year) const {
    double gdd = 0;
    for (int m = 4; m <= 10; ++m) {
        gdd += calcMonthlyGDD_Sine(months[m],
                                    Config::getDaysInMonth(m, target_year),
                                    Config::T_BASE);
    }
    return gdd;
}

double YearClimate::calcGDD_Evergreen(int target_year) const {
    double gdd = 0;
    for (int m = 1; m <= 12; ++m) {
        gdd += calcMonthlyGDD_Sine(months[m],
                                    Config::getDaysInMonth(m, target_year),
                                    Config::T_BASE);
    }
    return gdd;
}

double YearClimate::calcGDD(bool is_evergreen, int target_year) const {
    return is_evergreen ? calcGDD_Evergreen(target_year) : calcGDD_Deciduous(target_year);
}

ClimateSeries::ClimateSeries() : climate_id_(0) {}

void ClimateSeries::setClimateId(int id) { climate_id_ = id; }
int ClimateSeries::getClimateId() const { return climate_id_; }

void ClimateSeries::addMonthData(int year, const MonthlyClimate& data) {
    if (years_.find(year) == years_.end()) {
        years_[year].year = year;
        available_years_.push_back(year);
        std::sort(available_years_.begin(), available_years_.end());
    }
    years_[year][data.month] = data;
}

const YearClimate& ClimateSeries::getYear(int year) const {
    auto it = years_.find(year);
    if (it == years_.end()) {
        throw std::out_of_range("Year not found in climate series: " + std::to_string(year));
    }
    return it->second;
}

bool ClimateSeries::hasYear(int year) const {
    return years_.find(year) != years_.end();
}

const std::vector<int>& ClimateSeries::getAvailableYears() const {
    return available_years_;
}

std::vector<int> ClimateSeries::getYearsInWindow(int start_year, int window_size) const {
    std::vector<int> result;
    int end_year = start_year + window_size - 1;
    for (int y : available_years_) {
        if (y >= start_year && y <= end_year) {
            result.push_back(y);
        }
    }
    return result;
}

int ClimateSeries::getFirstYear() const {
    return available_years_.empty() ? 0 : available_years_.front();
}

int ClimateSeries::getLastYear() const {
    return available_years_.empty() ? 0 : available_years_.back();
}

int ClimateSeries::getRandomYear(std::mt19937& rng, int window_size) const {
    std::vector<int> window_years = getYearsInWindow(getFirstYear(), window_size);
    if (window_years.empty()) {
        return getFirstYear();
    }
    std::uniform_int_distribution<size_t> dist(0, window_years.size() - 1);
    return window_years[dist(rng)];
}

bool ClimateManager::loadFromFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open climate file: " + filename);
    }
    
    std::string line;
    std::getline(file, line);  // Skip header
    
    series_.clear();
    
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        
        std::istringstream iss(line);
        std::vector<std::string> tokens;
        std::string token;
        
        char delim = (line.find(',') != std::string::npos) ? ',' : '\t';
        while (std::getline(iss, token, delim)) {
            tokens.push_back(token);
        }
        
        if (tokens.size() < 8) continue;
        
        try {
            int climate_id = std::stoi(tokens[0]);
            int year = std::stoi(tokens[1]);
            
            MonthlyClimate mc;
            mc.month = std::stoi(tokens[2]);
            mc.temp_mean = std::stod(tokens[3]);
            mc.temp_min = std::stod(tokens[4]);
            mc.temp_max = std::stod(tokens[5]);
            mc.precipitation = std::max(0.0, std::stod(tokens[6]));
            mc.pet = std::stod(tokens[7]);
            
            if (series_.find(climate_id) == series_.end()) {
                series_[climate_id].setClimateId(climate_id);
            }
            
            series_[climate_id].addMonthData(year, mc);
            
        } catch (const std::exception&) {
            continue;
        }
    }
    
    return !series_.empty();
}

const ClimateSeries& ClimateManager::getSeries(int climate_id) const {
    auto it = series_.find(climate_id);
    if (it == series_.end()) {
        throw std::out_of_range("Climate series not found: " + std::to_string(climate_id));
    }
    return it->second;
}

const YearClimate& ClimateManager::getClimate(int climate_id, int year) const {
    return getSeries(climate_id).getYear(year);
}

const YearClimate& ClimateManager::getRandomYearClimate(int climate_id, std::mt19937& rng, 
                                                         int window_size) const {
    const auto& series = getSeries(climate_id);
    int random_year = series.getRandomYear(rng, window_size);
    return series.getYear(random_year);
}

const YearClimate& ClimateManager::getPreviousYearClimate(int climate_id, int year, 
                                                           std::mt19937& rng) const {
    auto key = std::make_pair(climate_id, year);
    
    // Check cache first
    auto cache_it = prev_year_cache_.find(key);
    if (cache_it != prev_year_cache_.end()) {
        return cache_it->second;
    }
    
    const auto& series = getSeries(climate_id);
    int prev_year = year - 1;
    
    if (series.hasYear(prev_year)) {
        prev_year_cache_[key] = series.getYear(prev_year);
    } else {
        // Use random year from base window as fallback
        prev_year_cache_[key] = getRandomYearClimate(climate_id, rng, Config::CLIMATE_BASE_WINDOW);
    }
    
    return prev_year_cache_[key];
}

bool ClimateManager::hasClimateId(int climate_id) const {
    return series_.find(climate_id) != series_.end();
}

std::vector<int> ClimateManager::getClimateIds() const {
    std::vector<int> ids;
    for (const auto& pair : series_) {
        ids.push_back(pair.first);
    }
    return ids;
}

// Winter Temperature Calculation
// T_winter = mean of Dec(n-1), Jan(n), Feb(n) monthly mean temperatures
double ClimateManager::calcWinterTemp(const YearClimate& prev_year, 
                                       const YearClimate& curr_year) {
    return (prev_year[12].temp_mean + curr_year[1].temp_mean + curr_year[2].temp_mean) / 3.0;
}

void ClimateManager::clearCache() const {
    prev_year_cache_.clear();
}