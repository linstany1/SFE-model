#include "HerbLayer.h"
#include "SpeciesParams.h"
#include "FloatUtil.h"
#include <cmath>
#include <algorithm>

HerbLayer::HerbLayer() : enabled_(true) {  // NEW: default enabled
    reset();
}

void HerbLayer::reset() {
    for (int u = 0; u < Config::GRID_DIM; ++u) {
        for (int v = 0; v < Config::GRID_DIM; ++v) {
            biomass[u][v] = Config::MIN_HERB_BIOMASS;
            al_top[u][v] = 1.0;
        }
    }
    // Note: reset() does NOT change enabled_ state
}

double HerbLayer::calcLAI(int u, int v) const {
    // LAI_herb = 5.84 × B_herb (empirical relationship based on our study)
    return 5.84 * biomass[u][v];
}

double HerbLayer::calcFoliageWeight(int u, int v) const {
    return 0.5 * biomass[u][v];
}

double HerbLayer::getAverageLAI() const {
    double total = 0;
    for (int u = 0; u < Config::GRID_DIM; ++u) {
        for (int v = 0; v < Config::GRID_DIM; ++v) {
            total += calcLAI(u, v);
        }
    }
    return total / (Config::GRID_DIM * Config::GRID_DIM);
}

double HerbLayer::getTotalBiomass() const {
    double total = 0;
    for (int u = 0; u < Config::GRID_DIM; ++u) {
        for (int v = 0; v < Config::GRID_DIM; ++v) {
            total += biomass[u][v];
        }
    }
    return total;
}

void HerbLayer::updateCellBiomass(int u, int v, double f_env, const SpeciesProfile& herb_sp) {
    // Calculate dynamic carrying capacity
    // K_pot = K_max * f_env * Area_cell (Area_cell = 1m² = 1.0)
    double K_pot = herb_sp.k_max_herb * f_env;
    
    // Calculate effective growth rate
    double r_eff = herb_sp.r_max_herb * f_env;
    
    // Current biomass
    double B = biomass[u][v];
    
    // Logistic growth: dB/dt = r * B * (1 - B/K)
    // B_t+1 = B_t + r_eff * (1 - B_t/K_pot) * B_t
    double delta_B = 0;
    if (K_pot > Config::EPSILON) {
        delta_B = r_eff * (1.0 - B / K_pot) * B;
    }
    
    // Update biomass with minimum threshold
    B += delta_B;
    B = std::max(B, Config::MIN_HERB_BIOMASS);
    B = std::max(B, 0.0); // Ensure non-negative
    
    biomass[u][v] = B;
}

void HerbLayer::updateAllCells(const SpeciesProfile& herb_sp, 
                    double f_temp, double f_drought,
                    const std::array<std::array<double, Config::GRID_DIM>, 
                                     Config::GRID_DIM>& light_grid) {
    for (int u = 0; u < Config::GRID_DIM; ++u) {
        for (int v = 0; v < Config::GRID_DIM; ++v) {
            // Calculate light response for this cell
            double f_light = calcLightResponse(light_grid[u][v], herb_sp.shade_tol);
            
            // Combined environmental factor
            double f_env = f_temp * f_drought * f_light;
            
            // Update biomass
            updateCellBiomass(u, v, f_env, herb_sp);
        }
    }
}

double HerbLayer::calcLightResponse(double AL, int shade_tol) {
    // Ensure non-negative light
    AL = std::max(0.0, AL);
    
    // f_intol (shade intolerant, ST=1): 1.60 * (1 - exp(-1.16 * (AL - 0.15)))
    double f_intol = 1.60 * (1.0 - std::exp(-1.16 * (AL - 0.15)));
    
    // f_tol (shade tolerant, ST=5): 1.0106 * (1 - exp(-4.6 * (AL - 0.03)))
    double f_tol = 1.0106 * (1.0 - std::exp(-4.6 * (AL - 0.03)));
    
    // Interpolate based on shade tolerance
    double weight = (shade_tol - 1) / 4.0;
    double f_raw = f_intol + weight * (f_tol - f_intol);
    
    // Clamp to [0, 1] only at the outer level
    return std::max(0.0, std::min(1.0, f_raw));
}

double HerbLayer::getTotalFoliageWeight() const {
    double total = 0;
    for (int u = 0; u < Config::GRID_DIM; ++u) {
        for (int v = 0; v < Config::GRID_DIM; ++v) {
            total += calcFoliageWeight(u, v);
        }
    }
    return total;
}

// NEW: Effective foliage weight for water competition
double HerbLayer::getEffectiveFoliageWeight() const {
    // When disabled, herb does not participate in water competition
    if (!enabled_) {
        return 0.0;
    }
    return getTotalFoliageWeight();
}

// MODIFIED: Light calculation for short saplings
double HerbLayer::calcLightForShortSapling(int u, int v, double sapling_height, double k_e) const {
    // NEW: When herb competition is disabled, short saplings receive
    // the same light as the top of herb layer (no attenuation)
    if (!enabled_) {
        return al_top[u][v];
    }
    
    // Original formula when enabled:
    // AL_sapling = AL_top_herb * exp(-k_e * LAI_herb * (H_herb - H_sap) / H_herb)
    double LAI_herb = calcLAI(u, v);
    double frac = (Config::HERB_HEIGHT - sapling_height) / Config::HERB_HEIGHT;
    return al_top[u][v] * std::exp(-k_e * LAI_herb * frac);
}

double HerbLayer::getBiomass(int u, int v) const {
    if (u >= 0 && u < Config::GRID_DIM && v >= 0 && v < Config::GRID_DIM) {
        return biomass[u][v];
    }
    return 0;
}

void HerbLayer::setBiomass(int u, int v, double value) {
    if (u >= 0 && u < Config::GRID_DIM && v >= 0 && v < Config::GRID_DIM) {
        biomass[u][v] = std::max(Config::MIN_HERB_BIOMASS, value);
    }
}

void HerbLayer::setTopLight(int u, int v, double value) {
    if (u >= 0 && u < Config::GRID_DIM && v >= 0 && v < Config::GRID_DIM) {
        al_top[u][v] = std::max(0.0, std::min(1.0, value));
    }
}

double HerbLayer::getTopLight(int u, int v) const {
    if (u >= 0 && u < Config::GRID_DIM && v >= 0 && v < Config::GRID_DIM) {
        return al_top[u][v];
    }
    return 1.0;
}

double HerbLayer::getMeanFloorLight() const {
    double total = 0;
    for (int u = 0; u < Config::GRID_DIM; ++u) {
        for (int v = 0; v < Config::GRID_DIM; ++v) {
            total += al_top[u][v];
        }
    }
    return total / (Config::GRID_DIM * Config::GRID_DIM);
}
