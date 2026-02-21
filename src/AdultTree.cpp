#include "AdultTree.h"
#include "SpeciesParams.h"
#include "FloatUtil.h"
#include <algorithm>

// Constructor
AdultTree::AdultTree() :
    id(0), species_id(0), x(0), y(0),
    dbh(0), height(Config::BREAST_HEIGHT), crown_base_height(0),
    crown_radius(0), leaf_area(0), fol_weight(0), potential_fol_weight(0),
    leaf_area_density(0),
    age(0), stress_years(0), f_env(1.0), available_light(1.0) {}

void AdultTree::updateGeometry(const SpeciesProfile& sp) {
    // NOTE: Height is set externally by GrowthModule::calcHeightFromDBH()
    // This method only updates derived geometry properties
    
    // Update crown radius using logarithmic allometry
    // R_max - DBH ; r_b r_a
    if (dbh > 0) {
        crown_radius = std::max(0.1, sp.r_a * std::pow(dbh, sp.r_b));
    } else {
        crown_radius = 0.1;
    }
    
    // ============================================================
    // [Phase 5d v2] Incremental foliage with crown condition factor
    // ============================================================
    // potential_fol_weight = allometric maximum × crown_condition
    // fol_weight = persistent state, grows by increment, clamped to potential
    //
    // Crown condition factor: actual_crown_length / open_crown_length
    // When crown is compressed by competition (HB pushed up), potential drops,
    // breaking the "hidden pump" where suppressed trees hold full allometric foliage.
    //
    // Feedback loop: low light → HB rises → crown_condition drops
    //   → potential & fol_weight drop → leaf_area drops → ECRB drops
    //   → less water competition → neighbors benefit
    if (dbh > 0) {
        double W_pot_allometric = sp.fw1 * std::pow(dbh, sp.fw2);
        
        // Crown condition: fraction of potential crown length actually realized
        double open_crown_len = height *  sp.alpha_base;
        double actual_crown_len = height - crown_base_height;
        double crown_cond = 1.0;
        if (open_crown_len > Config::EPSILON) {
            crown_cond = actual_crown_len / open_crown_len;
            crown_cond = std::max(0.05, std::min(1.0, crown_cond));
        }
        
        double W_pot_new = W_pot_allometric * crown_cond;
        
        if (potential_fol_weight <= Config::EPSILON) {
            // First initialization (e.g. initial stand loading)
            potential_fol_weight = W_pot_new;
            fol_weight = W_pot_new;
        } else {
            // Normal growth: increment from DBH growth, or clamp down if crown compressed
            double delta_W = std::max(0.0, W_pot_new - potential_fol_weight);
            potential_fol_weight = W_pot_new;
            // Key: when crown_cond drops, W_pot_new < old potential → delta_W = 0,
            // AND fol_weight is clamped to new (lower) ceiling by min()
            fol_weight = std::min(fol_weight + delta_W, potential_fol_weight);
        }
    } else {
        potential_fol_weight = 0;
        fol_weight = 0;
    }
    
    // Update leaf area
    // LA = SLA * FolW
    leaf_area = sp.SLA * fol_weight;
    
    // Update leaf area density
    updateLeafAreaDensity(sp);
}

void AdultTree::updateLeafAreaDensity(const SpeciesProfile& sp) {
    // Crown volume approximation for cone-like crown
    // V = ÃƒÆ’Ã‚ÂÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬ * R_maxÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â² * (H - HB) / (2*cs + 1)
    double crown_length = height - crown_base_height;
    if (crown_length > 0 && crown_radius > 0) {
        double crown_volume = M_PI * crown_radius * crown_radius * 
                              crown_length / (2.0 * sp.cs + 1.0);
        if (crown_volume > 0) {
            // ÃƒÆ’Ã‚ÂÃƒâ€šÃ‚Â = LA * (2*cs + 1) / (ÃƒÆ’Ã‚ÂÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬ * R_maxÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â² * (H - HB))
            leaf_area_density = leaf_area * (2.0 * sp.cs + 1.0) / 
                                (M_PI * crown_radius * crown_radius * crown_length);
        } else {
            leaf_area_density = 0;
        }
    } else {
        leaf_area_density = 0;
    }
}

// ============================================================
// NEW v2.1: Pruning foliage loss
// ============================================================
// When crown_base_height rises from HB_old to current value,
// the lost bottom volume fraction is computed from the crown
// shape integral. The retention ratio is:
//   retention = (1 - f)^(2*cs + 1)
// where f = (HB_new - HB_old) / (H - HB_old)

void AdultTree::applyPruningLoss(const SpeciesProfile& sp, double HB_old) {
    if (crown_base_height <= HB_old) return;  // No pruning occurred
    
    double crown_length_old = height - HB_old;
    if (crown_length_old <= Config::EPSILON) return;  // Degenerate crown
    
    double f = (crown_base_height - HB_old) / crown_length_old;
    f = std::min(f, 1.0);  // Safety clamp
    
    double retention = std::pow(1.0 - f, 2.0 * sp.cs + 1.0);
    
    fol_weight *= retention;
    leaf_area = sp.SLA * fol_weight;
}

void AdultTree::updateCrownBase(const SpeciesProfile& sp, double AL_top, double AL_HB) {
    // Open-grown baseline
    double HB_open = (1.0 - sp.alpha_base) * height;
    
    // Calculate target crown base from light conditions
    double z_star = crown_base_height; // Default: no change
    
    if (AL_HB < sp.l_min && sp.l_min < AL_top) {
        // Crown base needs to move up
        double ln_diff = FloatUtil::safeLog(AL_top) - FloatUtil::safeLog(AL_HB);
        if (std::abs(ln_diff) > Config::EPSILON) {
            z_star = crown_base_height + 
                     (height - crown_base_height) * 
                     (FloatUtil::safeLog(sp.l_min) - FloatUtil::safeLog(AL_HB)) / ln_diff;
        }
    }
    
    // Apply constraints:
    // 1. Crown base can only move up (pruning is irreversible)
    // 2. Must respect open-grown baseline
    // 3. Must maintain at least 10% crown length
    double new_HB = std::max({crown_base_height, HB_open, z_star});
    crown_base_height = std::min(Config::MAX_CROWN_RATIO * height, new_HB);
}

double AdultTree::calcBiomass(const SpeciesProfile& sp) const {
    if (dbh >= 5.0) {
        // Large tree: B = B0 * DBH^B1 * H^B2
        return sp.B0 * std::pow(dbh, sp.B1) * std::pow(height, sp.B2);
    } else {
        // Small tree (DBH < 5cm): B = b0 * DBH^b1 * H^b2
        return sp.b0_small * std::pow(dbh, sp.b1_small) * std::pow(height, sp.b2_small);
    }
}

double AdultTree::getCrownRadiusAtHeight(double z, double cs) const {
    if (z < crown_base_height || z > height) return 0;
    
    // R(z) = R_max * (1 - (z - HB)/(H - HB))^cs
    double crown_frac = (z - crown_base_height) / (height - crown_base_height);
    return crown_radius * std::pow(1.0 - crown_frac, cs);
}

double AdultTree::getEntryHeight(double horizontal_dist, double cs) const {
    if (horizontal_dist >= crown_radius) return -1; // Ray misses crown
    
    // Z_in = HB + (H - HB) * [1 - (d/R_max)^(1/cs)]
    double d_ratio = horizontal_dist / crown_radius;
    return crown_base_height + 
           (height - crown_base_height) * 
           (1.0 - std::pow(d_ratio, 1.0 / cs));
}

double AdultTree::getEffectiveCrownHeight() const {
    // Z_eff = H - 0.25 * (H - HB)
    return height - 0.25 * (height - crown_base_height);
}

double AdultTree::calcSeedProduction(const SpeciesProfile& sp) const {
    if (age < sp.maturity_age_yr) return 0;
    return sp.fecundity_f * leaf_area;
}

void AdultTree::updateStressCounter(double current_f_env, double stress_threshold) {
    f_env = current_f_env;
    
    // OPTIMIZED: Stress counter with recovery mechanism
    // Simulates carbon reserve buffering
    if (f_env < stress_threshold) {
        // Under stress: accumulate stress years
        stress_years++;
    } else {
        // Recovery period: gradually reduce stress years
        if (f_env >= Config::GOOD_RECOVERY_THRESHOLD) {
            // Excellent environment: fast recovery (reduce by 2)
            stress_years = std::max(0, stress_years - 2);
        } else {
            // Moderate environment: slow recovery (reduce by 1)
            stress_years = std::max(0, stress_years - 1);
        }
    }
}

void AdultTree::incrementAge() {
    age++;
}

bool AdultTree::canOverlapAt(double target_x, double target_y) const {
    double dx = target_x - x;
    double dy = target_y - y;
    return (dx * dx + dy * dy) <= (crown_radius * crown_radius);
}