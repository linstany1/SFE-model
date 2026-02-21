#include "SpeciesParams.h"
#include "Config.h"
#include <fstream>
#include <sstream>
#include <stdexcept>

// Constructor with default values
SpeciesProfile::SpeciesProfile() :
    species_id(0), species_name(""), is_evergreen(false),
    lambda1_m(0), lambda2_m(0), fecundity_f(0), maturity_age_yr(0),
    max_age_yr(100), h_max_m(0), dbh_max_cm(0),
    r_a(0), r_b(0), cs(1.0), alpha_base(0), SLA(0),
    s(0), g(0), gs(0), fw1(0), fw2(0), sapHD(1.0),
    deltaBD_pot(0.5), ya(0.2), adulta(-18.0), adultb(-0.035),
    adulta_low(-18.0), adultb_low(-0.035),    // [Phase 5b v1] defaults = mean curve
    adulta_high(-18.0), adultb_high(-0.035),
    B0(0), B1(0), B2(0), b0_small(0), b1_small(0), b2_small(0),
    shade_tol(3), dr_tol(0), dd_min(0), l_min(0), l_max(1.0),
    ngs_tmin_c(-100), ngs_tmax_c(100),
    stress_threshold(Config::DEFAULT_STRESS_THRESHOLD),
    k_max_herb(0), r_max_herb(0),
    extinction_ke(0.5),
    k_top_sapling(0.60), k_top_adult(0.40) {}

bool SpeciesParamsManager::loadFromFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open species params file: " + filename);
    }
    
    std::string line;
    // Skip header line
    std::getline(file, line);
    
    profiles_.clear();
    id_to_index_.clear();
    
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        
        SpeciesProfile sp;
        std::istringstream iss(line);
        std::string token;
        
        // Parse tab or comma separated values
        std::vector<std::string> tokens;
        char delim = (line.find('\t') != std::string::npos) ? '\t' : ',';
        
        while (std::getline(iss, token, delim)) {
            tokens.push_back(token);
        }
        
        // Minimum required columns
        if (tokens.size() < 30) continue;
        
        try {
            // =====================================================
            // CSV Column Order (must match species_params.csv):
            //  0: species_id
            //  1: species_name
            //  2: is_evergreen
            //  3: lambda1_m
            //  4: lambda2_m
            //  5: fecundity_f
            //  6: maturity_age_yr
            //  7: max_age_yr
            //  8: h_max_m
            //  9: dbh_max_cm
            // 10: r_a
            // 11: r_b
            // 12: cs
            // 13: alpha_base
            // 14: SLA
            // 15: s
            // 16: g
            // 17: gs
            // 18: fw1
            // 19: fw2
            // 20: sapHD
            // 21: deltaBD_pot    â† FIX: was previously read at end
            // 22: ya             â† FIX: was previously read at end
            // 23: adulta         â† FIX: was previously read at end
            // 24: adultb         â† FIX: was previously read at end
            // 25: B0
            // 26: B1
            // 27: B2
            // 28: b0_small
            // 29: b1_small
            // 30: b2_small
            // 31: shade_tol
            // 32: dr_tol
            // 33: dd_min
            // 34: l_min
            // 35: l_max
            // 36: ngs_tmin_c
            // 37: ngs_tmax_c
            // 38: stress_threshold
            // 39: k_max_herb
            // 40: r_max_herb
            // 41: extinction_ke
            // 42: k_top_sapling
            // 43: k_top_adult
            // [Phase 5b v1] New H-D envelope columns:
            // 44: adulta_low
            // 45: adultb_low
            // 46: adulta_high
            // 47: adultb_high
            // =====================================================
            
            int idx = 0;
            sp.species_id = std::stoi(tokens[idx++]);          // 0
            sp.species_name = tokens[idx++];                    // 1
            sp.is_evergreen = (std::stoi(tokens[idx++]) == 1); // 2
            sp.lambda1_m = std::stod(tokens[idx++]);           // 3
            sp.lambda2_m = std::stod(tokens[idx++]);           // 4
            sp.fecundity_f = std::stod(tokens[idx++]);         // 5
            sp.maturity_age_yr = std::stoi(tokens[idx++]);     // 6
            sp.max_age_yr = std::stoi(tokens[idx++]);          // 7
            sp.h_max_m = std::stod(tokens[idx++]);             // 8
            sp.dbh_max_cm = std::stod(tokens[idx++]);          // 9
            sp.r_a = std::stod(tokens[idx++]);                 // 10
            sp.r_b = std::stod(tokens[idx++]);                 // 11
            sp.cs = std::stod(tokens[idx++]);                  // 12
            sp.alpha_base = std::stod(tokens[idx++]);          // 13
            sp.SLA = std::stod(tokens[idx++]);                 // 14
            sp.s = std::stod(tokens[idx++]);                   // 15
            sp.g = std::stod(tokens[idx++]);                   // 16
            sp.gs = std::stod(tokens[idx++]);                  // 17 (deprecated, kept for compatibility)
            sp.fw1 = std::stod(tokens[idx++]);                 // 18
            sp.fw2 = std::stod(tokens[idx++]);                 // 19
            sp.sapHD = std::stod(tokens[idx++]);               // 20
            
            sp.deltaBD_pot = std::stod(tokens[idx++]);         // 21
            sp.ya = std::stod(tokens[idx++]);                  // 22
            sp.adulta = std::stod(tokens[idx++]);              // 23
            sp.adultb = std::stod(tokens[idx++]);              // 24
            
            sp.B0 = std::stod(tokens[idx++]);                  // 25
            sp.B1 = std::stod(tokens[idx++]);                  // 26
            sp.B2 = std::stod(tokens[idx++]);                  // 27
            sp.b0_small = std::stod(tokens[idx++]);            // 28
            sp.b1_small = std::stod(tokens[idx++]);            // 29
            sp.b2_small = std::stod(tokens[idx++]);            // 30
            sp.shade_tol = std::stoi(tokens[idx++]);           // 31
            sp.dr_tol = std::stod(tokens[idx++]);              // 32
            sp.dd_min = std::stod(tokens[idx++]);              // 33
            sp.l_min = std::stod(tokens[idx++]);               // 34
            sp.l_max = std::stod(tokens[idx++]);               // 35
            sp.ngs_tmin_c = std::stod(tokens[idx++]);          // 36
            sp.ngs_tmax_c = std::stod(tokens[idx++]);          // 37
            
            // Optional columns with safe fallbacks
            if (idx < static_cast<int>(tokens.size())) {
                sp.stress_threshold = std::stod(tokens[idx++]); // 38
            } else {
                sp.stress_threshold = Config::DEFAULT_STRESS_THRESHOLD;
            }
            
            if (idx < static_cast<int>(tokens.size())) {
                sp.k_max_herb = std::stod(tokens[idx++]);      // 39
            }
            if (idx < static_cast<int>(tokens.size())) {
                sp.r_max_herb = std::stod(tokens[idx++]);      // 40
            }
            if (idx < static_cast<int>(tokens.size())) {
                sp.extinction_ke = std::stod(tokens[idx++]);   // 41
            }
            if (idx < static_cast<int>(tokens.size())) {
                sp.k_top_sapling = std::stod(tokens[idx++]);   // 42
            }
            if (idx < static_cast<int>(tokens.size())) {
                sp.k_top_adult = std::stod(tokens[idx++]);     // 43
            }
            
            // [Phase 5b v1] H-D plasticity envelope parameters (columns 44-47)
            if (idx < static_cast<int>(tokens.size())) {
                sp.adulta_low = std::stod(tokens[idx++]);      // 44
            }
            if (idx < static_cast<int>(tokens.size())) {
                sp.adultb_low = std::stod(tokens[idx++]);      // 45
            }
            if (idx < static_cast<int>(tokens.size())) {
                sp.adulta_high = std::stod(tokens[idx++]);     // 46
            }
            if (idx < static_cast<int>(tokens.size())) {
                sp.adultb_high = std::stod(tokens[idx++]);     // 47
            }
            
            // Fallback: if envelope not provided, use mean curve for both
            if (sp.adulta_low == 0 && sp.adultb_low == 0) {
                sp.adulta_low = sp.adulta;
                sp.adultb_low = sp.adultb;
            }
            if (sp.adulta_high == 0 && sp.adultb_high == 0) {
                sp.adulta_high = sp.adulta;
                sp.adultb_high = sp.adultb;
            }
            
            // Clamp k_top values to valid range [0, 1]
            sp.k_top_sapling = std::max(0.0, std::min(1.0, sp.k_top_sapling));
            sp.k_top_adult = std::max(0.0, std::min(1.0, sp.k_top_adult));
            
            id_to_index_[sp.species_id] = profiles_.size();
            profiles_.push_back(sp);
            
        } catch (const std::exception&) {
            // Skip malformed lines
            continue;
        }
    }
    
    return !profiles_.empty();
}

const SpeciesProfile& SpeciesParamsManager::getById(int species_id) const {
    auto it = id_to_index_.find(species_id);
    if (it == id_to_index_.end()) {
        throw std::out_of_range("Species ID not found: " + std::to_string(species_id));
    }
    return profiles_[it->second];
}

const SpeciesProfile& SpeciesParamsManager::getByIndex(size_t index) const {
    if (index >= profiles_.size()) {
        throw std::out_of_range("Species index out of range");
    }
    return profiles_[index];
}

size_t SpeciesParamsManager::getNumTreeSpecies() const {
    size_t count = 0;
    for (const auto& sp : profiles_) {
        if (sp.isTree()) count++;
    }
    return count;
}

bool SpeciesParamsManager::hasSpecies(int species_id) const {
    return id_to_index_.find(species_id) != id_to_index_.end();
}

std::vector<int> SpeciesParamsManager::getTreeSpeciesIds() const {
    std::vector<int> ids;
    for (const auto& sp : profiles_) {
        if (sp.isTree()) {
            ids.push_back(sp.species_id);
        }
    }
    return ids;
}