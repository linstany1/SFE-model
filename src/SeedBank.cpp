#include "SeedBank.h"
#include <algorithm>
#include <cmath>

// ===========================================
// SeedBankEntry Implementation
// ===========================================

SeedBankEntry::SeedBankEntry() : fresh_seeds(0), old_seeds(0) {}

void SeedBankEntry::annualUpdate() {
    // Discard Age 1 seeds (they've reached max lifespan)
    // Promote Age 0 to Age 1
    old_seeds = fresh_seeds;
    // Reset fresh seeds for current year
    fresh_seeds = 0;
}

void SeedBankEntry::clear() {
    fresh_seeds = 0;
    old_seeds = 0;
}

double SeedBankEntry::total() const {
    return fresh_seeds + old_seeds;
}

void SeedBankEntry::addSeeds(double amount) {
    fresh_seeds += amount;
}

// ===========================================
// CellSeedBank Implementation
// ===========================================

SeedBankEntry& CellSeedBank::getEntry(int species_id) {
    return species_seeds_[species_id];
}

const SeedBankEntry* CellSeedBank::getEntryConst(int species_id) const {
    auto it = species_seeds_.find(species_id);
    if (it != species_seeds_.end()) {
        return &(it->second);
    }
    return nullptr;
}

bool CellSeedBank::hasSeeds(int species_id) const {
    auto it = species_seeds_.find(species_id);
    if (it == species_seeds_.end()) return false;
    return it->second.total() > 0;
}

double CellSeedBank::getTotalSeeds(int species_id) const {
    auto it = species_seeds_.find(species_id);
    if (it == species_seeds_.end()) return 0;
    return it->second.total();
}

void CellSeedBank::addSeeds(int species_id, double amount) {
    species_seeds_[species_id].addSeeds(amount);
}

void CellSeedBank::annualUpdate() {
    for (auto& pair : species_seeds_) {
        pair.second.annualUpdate();
    }
}

void CellSeedBank::clear() {
    for (auto& pair : species_seeds_) {
        pair.second.clear();
    }
}

std::vector<int> CellSeedBank::getSpeciesWithSeeds() const {
    std::vector<int> result;
    for (const auto& pair : species_seeds_) {
        if (pair.second.total() > 0) {
            result.push_back(pair.first);
        }
    }
    return result;
}

// ===========================================
// SeedBank Implementation
// ===========================================

CellSeedBank& SeedBank::getCell(int u, int v) {
    return cells_[u][v];
}

const CellSeedBank& SeedBank::getCellConst(int u, int v) const {
    return cells_[u][v];
}

void SeedBank::addSeeds(int u, int v, int species_id, double amount) {
    if (u >= 0 && u < Config::GRID_DIM && v >= 0 && v < Config::GRID_DIM) {
        cells_[u][v].addSeeds(species_id, amount);
    }
}

double SeedBank::getSeeds(int u, int v, int species_id) const {
    if (u >= 0 && u < Config::GRID_DIM && v >= 0 && v < Config::GRID_DIM) {
        return cells_[u][v].getTotalSeeds(species_id);
    }
    return 0;
}

double SeedBank::calcEstablishmentProb(int u, int v, int species_id) const {
    double N_total = getSeeds(u, v, species_id);
    return std::min(N_total / Config::N_UNLIMITED, 1.0);
}

void SeedBank::annualUpdate() {
    for (int u = 0; u < Config::GRID_DIM; ++u) {
        for (int v = 0; v < Config::GRID_DIM; ++v) {
            cells_[u][v].annualUpdate();
        }
    }
}

void SeedBank::clear() {
    for (int u = 0; u < Config::GRID_DIM; ++u) {
        for (int v = 0; v < Config::GRID_DIM; ++v) {
            cells_[u][v].clear();
        }
    }
}

double SeedBank::getTotalSeedsForSpecies(int species_id) const {
    double total = 0;
    for (int u = 0; u < Config::GRID_DIM; ++u) {
        for (int v = 0; v < Config::GRID_DIM; ++v) {
            total += cells_[u][v].getTotalSeeds(species_id);
        }
    }
    return total;
}

double SeedBank::getMeanSeedDensity(int species_id) const {
    return getTotalSeedsForSpecies(species_id) / 
           (Config::GRID_DIM * Config::GRID_DIM);
}
