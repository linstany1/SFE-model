#ifndef SEED_BANK_H
#define SEED_BANK_H

#include <array>
#include <vector>
#include <map>
#include "Config.h"

// Entry for a single species at a single cell
struct SeedBankEntry {
    double fresh_seeds;      // Current year arrivals (Age 0)
    double old_seeds;        // Previous year residual (Age 1)
    
    SeedBankEntry();
    
    // Annual update: age seeds, discard old, prepare for new
    void annualUpdate();
    
    // Clear all seeds (fire event)
    void clear();
    
    // Total viable seeds
    double total() const;
    
    // Add new arrivals
    void addSeeds(double amount);
};

// Seed bank for a single cell, managing multiple species
class CellSeedBank {
private:
    std::map<int, SeedBankEntry> species_seeds_;  // species_id -> seed entry
    
public:
    CellSeedBank() = default;
    
    // Get or create entry for a species
    SeedBankEntry& getEntry(int species_id);
    
    // Get entry (const)
    const SeedBankEntry* getEntryConst(int species_id) const;
    
    // Check if species has seeds
    bool hasSeeds(int species_id) const;
    
    // Get total seeds for species
    double getTotalSeeds(int species_id) const;
    
    // Add seeds for species
    void addSeeds(int species_id, double amount);
    
    // Perform annual update for all species
    void annualUpdate();
    
    // Clear all seeds (fire)
    void clear();
    
    // Get all species with seeds
    std::vector<int> getSpeciesWithSeeds() const;
};

// Plot-level seed bank (grid of CellSeedBanks)
class SeedBank {
private:
    std::array<std::array<CellSeedBank, Config::GRID_DIM>, Config::GRID_DIM> cells_;
    
public:
    SeedBank() = default;
    
    // Access cell seed bank
    CellSeedBank& getCell(int u, int v);
    const CellSeedBank& getCellConst(int u, int v) const;
    
    // Add seeds to a specific cell
    void addSeeds(int u, int v, int species_id, double amount);
    
    // Get total seeds at cell for species
    double getSeeds(int u, int v, int species_id) const;
    
    // Get seeds at cell (alias for compatibility)
    double getSeedsAtCell(int u, int v, int species_id) const {
        return getSeeds(u, v, species_id);
    }
    
    // Calculate establishment probability based on seed density
    double calcEstablishmentProb(int u, int v, int species_id) const;
    
    // Annual update for all cells
    void annualUpdate();
    
    // Clear all seeds (fire event)
    void clear();
    
    // Get total seeds for a species across all cells
    double getTotalSeedsForSpecies(int species_id) const;
    
    // Get mean seed density for a species
    double getMeanSeedDensity(int species_id) const;
};

#endif // SEED_BANK_H
