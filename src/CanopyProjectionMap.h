#ifndef CANOPY_PROJECTION_MAP_H
#define CANOPY_PROJECTION_MAP_H

#include <vector>
#include <list>
#include <array>
#include "Config.h"

// Forward declaration
class AdultTree;

// Canopy projection map for accelerating near-ground light calculations
// Records which trees' crowns could potentially shade each 1m×1m cell
class CanopyProjectionMap {
private:
    // For each cell (u,v), store pointers to trees whose crown projection covers it
    std::array<std::array<std::vector<AdultTree*>, Config::GRID_DIM>, Config::GRID_DIM> projection_;
    
public:
    CanopyProjectionMap();
    
    // Clear all projections (memory-efficient: keeps capacity)
    void clear();
    
    // Build projection map from a list of trees (CRITICAL: uses std::list)
    void build(std::list<AdultTree>& trees);
    
    // Add a single tree to the projection map (for incremental updates)
    void addTree(AdultTree* tree);
    
    // Remove a tree from the projection map by ID
    void removeTree(int tree_id);
    
    // Get trees that could shade a specific cell
    const std::vector<AdultTree*>& getTreesAtCell(int u, int v) const;
    
    // Get trees that could shade a specific point (using cell lookup)
    const std::vector<AdultTree*>& getTreesAtPoint(double x, double y) const;
    
    // Get count of potential shading trees at a cell
    size_t getShaderCount(int u, int v) const;
    
    // Check if any tree could potentially shade a cell
    bool hasPotentialShaders(int u, int v) const;
    
    // Get total projection entries (for debugging)
    size_t getTotalEntries() const;
    
    // Get maximum trees per cell (for debugging)
    size_t getMaxTreesPerCell() const;
};

#endif // CANOPY_PROJECTION_MAP_H
