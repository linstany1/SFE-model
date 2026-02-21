#include "CanopyProjectionMap.h"
#include "AdultTree.h"
#include <algorithm>
#include <cmath>

CanopyProjectionMap::CanopyProjectionMap() {}

void CanopyProjectionMap::clear() {
    // Memory-efficient: clear contents but keep capacity
    // DO NOT call shrink_to_fit() - reduces annual reallocation overhead
    for (int u = 0; u < Config::GRID_DIM; ++u) {
        for (int v = 0; v < Config::GRID_DIM; ++v) {
            projection_[u][v].clear();
            // Intentionally NOT calling: projection_[u][v].shrink_to_fit();
        }
    }
}

// [Phase 4b v1] Toroidal canopy projection: crowns near boundaries wrap around
// CRITICAL: Uses std::list to prevent pointer invalidation
void CanopyProjectionMap::build(std::list<AdultTree>& trees) {
    clear();
    
    for (auto& tree : trees) {
        // Crown radius defines the bounding box extent
        double cr = tree.crown_radius;
        
        // Scan all integer cell offsets within the crown radius + margin
        int cell_range = static_cast<int>(std::ceil(cr + 0.707));
        
        int u_center = static_cast<int>(std::floor(tree.x));
        int v_center = static_cast<int>(std::floor(tree.y));
        
        for (int du = -cell_range; du <= cell_range; ++du) {
            for (int dv = -cell_range; dv <= cell_range; ++dv) {
                // Toroidal wrap of cell index
                int u = (u_center + du % Config::GRID_DIM + Config::GRID_DIM) 
                        % Config::GRID_DIM;
                int v = (v_center + dv % Config::GRID_DIM + Config::GRID_DIM) 
                        % Config::GRID_DIM;
                
                // Cell center in toroidal space
                double cell_x = u + 0.5;
                double cell_y = v + 0.5;
                
                // Toroidal distance from tree to cell center
                double dist = Config::toroidalDist(cell_x, cell_y, tree.x, tree.y);
                
                // Include if cell center within crown radius + half-cell diagonal
                if (dist <= cr + 0.707) {
                    projection_[u][v].push_back(&tree);
                }
            }
        }
    }
}

// [Phase 4b v1] Toroidal addTree — same logic as build() for single tree
void CanopyProjectionMap::addTree(AdultTree* tree) {
    if (!tree) return;
    
    double cr = tree->crown_radius;
    int cell_range = static_cast<int>(std::ceil(cr + 0.707));
    
    int u_center = static_cast<int>(std::floor(tree->x));
    int v_center = static_cast<int>(std::floor(tree->y));
    
    for (int du = -cell_range; du <= cell_range; ++du) {
        for (int dv = -cell_range; dv <= cell_range; ++dv) {
            int u = (u_center + du % Config::GRID_DIM + Config::GRID_DIM) 
                    % Config::GRID_DIM;
            int v = (v_center + dv % Config::GRID_DIM + Config::GRID_DIM) 
                    % Config::GRID_DIM;
            
            double cell_x = u + 0.5;
            double cell_y = v + 0.5;
            double dist = Config::toroidalDist(cell_x, cell_y, tree->x, tree->y);
            
            if (dist <= cr + 0.707) {
                projection_[u][v].push_back(tree);
            }
        }
    }
}

void CanopyProjectionMap::removeTree(int tree_id) {
    for (int u = 0; u < Config::GRID_DIM; ++u) {
        for (int v = 0; v < Config::GRID_DIM; ++v) {
            auto& trees = projection_[u][v];
            trees.erase(
                std::remove_if(trees.begin(), trees.end(),
                    [tree_id](const AdultTree* t) { 
                        return t && t->id == tree_id; 
                    }),
                trees.end()
            );
        }
    }
}

const std::vector<AdultTree*>& CanopyProjectionMap::getTreesAtCell(int u, int v) const {
    static std::vector<AdultTree*> empty;
    if (u >= 0 && u < Config::GRID_DIM && v >= 0 && v < Config::GRID_DIM) {
        return projection_[u][v];
    }
    return empty;
}

const std::vector<AdultTree*>& CanopyProjectionMap::getTreesAtPoint(double x, double y) const {
    int u = static_cast<int>(x);
    int v = static_cast<int>(y);
    return getTreesAtCell(u, v);
}

size_t CanopyProjectionMap::getShaderCount(int u, int v) const {
    if (u >= 0 && u < Config::GRID_DIM && v >= 0 && v < Config::GRID_DIM) {
        return projection_[u][v].size();
    }
    return 0;
}

bool CanopyProjectionMap::hasPotentialShaders(int u, int v) const {
    return getShaderCount(u, v) > 0;
}

size_t CanopyProjectionMap::getTotalEntries() const {
    size_t count = 0;
    for (int u = 0; u < Config::GRID_DIM; ++u) {
        for (int v = 0; v < Config::GRID_DIM; ++v) {
            count += projection_[u][v].size();
        }
    }
    return count;
}

size_t CanopyProjectionMap::getMaxTreesPerCell() const {
    size_t max_count = 0;
    for (int u = 0; u < Config::GRID_DIM; ++u) {
        for (int v = 0; v < Config::GRID_DIM; ++v) {
            max_count = std::max(max_count, projection_[u][v].size());
        }
    }
    return max_count;
}