#ifndef SPATIAL_HASH_H
#define SPATIAL_HASH_H

#include <vector>
#include <list>
#include <array>
#include "Config.h"

// Forward declaration
class AdultTree;

class SpatialHash {
private:
    // 2D array of buckets, each containing tree pointers
    std::array<std::array<std::vector<AdultTree*>, Config::HASH_DIM>, Config::HASH_DIM> buckets_;
    
    // Convert coordinate to bucket index
    int coordToBucket(double coord) const;
    
public:
    SpatialHash();
    
    // Clear all buckets (memory-efficient: keeps capacity)
    void clear();
    
    // Insert a tree into the hash
    void insert(AdultTree* tree);
    
    // Remove a tree from the hash (by ID)
    void remove(int tree_id);
    
    // Get all trees in a specific bucket
    const std::vector<AdultTree*>& getBucket(int bi, int bj) const;
    
    // Get nearby trees (in 3x3 neighborhood of buckets)
    std::vector<AdultTree*> getNearbyTrees(double x, double y) const;
    
    // Get trees within a radius
    std::vector<AdultTree*> getTreesInRadius(double x, double y, double radius) const;
    
    // Get potential shaders for a point
    std::vector<AdultTree*> getPotentialShaders(double x, double y, double max_crown_radius) const;
    
    // Rebuild the hash from a list of trees (CRITICAL: uses std::list)
    void rebuild(std::list<AdultTree>& trees);
    
    // Get total number of trees in hash
    size_t getTotalCount() const;
    
    // Debug: print bucket occupancy
    void printOccupancy() const;
};

#endif // SPATIAL_HASH_H
