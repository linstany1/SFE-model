#include "SpatialHash.h"
#include "AdultTree.h"
#include <algorithm>
#include <cmath>
#include <cstdio>

SpatialHash::SpatialHash() {}

int SpatialHash::coordToBucket(double coord) const {
    int bucket = static_cast<int>(coord / Config::HASH_BUCKET_SIZE);
    return std::max(0, std::min(Config::HASH_DIM - 1, bucket));
}

void SpatialHash::clear() {
    // Memory-efficient: clear contents but keep capacity
    // DO NOT call shrink_to_fit() - reduces annual reallocation overhead
    for (int i = 0; i < Config::HASH_DIM; ++i) {
        for (int j = 0; j < Config::HASH_DIM; ++j) {
            buckets_[i][j].clear();
            // Intentionally NOT calling: buckets_[i][j].shrink_to_fit();
        }
    }
}

void SpatialHash::insert(AdultTree* tree) {
    if (!tree) return;
    int bi = coordToBucket(tree->x);
    int bj = coordToBucket(tree->y);
    buckets_[bi][bj].push_back(tree);
}

void SpatialHash::remove(int tree_id) {
    for (int i = 0; i < Config::HASH_DIM; ++i) {
        for (int j = 0; j < Config::HASH_DIM; ++j) {
            auto& bucket = buckets_[i][j];
            bucket.erase(
                std::remove_if(bucket.begin(), bucket.end(),
                    [tree_id](const AdultTree* t) { 
                        return t && t->id == tree_id; 
                    }),
                bucket.end()
            );
        }
    }
}

const std::vector<AdultTree*>& SpatialHash::getBucket(int bi, int bj) const {
    static std::vector<AdultTree*> empty;
    if (bi >= 0 && bi < Config::HASH_DIM && bj >= 0 && bj < Config::HASH_DIM) {
        return buckets_[bi][bj];
    }
    return empty;
}

std::vector<AdultTree*> SpatialHash::getNearbyTrees(double x, double y) const {
    std::vector<AdultTree*> result;
    
    int bi = coordToBucket(x);
    int bj = coordToBucket(y);
    
    // [Phase 4a v1] Toroidal 3x3 neighborhood: wraps around boundaries
    for (int di = -1; di <= 1; ++di) {
        for (int dj = -1; dj <= 1; ++dj) {
            int ni = (bi + di + Config::HASH_DIM) % Config::HASH_DIM;
            int nj = (bj + dj + Config::HASH_DIM) % Config::HASH_DIM;
            
            for (AdultTree* tree : buckets_[ni][nj]) {
                result.push_back(tree);
            }
        }
    }
    
    return result;
}

std::vector<AdultTree*> SpatialHash::getTreesInRadius(double x, double y, double radius) const {
    std::vector<AdultTree*> result;
    double radius_sq = radius * radius;
    
    // [Phase 4a v1] Toroidal bucket scan: wrap indices instead of clipping
    int bi_center = coordToBucket(x);
    int bj_center = coordToBucket(y);
    
    // How many buckets does the radius span? (+1 for safety)
    int bucket_range = static_cast<int>(std::ceil(radius / Config::HASH_BUCKET_SIZE)) + 1;
    
    for (int di = -bucket_range; di <= bucket_range; ++di) {
        for (int dj = -bucket_range; dj <= bucket_range; ++dj) {
            int bi = (bi_center + di + Config::HASH_DIM * (bucket_range + 1))
                     % Config::HASH_DIM;
            int bj = (bj_center + dj + Config::HASH_DIM * (bucket_range + 1))
                     % Config::HASH_DIM;
            
            for (AdultTree* tree : buckets_[bi][bj]) {
                // Toroidal distance check
                double dx = Config::toroidalDelta(tree->x, x);
                double dy = Config::toroidalDelta(tree->y, y);
                if (dx * dx + dy * dy <= radius_sq) {
                    result.push_back(tree);
                }
            }
        }
    }
    
    return result;
}

std::vector<AdultTree*> SpatialHash::getPotentialShaders(double x, double y, 
                                                          double max_crown_radius) const {
    return getTreesInRadius(x, y, max_crown_radius + Config::HASH_BUCKET_SIZE);
}

// CRITICAL: Uses std::list to prevent pointer invalidation
void SpatialHash::rebuild(std::list<AdultTree>& trees) {
    clear();
    for (auto& tree : trees) {
        insert(&tree);
    }
}

size_t SpatialHash::getTotalCount() const {
    size_t count = 0;
    for (int i = 0; i < Config::HASH_DIM; ++i) {
        for (int j = 0; j < Config::HASH_DIM; ++j) {
            count += buckets_[i][j].size();
        }
    }
    return count;
}

void SpatialHash::printOccupancy() const {
    for (int i = 0; i < Config::HASH_DIM; ++i) {
        for (int j = 0; j < Config::HASH_DIM; ++j) {
            printf("%3zu ", buckets_[i][j].size());
        }
        printf("\n");
    }
}