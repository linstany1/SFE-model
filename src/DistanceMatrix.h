#ifndef DISTANCE_MATRIX_H
#define DISTANCE_MATRIX_H

#include <vector>
#include <map>
#include <string>

// Manages the distance matrix between source and target plots
// -99 indicates physical isolation (no dispersal possible)
class DistanceMatrix {
private:
    std::map<int, std::map<int, double>> matrix_;
    std::vector<int> source_ids_;
    std::vector<int> target_ids_;
    static constexpr double ISOLATION_VALUE = -99.0;
    
public:
    DistanceMatrix();
    
    bool loadFromFile(const std::string& filename);
    double getDistance(int source_id, int target_id) const;
    bool areConnected(int source_id, int target_id) const;
    static bool isIsolated(double distance);
    
    const std::vector<int>& getSourceIds() const;
    const std::vector<int>& getTargetIds() const;
    std::vector<int> getConnectedSources(int target_id) const;
    bool hasSource(int source_id) const;
    size_t getNumSources() const;
    size_t getNumTargets() const;
};

#endif // DISTANCE_MATRIX_H
