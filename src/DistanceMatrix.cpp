#include "DistanceMatrix.h"
#include <fstream>
#include <sstream>
#include <stdexcept>

DistanceMatrix::DistanceMatrix() {}

bool DistanceMatrix::loadFromFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open distance matrix file: " + filename);
    }
    
    std::string line;
    std::getline(file, line);
    std::istringstream header_iss(line);
    std::string token;
    
    char delim = (line.find(',') != std::string::npos) ? ',' : '\t';
    
    std::getline(header_iss, token, delim);  // Skip source_id header
    
    target_ids_.clear();
    while (std::getline(header_iss, token, delim)) {
        size_t pos = token.find('_');
        if (pos != std::string::npos) {
            try {
                int target_id = std::stoi(token.substr(pos + 1));
                target_ids_.push_back(target_id);
            } catch (...) {
                try {
                    int target_id = std::stoi(token);
                    target_ids_.push_back(target_id);
                } catch (...) {
                    continue;
                }
            }
        } else {
            try {
                int target_id = std::stoi(token);
                target_ids_.push_back(target_id);
            } catch (...) {
                continue;
            }
        }
    }
    
    matrix_.clear();
    source_ids_.clear();
    
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        
        std::istringstream iss(line);
        std::vector<std::string> tokens;
        
        while (std::getline(iss, token, delim)) {
            tokens.push_back(token);
        }
        
        if (tokens.size() < 2) continue;
        
        try {
            int source_id = std::stoi(tokens[0]);
            source_ids_.push_back(source_id);
            
            for (size_t i = 1; i < tokens.size() && i - 1 < target_ids_.size(); ++i) {
                double distance = std::stod(tokens[i]);
                int target_id = target_ids_[i - 1];
                matrix_[source_id][target_id] = distance;
            }
        } catch (const std::exception&) {
            continue;
        }
    }
    
    return !matrix_.empty();
}

double DistanceMatrix::getDistance(int source_id, int target_id) const {
    auto source_it = matrix_.find(source_id);
    if (source_it == matrix_.end()) {
        return ISOLATION_VALUE;
    }
    
    auto target_it = source_it->second.find(target_id);
    if (target_it == source_it->second.end()) {
        return ISOLATION_VALUE;
    }
    
    return target_it->second;
}

bool DistanceMatrix::areConnected(int source_id, int target_id) const {
    double dist = getDistance(source_id, target_id);
    return dist >= 0;
}

bool DistanceMatrix::isIsolated(double distance) {
    return distance < 0;
}

const std::vector<int>& DistanceMatrix::getSourceIds() const {
    return source_ids_;
}

const std::vector<int>& DistanceMatrix::getTargetIds() const {
    return target_ids_;
}

std::vector<int> DistanceMatrix::getConnectedSources(int target_id) const {
    std::vector<int> result;
    for (int source_id : source_ids_) {
        if (areConnected(source_id, target_id)) {
            result.push_back(source_id);
        }
    }
    return result;
}

bool DistanceMatrix::hasSource(int source_id) const {
    return matrix_.find(source_id) != matrix_.end();
}

size_t DistanceMatrix::getNumSources() const {
    return source_ids_.size();
}

size_t DistanceMatrix::getNumTargets() const {
    return target_ids_.size();
}
