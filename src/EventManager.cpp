#include "EventManager.h"
#include <fstream>
#include <sstream>

DisturbanceEvent::DisturbanceEvent() : year(0), plot_id(0), type(""), intensity(0) {}

DisturbanceEvent::DisturbanceEvent(int y, int p, const std::string& t, double i) 
    : year(y), plot_id(p), type(t), intensity(i) {}

bool DisturbanceEvent::isFire() const { return type == "Fire"; }
bool DisturbanceEvent::isFullIntensity() const { return intensity >= 1.0; }

EventManager::EventManager() {}

bool EventManager::loadFromFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        return true;  // OK if no events file
    }
    
    std::string line;
    std::getline(file, line);  // Skip header
    
    all_events_.clear();
    events_by_year_.clear();
    events_by_plot_.clear();
    
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        
        std::istringstream iss(line);
        std::vector<std::string> tokens;
        std::string token;
        
        char delim = (line.find(',') != std::string::npos) ? ',' : '\t';
        while (std::getline(iss, token, delim)) {
            tokens.push_back(token);
        }
        
        if (tokens.size() < 4) continue;
        
        try {
            DisturbanceEvent event;
            event.year = std::stoi(tokens[0]);
            event.plot_id = std::stoi(tokens[1]);
            event.type = tokens[2];
            event.intensity = std::stod(tokens[3]);
            
            all_events_.push_back(event);
            events_by_year_[event.year].push_back(event);
            events_by_plot_[event.plot_id].push_back(event);
            
        } catch (const std::exception&) {
            continue;
        }
    }
    
    return true;
}

bool EventManager::hasEvent(int year, int plot_id) const {
    auto year_it = events_by_year_.find(year);
    if (year_it == events_by_year_.end()) return false;
    
    for (const auto& event : year_it->second) {
        if (event.plot_id == plot_id) return true;
    }
    return false;
}

const DisturbanceEvent* EventManager::getEvent(int year, int plot_id) const {
    auto year_it = events_by_year_.find(year);
    if (year_it == events_by_year_.end()) return nullptr;
    
    for (const auto& event : year_it->second) {
        if (event.plot_id == plot_id) return &event;
    }
    return nullptr;
}

const std::vector<DisturbanceEvent>& EventManager::getEventsForYear(int year) const {
    static std::vector<DisturbanceEvent> empty;
    auto it = events_by_year_.find(year);
    return (it != events_by_year_.end()) ? it->second : empty;
}

const std::vector<DisturbanceEvent>& EventManager::getEventsForPlot(int plot_id) const {
    static std::vector<DisturbanceEvent> empty;
    auto it = events_by_plot_.find(plot_id);
    return (it != events_by_plot_.end()) ? it->second : empty;
}

bool EventManager::hadFireAtPlot(int plot_id) const {
    auto it = events_by_plot_.find(plot_id);
    if (it == events_by_plot_.end()) return false;
    
    for (const auto& event : it->second) {
        if (event.isFire()) return true;
    }
    return false;
}

std::vector<int> EventManager::getEventYears() const {
    std::vector<int> years;
    for (const auto& pair : events_by_year_) {
        years.push_back(pair.first);
    }
    return years;
}

size_t EventManager::getNumEvents() const {
    return all_events_.size();
}

bool EventManager::hasEvents() const {
    return !all_events_.empty();
}
