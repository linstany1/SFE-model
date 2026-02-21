#ifndef EVENT_MANAGER_H
#define EVENT_MANAGER_H

#include <vector>
#include <map>
#include <string>

struct DisturbanceEvent {
    int year;
    int plot_id;
    std::string type;
    double intensity;
    
    DisturbanceEvent();
    DisturbanceEvent(int y, int p, const std::string& t, double i);
    
    bool isFire() const;
    bool isFullIntensity() const;
};

class EventManager {
private:
    std::map<int, std::vector<DisturbanceEvent>> events_by_year_;
    std::map<int, std::vector<DisturbanceEvent>> events_by_plot_;
    std::vector<DisturbanceEvent> all_events_;
    
public:
    EventManager();
    
    bool loadFromFile(const std::string& filename);
    bool hasEvent(int year, int plot_id) const;
    const DisturbanceEvent* getEvent(int year, int plot_id) const;
    const std::vector<DisturbanceEvent>& getEventsForYear(int year) const;
    const std::vector<DisturbanceEvent>& getEventsForPlot(int plot_id) const;
    bool hadFireAtPlot(int plot_id) const;
    std::vector<int> getEventYears() const;
    size_t getNumEvents() const;
    bool hasEvents() const;
};

#endif // EVENT_MANAGER_H
