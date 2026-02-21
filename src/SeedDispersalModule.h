#ifndef SEED_DISPERSAL_MODULE_H
#define SEED_DISPERSAL_MODULE_H

#include <vector>
#include <list>
#include <map>
#include <set>
#include <string>
#include "Config.h"
#include "SpeciesParams.h"

class AdultTree;
class SeedBank;
class DistanceMatrix;

struct ExternalSourceInfo {
    int plot_id;
    int species_id;
    double composition_ratio;
    
    ExternalSourceInfo();
    ExternalSourceInfo(int p, int s, double r);
};

class SeedDispersalModule {
private:
    const SpeciesParamsManager& species_params_;
    const DistanceMatrix& distance_matrix_;
    std::vector<ExternalSourceInfo> external_sources_;
    
public:
    SeedDispersalModule(const SpeciesParamsManager& sp, const DistanceMatrix& dm);
    
    bool loadExternalSources(const std::string& filename);
    
    // Dynamic b parameter based on lambda2
    static double dispersalKernel(double distance, double lambda1, double lambda2);
    
    static double calcTreeSeedProduction(const AdultTree& tree, const SpeciesProfile& sp);
    static double calcPlotSourceStrength(const std::list<AdultTree>& trees);
    
    void calcWithinPlotDispersal(const std::list<AdultTree>& trees,
                                  SeedBank& seed_bank,
                                  int target_plot_id) const;
    
    void calcBetweenPlotDispersal(int source_plot_id,
                                   const std::list<AdultTree>& source_trees,
                                   int target_plot_id,
                                   SeedBank& seed_bank) const;
    
    void calcExternalSourceDispersal(int target_plot_id,
                                      SeedBank& seed_bank) const;
    
    void applySaturatedSeedRain(SeedBank& seed_bank,
                                 const std::vector<int>& species_ids) const;
    
    static double calcSeedProb(double N_total);
    
    std::vector<int> getExternalSpeciesIds() const;
    std::vector<int> getSourcesForTarget(int target_plot_id) const;
    double getTotalInternalSeeds(int species_id, const SeedBank& seed_bank) const;
    
    // Get regional species pool for a target plot
    // Based on: 1) Species with composition_ratio > 0 in connected external sources
    //           2) Species available through connected target plots
    std::vector<int> getRegionalSpeciesPool(int target_plot_id, 
                                            const std::vector<int>& target_plot_ids) const;
};

#endif // SEED_DISPERSAL_MODULE_H
