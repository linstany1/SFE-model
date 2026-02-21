/**
 * SFE-Model: Spatially Explicit Subalpine Forest-Grassland Ecotone Model
 * 
 * A hybrid individual-based / grid-based forest succession model for
 * simulating vegetation dynamics at alpine treeline ecotones.
 * 
 * Features:
 * - Individual-based simulation of adult trees and saplings
 * - Grid-based herbaceous layer (1m x 1m cells)
 * - 3D light competition with Beer-Lambert extinction
 * - Two-layer soil water balance with asymmetric competition
 * - Dual-exponential seed dispersal kernel
 * - Fire disturbance and recovery dynamics
 * - NGS temperature filtering for establishment
 * - Leap year support in GDD calculations
 * - Parameterized stress threshold per species
 * - New mortality formulas (sigmoid + stress smoothing)
 * - Optional herb competition control (--no-herb flag)
 * 
 * Standard: C++17
 */

#include <iostream>
#include <string>
#include <filesystem>
#include <cstdlib>

#include "SimulationController.h"

namespace fs = std::filesystem;

// Print usage information
void printUsage(const char* program_name) {
    std::cout << "Usage: " << program_name << " [options]\n\n";
    std::cout << "Options:\n";
    std::cout << "  -d <dir>     Data directory (default: ./data)\n";
    std::cout << "  -o <dir>     Output directory (default: ./output)\n";
    std::cout << "  -s <seed>    Random seed (default: 42)\n";
    std::cout << "  --spinup <n> Spin-up years (default: 500)\n";
    std::cout << "  --sat <n>    Saturation years (default: 50)\n";
    std::cout << "  --start <y>  Transient start year (default: 1901)\n";
    std::cout << "  --end <y>    Transient end year (default: 2020)\n";
    std::cout << "  --interval <n> Output interval years (default: 5)\n";
    std::cout << "  --spinup-output [n]  Enable spin-up output with optional interval (default: 50)\n";
    std::cout << "  --debug      Enable debug output\n";
    std::cout << "  --no-herb    Disable herb layer competition with trees/saplings\n";
    std::cout << "  -h, --help   Show this help message\n";
}

int main(int argc, char* argv[]) {
    // Default configuration
    SimulationConfig config;
    std::string data_dir = "./data";
    
    // Parse command line arguments
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        
        if (arg == "-h" || arg == "--help") {
            printUsage(argv[0]);
            return 0;
        }
        else if (arg == "-d" && i + 1 < argc) {
            data_dir = argv[++i];
        }
        else if (arg == "-o" && i + 1 < argc) {
            config.output_dir = argv[++i];
        }
        else if (arg == "-s" && i + 1 < argc) {
            config.random_seed = static_cast<unsigned int>(std::stoul(argv[++i]));
        }
        else if (arg == "--spinup" && i + 1 < argc) {
            config.spin_up_years = std::stoi(argv[++i]);
        }
        else if (arg == "--sat" && i + 1 < argc) {
            config.saturation_years = std::stoi(argv[++i]);
        }
        else if (arg == "--start" && i + 1 < argc) {
            config.transient_start_year = std::stoi(argv[++i]);
        }
        else if (arg == "--end" && i + 1 < argc) {
            config.transient_end_year = std::stoi(argv[++i]);
        }
        else if (arg == "--interval" && i + 1 < argc) {
            config.output_interval = std::stoi(argv[++i]);
        }
        else if (arg == "--debug") {
            config.debug_mode = true;
        }
        else if (arg == "--spinup-output") {
            config.enable_spinup_output = true;
            // Check if next argument is a number (optional interval)
            if (i + 1 < argc) {
                std::string next_arg = argv[i + 1];
                // Check if it's a number (not another flag)
                if (!next_arg.empty() && next_arg[0] != '-') {
                    try {
                        config.spinup_output_interval = std::stoi(next_arg);
                        ++i;  // Consume the interval argument
                    } catch (...) {
                        // Not a number, keep default interval
                    }
                }
            }
        }
        else if (arg == "--no-herb") {
            config.enable_herb_competition = false;
        }
    }
    
    // Print header
    std::cout << "========================================\n";
    std::cout << "  SFE-Model v1.0 (Final)\n";
    std::cout << "  Subalpine Forest-Grassland Ecotone Model\n";
    std::cout << "========================================\n\n";
    
    // Print configuration
    std::cout << "Configuration:\n";
    std::cout << "  Data directory: " << data_dir << "\n";
    std::cout << "  Output directory: " << config.output_dir << "\n";
    std::cout << "  Random seed: " << config.random_seed << "\n";
    std::cout << "  Spin-up years: " << config.spin_up_years << "\n";
    std::cout << "  Saturation years: " << config.saturation_years << "\n";
    std::cout << "  Transient period: " << config.transient_start_year 
              << "-" << config.transient_end_year << "\n";
    std::cout << "  Output interval: " << config.output_interval << " years\n";
    std::cout << "  Spin-up output: " << (config.enable_spinup_output ? "ON" : "OFF");
    if (config.enable_spinup_output) {
        std::cout << " (interval: " << config.spinup_output_interval << " years)";
    }
    std::cout << "\n";
    std::cout << "  Debug mode: " << (config.debug_mode ? "ON" : "OFF") << "\n";
    std::cout << "  Herb competition: " << (config.enable_herb_competition ? "ON" : "OFF") << "\n\n";
    
    // Create output directory if it doesn't exist
    try {
        fs::create_directories(config.output_dir);
    } catch (const fs::filesystem_error& e) {
        std::cerr << "Error creating output directory: " << e.what() << "\n";
        return 1;
    }
    
    // Define input file paths
    std::string species_file = data_dir + "/species_params.csv";
    std::string site_file = data_dir + "/site.csv";
    std::string climate_file = data_dir + "/climate.csv";
    std::string distance_file = data_dir + "/distance_matrix.csv";
    std::string events_file = data_dir + "/events.csv";
    std::string seed_source_file = data_dir + "/seed_source_composition.csv";
    
    // Create and initialize simulation controller
    SimulationController controller;
    
    if (!controller.initialize(config, species_file, site_file, climate_file,
                                distance_file, events_file, seed_source_file)) {
        std::cerr << "Failed to initialize simulation. Check input files.\n";
        return 1;
    }
    
    // Run simulation
    try {
        controller.run();
    } catch (const std::exception& e) {
        std::cerr << "Simulation error: " << e.what() << "\n";
        return 1;
    }
    
    std::cout << "\nOutput files written to: " << config.output_dir << "\n";
    
    return 0;
}