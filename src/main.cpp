#include "block_lu_simulator.hpp"
#include "block_lu_simulator_3d.hpp"
#include "config_3d.hpp"
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include <cmath>

// Helper function to check if a number is a perfect square
bool isPerfectSquare(uint32_t n) {
    if (n == 0) return false;
    uint32_t root = static_cast<uint32_t>(std::sqrt(n));
    return root * root == n;
}

// Helper function to get integer square root
uint32_t intSqrt(uint32_t n) {
    return static_cast<uint32_t>(std::sqrt(n));
}

void printUsage(const char* program_name) {
    std::cout << "Usage: " << program_name << " [options]\n\n";
    std::cout << "Options:\n";
    std::cout << "  -n <size>      Matrix size (power of 2, default: 16)\n";
    std::cout << "  -p <size>      PE array size (power of 2, default: 4)\n";
    std::cout << "  -b <size>      Block size (power of 2, default: equals PE array size)\n";
    std::cout << "  -mac <cycles>  MAC latency in cycles (default: 1)\n";
    std::cout << "  -div <cycles>  Division latency in cycles (default: 10)\n";
    std::cout << "  -load <cycles> Memory load delay in cycles (default: 5)\n";
    std::cout << "  -write <cycles> Memory write delay in cycles (default: 5)\n";
    std::cout << "  -seed <value>  Random seed for matrix generation (default: 42)\n";
    std::cout << "  -v             Enable verbose mode (cycle-by-cycle trace)\n";
    std::cout << "  -h, --help     Show this help message\n\n";
    std::cout << "3D Options:\n";
    std::cout << "  -3d            Enable 3D simulation mode\n";
    std::cout << "  -z <layers>    Number of Z layers (default: 4)\n";
    std::cout << "  -tsv <cycles>  TSV latency per hop in cycles (default: 1)\n\n";
    std::cout << "Comparison Modes:\n";
    std::cout << "  -compare       Run both 2D and 3D with same PE array, print comparison\n";
    std::cout << "  -equal-pe      Equal PE budget comparison: 2D (p x p) vs 3D ((p/sqrt(Z)) x (p/sqrt(Z)) x Z)\n";
    std::cout << "                 Requires: Z is a perfect square and sqrt(Z) divides p\n\n";
    std::cout << "Examples:\n";
    std::cout << "  " << program_name << " -n 32 -p 8 -mac 1 -div 10 -v\n";
    std::cout << "  " << program_name << " -n 32 -p 8 -3d -z 4\n";
    std::cout << "  " << program_name << " -n 32 -p 8 -compare -z 4\n";
    std::cout << "  " << program_name << " -n 64 -p 8 -equal-pe -z 4   # 2D: 8x8=64 PEs, 3D: 4x4x4=64 PEs\n\n";
}

void printComparison(const SimStats& stats_2d, const SimStats3D& stats_3d, 
                     uint32_t num_layers) {
    std::cout << "\n";
    std::cout << "========================================================================\n";
    std::cout << "                    2D vs 3D Comparison Results                         \n";
    std::cout << "========================================================================\n";
    std::cout << "\n";
    
    // Total cycles
    double speedup = static_cast<double>(stats_2d.total_cycles) / 
                     std::max(stats_3d.total_cycles, (uint64_t)1);
    std::cout << "Total Cycles:\n";
    std::cout << "    2D:          " << std::setw(12) << stats_2d.total_cycles << "\n";
    std::cout << "    3D (" << num_layers << "L):     " << std::setw(12) << stats_3d.total_cycles << "\n";
    std::cout << "    Speedup:     " << std::fixed << std::setprecision(2) << std::setw(12) 
              << speedup << "x\n\n";
    
    // Case breakdown
    std::cout << "Case Breakdown (cycles):\n";
    std::cout << "                          2D            3D        Speedup\n";
    
    auto calc_speedup = [](uint64_t a, uint64_t b) {
        return static_cast<double>(a) / std::max(b, (uint64_t)1);
    };
    
    std::cout << "    Case 1 (diag):   " << std::setw(10) << stats_2d.case1_cycles 
              << "    " << std::setw(10) << stats_3d.case1_cycles 
              << "    " << std::setw(6) << calc_speedup(stats_2d.case1_cycles, stats_3d.case1_cycles) << "x\n";
    std::cout << "    Case 2 (horiz):  " << std::setw(10) << stats_2d.case2_cycles 
              << "    " << std::setw(10) << stats_3d.case2_cycles 
              << "    " << std::setw(6) << calc_speedup(stats_2d.case2_cycles, stats_3d.case2_cycles) << "x\n";
    std::cout << "    Case 3 (vert):   " << std::setw(10) << stats_2d.case3_cycles 
              << "    " << std::setw(10) << stats_3d.case3_cycles 
              << "    " << std::setw(6) << calc_speedup(stats_2d.case3_cycles, stats_3d.case3_cycles) << "x\n";
    std::cout << "    Case 4 (trail):  " << std::setw(10) << stats_2d.case4_cycles 
              << "    " << std::setw(10) << stats_3d.case4_cycles 
              << "    " << std::setw(6) << calc_speedup(stats_2d.case4_cycles, stats_3d.case4_cycles) << "x\n\n";
    
    // Utilization
    double util_2d = (stats_2d.total_pe_possible_cycles > 0) ?
        (100.0 * stats_2d.total_pe_active_cycles / stats_2d.total_pe_possible_cycles) : 0.0;
    double util_3d = (stats_3d.total_pe_possible_cycles > 0) ?
        (100.0 * stats_3d.total_pe_active_cycles / stats_3d.total_pe_possible_cycles) : 0.0;
    
    std::cout << "PE Utilization:\n";
    std::cout << "    2D:          " << std::setw(12) << std::setprecision(2) << util_2d << "%\n";
    std::cout << "    3D:          " << std::setw(12) << util_3d << "%\n\n";
    
    // Operations
    std::cout << "Operations:\n";
    std::cout << "    MACs:        " << std::setw(12) << stats_2d.mac_operations 
              << " (2D)  " << std::setw(12) << stats_3d.mac_operations << " (3D)\n";
    std::cout << "    DIVs:        " << std::setw(12) << stats_2d.div_operations 
              << " (2D)  " << std::setw(12) << stats_3d.div_operations << " (3D)\n\n";
    
    // TSV overhead
    std::cout << "3D-Specific Overhead:\n";
    std::cout << "    TSV Transfer Cycles: " << std::setw(10) << stats_3d.tsv_transfer_cycles << "\n";
    double tsv_overhead = (stats_3d.total_cycles > 0) ?
        (100.0 * stats_3d.tsv_transfer_cycles / stats_3d.total_cycles) : 0.0;
    std::cout << "    TSV Overhead:        " << std::setw(10) << std::setprecision(1) 
              << tsv_overhead << "%\n";
    
    std::cout << "\n========================================================================\n";
}

void printEqualPEComparison(const SimStats& stats_2d, const SimStats3D& stats_3d,
                            uint32_t pe_2d, uint32_t pe_3d, uint32_t block_2d, uint32_t block_3d,
                            uint32_t num_layers, uint32_t total_pes) {
    std::cout << "\n";
    std::cout << "========================================================================\n";
    std::cout << "           Equal PE Budget: 2D vs 3D Comparison Results                \n";
    std::cout << "========================================================================\n";
    std::cout << "\n";
    
    // PE configuration summary
    std::cout << "Silicon Budget (Total PEs): " << total_pes << "\n\n";
    std::cout << "PE Array Configuration:\n";
    std::cout << "    2D: " << pe_2d << " x " << pe_2d << " = " << (pe_2d * pe_2d) << " PEs\n";
    std::cout << "    3D: " << pe_3d << " x " << pe_3d << " x " << num_layers 
              << " = " << (pe_3d * pe_3d * num_layers) << " PEs\n\n";
    std::cout << "Block Size:\n";
    std::cout << "    2D: " << block_2d << " x " << block_2d << "\n";
    std::cout << "    3D: " << block_3d << " x " << block_3d << "\n\n";
    
    // Total cycles
    double speedup = static_cast<double>(stats_2d.total_cycles) / 
                     std::max(stats_3d.total_cycles, (uint64_t)1);
    std::cout << "Total Cycles:\n";
    std::cout << "    2D:          " << std::setw(12) << stats_2d.total_cycles << "\n";
    std::cout << "    3D (" << num_layers << "L):     " << std::setw(12) << stats_3d.total_cycles << "\n";
    std::cout << "    Speedup:     " << std::fixed << std::setprecision(2) << std::setw(12) 
              << speedup << "x";
    if (speedup > 1.0) {
        std::cout << " (3D is faster)";
    } else if (speedup < 1.0) {
        std::cout << " (2D is faster)";
    } else {
        std::cout << " (equal)";
    }
    std::cout << "\n\n";
    
    // Case breakdown
    std::cout << "Case Breakdown (cycles):\n";
    std::cout << "                          2D            3D        Speedup\n";
    
    auto calc_speedup = [](uint64_t a, uint64_t b) {
        return static_cast<double>(a) / std::max(b, (uint64_t)1);
    };
    
    std::cout << "    Case 1 (diag):   " << std::setw(10) << stats_2d.case1_cycles 
              << "    " << std::setw(10) << stats_3d.case1_cycles 
              << "    " << std::setw(6) << calc_speedup(stats_2d.case1_cycles, stats_3d.case1_cycles) << "x\n";
    std::cout << "    Case 2 (horiz):  " << std::setw(10) << stats_2d.case2_cycles 
              << "    " << std::setw(10) << stats_3d.case2_cycles 
              << "    " << std::setw(6) << calc_speedup(stats_2d.case2_cycles, stats_3d.case2_cycles) << "x\n";
    std::cout << "    Case 3 (vert):   " << std::setw(10) << stats_2d.case3_cycles 
              << "    " << std::setw(10) << stats_3d.case3_cycles 
              << "    " << std::setw(6) << calc_speedup(stats_2d.case3_cycles, stats_3d.case3_cycles) << "x\n";
    std::cout << "    Case 4 (trail):  " << std::setw(10) << stats_2d.case4_cycles 
              << "    " << std::setw(10) << stats_3d.case4_cycles 
              << "    " << std::setw(6) << calc_speedup(stats_2d.case4_cycles, stats_3d.case4_cycles) << "x\n\n";
    
    // Utilization
    double util_2d = (stats_2d.total_pe_possible_cycles > 0) ?
        (100.0 * stats_2d.total_pe_active_cycles / stats_2d.total_pe_possible_cycles) : 0.0;
    double util_3d = (stats_3d.total_pe_possible_cycles > 0) ?
        (100.0 * stats_3d.total_pe_active_cycles / stats_3d.total_pe_possible_cycles) : 0.0;
    
    std::cout << "PE Utilization:\n";
    std::cout << "    2D:          " << std::setw(12) << std::setprecision(2) << util_2d << "%\n";
    std::cout << "    3D:          " << std::setw(12) << util_3d << "%\n\n";
    
    // Operations
    std::cout << "Operations:\n";
    std::cout << "    MACs:        " << std::setw(12) << stats_2d.mac_operations 
              << " (2D)  " << std::setw(12) << stats_3d.mac_operations << " (3D)\n";
    std::cout << "    DIVs:        " << std::setw(12) << stats_2d.div_operations 
              << " (2D)  " << std::setw(12) << stats_3d.div_operations << " (3D)\n\n";
    
    // Number of blocks (affects parallelism potential)
    std::cout << "Algorithm Granularity:\n";
    std::cout << "    2D blocks: n/b = n/" << block_2d << " blocks per dimension\n";
    std::cout << "    3D blocks: n/b = n/" << block_3d << " blocks per dimension\n";
    
    // TSV overhead
    std::cout << "3D-Specific Overhead:\n";
    std::cout << "    TSV Transfer Cycles: " << std::setw(10) << stats_3d.tsv_transfer_cycles << "\n";
    double tsv_overhead = (stats_3d.total_cycles > 0) ?
        (100.0 * stats_3d.tsv_transfer_cycles / stats_3d.total_cycles) : 0.0;
    std::cout << "    TSV Overhead:        " << std::setw(10) << std::setprecision(1) 
              << tsv_overhead << "%\n";
    
    // Summary
    std::cout << "\n------------------------------------------------------------------------\n";

    std::cout << "========================================================================\n";
}

int main(int argc, char* argv[]) {
    SimConfig config;
    uint32_t seed = 42;
    
    // 3D-specific options
    bool use_3d = false;
    bool compare_mode = false;
    bool equal_pe_mode = false;
    uint32_t num_layers = 4;
    uint32_t tsv_latency = 1;
    
    // Parse command line arguments
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-n") == 0 && i + 1 < argc) {
            config.matrix_size = std::atoi(argv[++i]);
        }
        else if (strcmp(argv[i], "-p") == 0 && i + 1 < argc) {
            config.pe_array_size = std::atoi(argv[++i]);
            // Default: block_size = pe_array_size unless explicitly set
            config.block_size = config.pe_array_size;
        }
        else if (strcmp(argv[i], "-b") == 0 && i + 1 < argc) {
            config.block_size = std::atoi(argv[++i]);
        }
        else if (strcmp(argv[i], "-mac") == 0 && i + 1 < argc) {
            config.mac_latency = std::atoi(argv[++i]);
        }
        else if (strcmp(argv[i], "-div") == 0 && i + 1 < argc) {
            config.div_latency = std::atoi(argv[++i]);
        }
        else if (strcmp(argv[i], "-load") == 0 && i + 1 < argc) {
            config.mem_load_delay = std::atoi(argv[++i]);
        }
        else if (strcmp(argv[i], "-write") == 0 && i + 1 < argc) {
            config.mem_write_delay = std::atoi(argv[++i]);
        }
        else if (strcmp(argv[i], "-seed") == 0 && i + 1 < argc) {
            seed = std::atoi(argv[++i]);
        }
        else if (strcmp(argv[i], "-v") == 0) {
            config.verbose = true;
        }
        else if (strcmp(argv[i], "-3d") == 0) {
            use_3d = true;
        }
        else if (strcmp(argv[i], "-z") == 0 && i + 1 < argc) {
            num_layers = std::atoi(argv[++i]);
            use_3d = true;  // Implied
        }
        else if (strcmp(argv[i], "-tsv") == 0 && i + 1 < argc) {
            tsv_latency = std::atoi(argv[++i]);
            use_3d = true;  // Implied
        }
        else if (strcmp(argv[i], "-compare") == 0) {
            compare_mode = true;
        }
        else if (strcmp(argv[i], "-equal-pe") == 0) {
            equal_pe_mode = true;
        }
        else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            printUsage(argv[0]);
            return 0;
        }
        else {
            std::cerr << "Unknown option: " << argv[i] << "\n";
            printUsage(argv[0]);
            return 1;
        }
    }
    
    // Validate configuration
    if (!config.validate()) {
        std::cerr << "Error: Invalid configuration.\n";
        std::cerr << "Requirements:\n";
        std::cerr << "  - Matrix size must be a power of 2\n";
        std::cerr << "  - PE array size must be a power of 2\n";
        std::cerr << "  - Block size must be a power of 2\n";
        std::cerr << "  - Matrix size >= PE array size\n";
        std::cerr << "  - Block size <= PE array size\n";
        std::cerr << "  - Matrix size must be divisible by block size\n";
        return 1;
    }
    
    // Validate 3D-specific parameters
    if ((use_3d || compare_mode) && num_layers == 0) {
        std::cerr << "Error: Number of layers must be > 0\n";
        return 1;
    }
    
    // Equal PE comparison mode
    if (equal_pe_mode) {
        // Validate constraints for equal PE mode
        if (!isPerfectSquare(num_layers)) {
            std::cerr << "Error: For -equal-pe mode, Z (number of layers) must be a perfect square.\n";
            std::cerr << "  Valid values: 1, 4, 9, 16, 25, ...\n";
            std::cerr << "  Current value: Z = " << num_layers << "\n";
            return 1;
        }
        
        uint32_t sqrt_z = intSqrt(num_layers);
        if (config.pe_array_size % sqrt_z != 0) {
            std::cerr << "Error: For -equal-pe mode, sqrt(Z) must divide PE array size (p).\n";
            std::cerr << "  p = " << config.pe_array_size << ", Z = " << num_layers 
                      << ", sqrt(Z) = " << sqrt_z << "\n";
            std::cerr << "  " << config.pe_array_size << " % " << sqrt_z << " = " 
                      << (config.pe_array_size % sqrt_z) << " (must be 0)\n";
            return 1;
        }
        
        // Calculate 3D configuration
        uint32_t pe_2d = config.pe_array_size;  // 2D PE array: p x p
        uint32_t block_2d = config.block_size;  // 2D block size: b
        uint32_t pe_3d = config.pe_array_size / sqrt_z;  // 3D PE array: (p/√Z) x (p/√Z)
        uint32_t block_3d = pe_3d;  // 3D block size: same as 3D PE array size
        uint32_t total_pes = pe_2d * pe_2d;  // Total PE budget
        
        // Validate 3D block size is compatible with matrix
        if (config.matrix_size % block_3d != 0) {
            std::cerr << "Error: Matrix size must be divisible by 3D block size.\n";
            std::cerr << "  Matrix size = " << config.matrix_size << ", 3D block size = " << block_3d << "\n";
            return 1;
        }
        
        std::cout << "========================================\n";
        std::cout << "Equal PE Budget: 2D vs 3D Comparison\n";
        std::cout << "========================================\n\n";
        
        std::cout << "Configuration:\n";
        std::cout << "  Matrix size:       " << config.matrix_size << " x " << config.matrix_size << "\n";
        std::cout << "  Total PE budget:   " << total_pes << " PEs\n";
        std::cout << "  2D PE array:       " << pe_2d << " x " << pe_2d << " = " << (pe_2d * pe_2d) << " PEs\n";
        std::cout << "  2D block size:     " << block_2d << " x " << block_2d << "\n";
        std::cout << "  3D PE array:       " << pe_3d << " x " << pe_3d << " x " << num_layers 
                  << " = " << (pe_3d * pe_3d * num_layers) << " PEs\n";
        std::cout << "  3D block size:     " << block_3d << " x " << block_3d << "\n";
        std::cout << "  TSV latency:       " << tsv_latency << " cycle(s)/hop\n\n";
        
        // Run 2D simulation
        std::cout << "=== Running 2D Simulation ===\n";
        SimConfig config_2d = config;
        config_2d.pe_array_size = pe_2d;
        config_2d.block_size = block_2d;
        BlockLUSimulator sim_2d(config_2d);
        sim_2d.initializeRandom(seed);
        sim_2d.run();
        bool verified_2d = sim_2d.verify();
        std::cout << "Verification: " << (verified_2d ? "PASSED" : "FAILED") << "\n\n";
        
        // Run 3D simulation with reduced PE array
        std::cout << "=== Running 3D Simulation (" << pe_3d << "x" << pe_3d << "x" << num_layers << ") ===\n";
        SimConfig config_3d_base;
        config_3d_base.matrix_size = config.matrix_size;
        config_3d_base.pe_array_size = pe_3d;
        config_3d_base.block_size = block_3d;
        config_3d_base.mac_latency = config.mac_latency;
        config_3d_base.div_latency = config.div_latency;
        config_3d_base.mem_load_delay = config.mem_load_delay;
        config_3d_base.mem_write_delay = config.mem_write_delay;
        config_3d_base.verbose = config.verbose;
        
        SimConfig3D config_3d(config_3d_base);
        config_3d.num_layers = num_layers;
        config_3d.tsv_latency = tsv_latency;
        BlockLUSimulator3D sim_3d(config_3d);
        sim_3d.run();
        bool verified_3d = sim_3d.verify();
        std::cout << "Verification: " << (verified_3d ? "PASSED" : "FAILED") << "\n";
        
        // Print comparison
        printEqualPEComparison(sim_2d.getStats(), sim_3d.getStats(),
                               pe_2d, pe_3d, block_2d, block_3d, num_layers, total_pes);
        
        return (verified_2d && verified_3d) ? 0 : 1;
    }
    
    if (compare_mode) {
        // Run both 2D and 3D simulations for comparison
        std::cout << "========================================\n";
        std::cout << "2D vs 3D Block LU Decomposition Comparison\n";
        std::cout << "========================================\n\n";
        
        std::cout << "Configuration:\n";
        std::cout << "  Matrix size:   " << config.matrix_size << " x " << config.matrix_size << "\n";
        std::cout << "  PE array size: " << config.pe_array_size << " x " << config.pe_array_size << "\n";
        std::cout << "  Block size:    " << config.block_size << " x " << config.block_size << "\n";
        std::cout << "  Z layers:      " << num_layers << "\n";
        std::cout << "  TSV latency:   " << tsv_latency << " cycle(s)/hop\n\n";
        
        // Run 2D simulation
        std::cout << "=== Running 2D Simulation ===\n";
        BlockLUSimulator sim_2d(config);
        sim_2d.initializeRandom(seed);
        sim_2d.run();
        bool verified_2d = sim_2d.verify();
        std::cout << "Verification: " << (verified_2d ? "PASSED" : "FAILED") << "\n\n";
        
        // Run 3D simulation
        std::cout << "=== Running 3D Simulation (" << num_layers << " layers) ===\n";
        SimConfig3D config_3d(config);
        config_3d.num_layers = num_layers;
        config_3d.tsv_latency = tsv_latency;
        BlockLUSimulator3D sim_3d(config_3d);
        sim_3d.run();
        bool verified_3d = sim_3d.verify();
        std::cout << "Verification: " << (verified_3d ? "PASSED" : "FAILED") << "\n";
        
        // Print comparison
        printComparison(sim_2d.getStats(), sim_3d.getStats(), num_layers);
        
        return (verified_2d && verified_3d) ? 0 : 1;
    }
    else if (use_3d) {
        // Run 3D simulation only
        std::cout << "========================================\n";
        std::cout << "3D Block LU Decomposition Simulator\n";
        std::cout << "========================================\n\n";
        
        SimConfig3D config_3d(config);
        config_3d.num_layers = num_layers;
        config_3d.tsv_latency = tsv_latency;
        config_3d.print();
        
        BlockLUSimulator3D simulator(config_3d);
        simulator.run();
        
        // Verify result
        std::cout << "=== Verification ===\n";
        bool verified = simulator.verify();
        if (verified) {
            std::cout << "Result: PASSED\n\n";
        } else {
            std::cout << "Result: FAILED\n\n";
        }
        
        // Print results
        simulator.printResults();
        
        return verified ? 0 : 1;
    }
    else {
        // Run 2D simulation only (original behavior)
        std::cout << "========================================\n";
        std::cout << "2D Block LU Decomposition Simulator\n";
        std::cout << "========================================\n\n";
        
        BlockLUSimulator simulator(config);
        simulator.initializeRandom(seed);
        
        simulator.run();
        
        // Verify result
        std::cout << "=== Verification ===\n";
        bool verified = simulator.verify();
        if (verified) {
            std::cout << "Result: PASSED\n\n";
        } else {
            std::cout << "Result: FAILED\n\n";
        }
        
        // Print results
        simulator.printResults();
        
        return verified ? 0 : 1;
    }
}
