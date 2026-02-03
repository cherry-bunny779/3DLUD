#include "block_lu_simulator.hpp"
#include <iostream>
#include <cstdlib>
#include <cstring>

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
    std::cout << "Example:\n";
    std::cout << "  " << program_name << " -n 32 -p 8 -mac 1 -div 10 -v\n\n";
}

int main(int argc, char* argv[]) {
    SimConfig config;
    uint32_t seed = 42;
    
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
    
    // Create and run simulator
    BlockLUSimulator simulator(config);
    simulator.initializeRandom(seed);
    
    std::cout << "========================================\n";
    std::cout << "2D Block LU Decomposition Simulator\n";
    std::cout << "========================================\n\n";
    
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
