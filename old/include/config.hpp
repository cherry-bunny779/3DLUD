#ifndef CONFIG_HPP
#define CONFIG_HPP

#include <cstdint>
#include <string>

/**
 * Configuration structure for the Block LU Decomposition Simulator
 * All dimensions must be powers of 2
 */
struct SimConfig {
    // Matrix dimensions (N x N)
    uint32_t matrix_size;
    
    // PE array dimensions (P x P), must be <= matrix_size
    uint32_t pe_array_size;
    
    // Block size for block LU decomposition
    // Default assumption: block_size = pe_array_size
    uint32_t block_size;
    
    // Latency parameters (in cycles)
    uint32_t mac_latency;       // Multiply-accumulate latency
    uint32_t div_latency;       // Division latency
    uint32_t mem_load_delay;    // Memory to PE load delay
    uint32_t mem_write_delay;   // PE to memory write delay
    
    // Output control
    bool verbose;               // Enable cycle-by-cycle trace
    
    // Default constructor with reasonable defaults
    SimConfig() 
        : matrix_size(16)
        , pe_array_size(4)
        , block_size(4)
        , mac_latency(1)
        , div_latency(10)
        , mem_load_delay(5)
        , mem_write_delay(5)
        , verbose(false)
    {}
    
    // Validation
    bool validate() const {
        // Check power of 2
        auto isPowerOf2 = [](uint32_t n) { return n > 0 && (n & (n - 1)) == 0; };
        
        if (!isPowerOf2(matrix_size)) return false;
        if (!isPowerOf2(pe_array_size)) return false;
        if (!isPowerOf2(block_size)) return false;
        if (matrix_size < pe_array_size) return false;
        if (block_size > pe_array_size) return false;
        if (matrix_size % block_size != 0) return false;
        
        return true;
    }
    
    // Get number of blocks along one dimension
    uint32_t getNumBlocks() const {
        return matrix_size / block_size;
    }
    
    // Print configuration
    void print() const;
};

/**
 * Statistics collected during simulation
 */
struct SimStats {
    // Total cycles
    uint64_t total_cycles;
    
    // Cycles breakdown by case
    uint64_t case1_cycles;  // Diagonal block LU
    uint64_t case2_cycles;  // Horizontal U update
    uint64_t case3_cycles;  // Vertical L update
    uint64_t case4_cycles;  // Trailing matrix update
    
    // Memory cycles
    uint64_t memory_load_cycles;
    uint64_t memory_write_cycles;
    
    // PE utilization
    uint64_t total_pe_active_cycles;   // Sum of active cycles across all PEs
    uint64_t total_pe_possible_cycles; // Total possible PE-cycles
    
    // Operation counts
    uint64_t mac_operations;
    uint64_t div_operations;
    
    SimStats()
        : total_cycles(0)
        , case1_cycles(0)
        , case2_cycles(0)
        , case3_cycles(0)
        , case4_cycles(0)
        , memory_load_cycles(0)
        , memory_write_cycles(0)
        , total_pe_active_cycles(0)
        , total_pe_possible_cycles(0)
        , mac_operations(0)
        , div_operations(0)
    {}
    
    double getPEUtilization() const {
        if (total_pe_possible_cycles == 0) return 0.0;
        return static_cast<double>(total_pe_active_cycles) / total_pe_possible_cycles * 100.0;
    }
    
    void print() const;
};

/**
 * PE state enumeration
 */
enum class PEState {
    IDLE,
    LOADING,      // Loading data from memory
    COMPUTING,    // Performing MAC or DIV
    WRITING,      // Writing result to memory
    WAITING       // Waiting for data dependency
};

/**
 * Convert PEState to string for verbose output
 */
inline std::string peStateToString(PEState state) {
    switch (state) {
        case PEState::IDLE:      return "IDLE";
        case PEState::LOADING:   return "LOAD";
        case PEState::COMPUTING: return "COMP";
        case PEState::WRITING:   return "WRIT";
        case PEState::WAITING:   return "WAIT";
        default:                 return "UNKN";
    }
}

#endif // CONFIG_HPP
