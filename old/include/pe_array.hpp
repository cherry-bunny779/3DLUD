#ifndef PE_ARRAY_HPP
#define PE_ARRAY_HPP

#include "config.hpp"
#include "pe.hpp"
#include <vector>
#include <functional>

/**
 * PE Array class
 * Manages a 2D grid of Processing Elements
 */
class PEArray {
public:
    // Configuration
    const SimConfig& config;
    
    // PE grid (pe_array_size x pe_array_size)
    std::vector<std::vector<ProcessingElement>> pes;
    
    // Current simulation cycle
    uint64_t current_cycle;
    
    // Constructor
    PEArray(const SimConfig& cfg);
    
    // Reset all PEs
    void reset();
    
    // Get PE at position
    ProcessingElement& getPE(uint32_t row, uint32_t col);
    const ProcessingElement& getPE(uint32_t row, uint32_t col) const;
    
    // Advance all PEs by one cycle
    void tick();
    
    // Check if all PEs are idle
    bool allIdle() const;
    
    // Check if all PEs have completed current operation
    bool allComplete() const;
    
    // Wait until all PEs are idle, return cycles waited
    uint64_t waitUntilIdle();
    
    // Get memory load delay for a PE position (dummy function for now)
    uint32_t getMemoryLoadDelay(uint32_t pe_row, uint32_t pe_col) const;
    
    // Get memory write delay for a PE position (dummy function for now)
    uint32_t getMemoryWriteDelay(uint32_t pe_row, uint32_t pe_col) const;
    
    // Collect utilization statistics
    void collectStats(SimStats& stats) const;
    
    // Print current state (for verbose mode)
    void printState() const;
    
    // Count active PEs
    uint32_t countActivePEs() const;
};

#endif // PE_ARRAY_HPP
