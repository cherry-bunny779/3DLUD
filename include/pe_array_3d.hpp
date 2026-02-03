#ifndef PE_ARRAY_3D_HPP
#define PE_ARRAY_3D_HPP

#include "config_3d.hpp"
#include "pe_3d.hpp"
#include <vector>

/**
 * 3D PE Array class
 * Manages a 3D grid of Processing Elements (X × Y × Z)
 */
class PEArray3D {
public:
    // Configuration
    const SimConfig3D& config;
    
    // PE grid: [layer][row][col]
    std::vector<std::vector<std::vector<ProcessingElement3D>>> pes;
    
    // Current simulation cycle
    uint64_t current_cycle;
    
    // Constructor
    PEArray3D(const SimConfig3D& cfg);
    
    // Reset all PEs
    void reset();
    
    // Get PE at position
    ProcessingElement3D& getPE(uint32_t layer, uint32_t row, uint32_t col);
    const ProcessingElement3D& getPE(uint32_t layer, uint32_t row, uint32_t col) const;
    
    // Advance all PEs by one cycle
    void tick();
    
    // Check if all PEs are idle (across all layers)
    bool allIdle() const;
    
    // Check if all PEs on a specific layer are idle
    bool layerIdle(uint32_t layer) const;
    
    // Wait until all PEs are idle, return cycles waited
    uint64_t waitUntilIdle();
    
    // Wait until a specific layer is idle
    uint64_t waitUntilLayerIdle(uint32_t layer);
    
    // Get memory load delay for a PE position
    uint32_t getMemoryLoadDelay(uint32_t layer, uint32_t pe_row, uint32_t pe_col) const;
    
    // Get memory write delay for a PE position
    uint32_t getMemoryWriteDelay(uint32_t layer, uint32_t pe_row, uint32_t pe_col) const;
    
    // TSV operations
    
    // Broadcast data from one layer to all other layers via TSV
    // Returns cycles needed for broadcast to complete
    uint64_t broadcastToAllLayers(uint32_t src_layer, uint32_t row, uint32_t col, float data);
    
    // Send data from one PE to corresponding PE on another layer
    // Returns cycles needed
    uint64_t sendToLayer(uint32_t src_layer, uint32_t dst_layer, 
                         uint32_t row, uint32_t col, float data);
    
    // Aggregate data from all layers to one layer (e.g., sum partial results)
    // Returns cycles needed
    uint64_t aggregateToLayer(uint32_t dst_layer, uint32_t row, uint32_t col);
    
    // Collect utilization statistics
    void collectStats(SimStats3D& stats) const;
    
    // Print current state (for verbose mode)
    void printState() const;
    
    // Print state of a specific layer
    void printLayerState(uint32_t layer) const;
    
    // Count active PEs across all layers
    uint32_t countActivePEs() const;
    
    // Count active PEs on a specific layer
    uint32_t countActivePEsOnLayer(uint32_t layer) const;
};

#endif // PE_ARRAY_3D_HPP
