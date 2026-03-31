#ifndef CONFIG_3D_HPP
#define CONFIG_3D_HPP

#include "config.hpp"
#include <cstdint>
#include <vector>

/**
 * Extended configuration for 3D Block LU Decomposition Simulator
 * Inherits from SimConfig and adds Z-dimension parameters
 */
struct SimConfig3D : public SimConfig {
    // 3D-specific parameters
    uint32_t num_layers;        // Z dimension (number of stacked layers)
    uint32_t tsv_latency;       // Latency per TSV hop (cycles)
    uint32_t pipeline_version;  // 1 = v1 (arbitrary block assignment), 2 = v2 (row-parallel)
    
    // Default constructor
    SimConfig3D()
        : SimConfig()
        , num_layers(4)         // Default: 4 layers (as in paper)
        , tsv_latency(1)        // Default: 1 cycle per hop
        , pipeline_version(2)   // Default: v2 (new row-parallel design)
    {}
    
    // Constructor from base config
    SimConfig3D(const SimConfig& base)
        : SimConfig(base)
        , num_layers(4)
        , tsv_latency(1)
        , pipeline_version(2)
    {}
    
    // Validation (extends base validation)
    bool validate() const {
        if (!SimConfig::validate()) return false;
        if (num_layers == 0) return false;
        // num_layers doesn't need to be power of 2
        // Algorithm handles uneven distribution
        return true;
    }
    
    // Get TSV latency from layer src to layer dst
    uint32_t getTSVLatency(uint32_t src_layer, uint32_t dst_layer) const {
        if (src_layer == dst_layer) return 0;
        uint32_t hops = (src_layer > dst_layer) ? 
                        (src_layer - dst_layer) : 
                        (dst_layer - src_layer);
        return hops * tsv_latency;
    }
    
    // Get rows assigned to a specific layer for row distribution
    // Returns (start_row, num_rows) for the given layer within a block
    std::pair<uint32_t, uint32_t> getLayerRowRange(uint32_t layer, uint32_t block_size_param) const {
        uint32_t rows_per_layer = block_size_param / num_layers;
        uint32_t remainder = block_size_param % num_layers;
        
        uint32_t start_row;
        uint32_t num_rows;
        
        if (layer < remainder) {
            // Layers 0 to (remainder-1) get one extra row
            num_rows = rows_per_layer + 1;
            start_row = layer * num_rows;
        } else {
            // Remaining layers get rows_per_layer rows
            num_rows = rows_per_layer;
            start_row = remainder * (rows_per_layer + 1) + (layer - remainder) * rows_per_layer;
        }
        
        return {start_row, num_rows};
    }
    
    // Check if layer distribution is balanced
    bool isBalancedDistribution() const {
        return (block_size % num_layers) == 0;
    }
    
    // Print 3D configuration
    void print() const;
};

/**
 * Extended statistics for 3D simulation
 */
struct SimStats3D : public SimStats {
    // 3D-specific statistics
    uint64_t tsv_transfer_cycles;   // Total cycles spent on TSV transfers
    uint64_t tsv_transfers;         // Number of TSV transfer operations
    
    // Per-layer statistics
    std::vector<uint64_t> layer_active_cycles;
    std::vector<uint64_t> layer_idle_cycles;
    
    SimStats3D(uint32_t num_layers = 4)
        : SimStats()
        , tsv_transfer_cycles(0)
        , tsv_transfers(0)
        , layer_active_cycles(num_layers, 0)
        , layer_idle_cycles(num_layers, 0)
    {}
    
    // Get utilization for a specific layer
    double getLayerUtilization(uint32_t layer) const {
        if (layer >= layer_active_cycles.size()) return 0.0;
        uint64_t total = layer_active_cycles[layer] + layer_idle_cycles[layer];
        if (total == 0) return 0.0;
        return static_cast<double>(layer_active_cycles[layer]) / total * 100.0;
    }
    
    // Print 3D statistics
    void print() const;
};

#endif // CONFIG_3D_HPP
