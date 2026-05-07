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

    // ========================================================================
    // Bandwidth contention model
    // ========================================================================
    //
    // The 2D Partitioned architecture has 4 PE regions sharing a SINGLE memory
    // bus, so when multiple regions try to load/store simultaneously they must
    // serialize. The 3D architecture gives each layer its own memory port, so
    // there is no contention between layers (only TSV overhead for U broadcast).
    //
    // We model this by multiplying Case 4's per-layer cycles by a slowdown
    // factor derived from (memory demand) / (memory supply). See the
    // getBandwidthSlowdownFactor() method below for the demand/supply formulas.
    //
    // Toggles:
    //   is_2d_partitioned          - architecture flag (set by compare-part)
    //   enable_bandwidth_model     - master switch (false => slowdown is always 1.0)
    //   memory_bandwidth_per_layer - elements/cycle per layer (informational)
    //   contention_mode            - 0=ideal, 1=realistic (default), 2=pessimistic
    bool     is_2d_partitioned;
    bool     enable_bandwidth_model;
    uint32_t memory_bandwidth_per_layer;
    uint32_t contention_mode;

    // Default constructor
    SimConfig3D()
        : SimConfig()
        , num_layers(4)         // Default: 4 layers (as in paper)
        , tsv_latency(1)        // Default: 1 cycle per hop
        , pipeline_version(2)   // Default: v2 (new row-parallel design)
        , is_2d_partitioned(false)
        , enable_bandwidth_model(false)  // Off by default for backward compatibility
        , memory_bandwidth_per_layer(0)  // 0 = use block_size as default
        , contention_mode(1)             // Realistic contention by default
    {}

    // Constructor from base config
    SimConfig3D(const SimConfig& base)
        : SimConfig(base)
        , num_layers(4)
        , tsv_latency(1)
        , pipeline_version(2)
        , is_2d_partitioned(false)
        , enable_bandwidth_model(false)
        , memory_bandwidth_per_layer(0)
        , contention_mode(1)
    {}
    
    // Validation (extends base validation)
    bool validate() const {
        if (!SimConfig::validate()) return false;
        if (num_layers == 0) return false;
        if (contention_mode > 2) return false;
        // num_layers doesn't need to be power of 2
        // Algorithm handles uneven distribution
        return true;
    }

    // ------------------------------------------------------------------------
    // Bandwidth contention slowdown factor for Case 4
    // ------------------------------------------------------------------------
    // Returns a multiplier applied to per-layer Case 4 cycle counts.
    //
    // Demand/supply analysis (b = block_size):
    //
    //   2D Partitioned: a single ~2b-wide memory edge serves all 4 regions.
    //     Per region during streaming: A_in (b) + A_out (b) = 2b.
    //     U is read once and wire-distributed (free under realistic mode).
    //     - Realistic (mode 1) demand:    R * 2b + b   /  2b
    //     - Pessimistic (mode 2) demand:  R * 3b       /  2b   (no U broadcast)
    //   where R = num_active_regions.
    //
    //   3D New V2: each layer has its own b-wide port.
    //     Per layer: A_in + A_out + (U via TSV from layer 0) = 3b on layer 0,
    //     2b on the rest. Synchronization on the U stream pegs everyone at the
    //     layer-0 ratio, so we report a flat 3.0x slowdown (still beats 2D
    //     because 4 layers run in parallel).
    //
    // Mode 0 returns 1.0 unconditionally (legacy / ideal-broadcast behavior).
    double getBandwidthSlowdownFactor(uint32_t num_active_regions) const {
        if (!enable_bandwidth_model || num_active_regions == 0) {
            return 1.0;
        }
        if (contention_mode == 0) {
            return 1.0;  // Ideal: no contention
        }

        const uint32_t b = block_size;

        if (is_2d_partitioned) {
            const uint32_t total_bandwidth = 2 * b;  // 2D shared bus
            uint32_t demand = 0;
            if (contention_mode == 1) {
                // Realistic: U broadcast free, L/A per-region serialized
                demand = num_active_regions * 2 * b + b;
            } else {
                // Pessimistic (mode 2): full per-region serialization
                demand = num_active_regions * 3 * b;
            }
            const double slowdown =
                static_cast<double>(demand) / static_cast<double>(total_bandwidth);
            return slowdown < 1.0 ? 1.0 : slowdown;
        }

        // 3D path: independent per-layer ports, but layer 0 is the U source so
        // its 3b/b = 3 ratio dominates synchronization.
        return 3.0;
    }
    // ------------------------------------------------------------------------
    
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
