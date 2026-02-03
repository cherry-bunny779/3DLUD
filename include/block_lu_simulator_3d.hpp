#ifndef BLOCK_LU_SIMULATOR_3D_HPP
#define BLOCK_LU_SIMULATOR_3D_HPP

#include "config_3d.hpp"
#include "pe_array_3d.hpp"
#include "memory.hpp"

/**
 * 3D Block LU Decomposition Simulator
 * 
 * Implements LU decomposition using a 3D systolic array architecture.
 * 
 * Two mapping schemes are available:
 * 
 * 1. ROW DISTRIBUTION (legacy, unused):
 *    - Distributes rows within a single block across Z layers
 *    - Low PE utilization when b/Z is small
 *    - Methods suffixed with _row
 * 
 * 2. BLOCK DISTRIBUTION (active):
 *    - Distributes entire blocks across Z layers
 *    - Exploits inter-block independence in Cases 2, 3, 4
 *    - Full PE array utilization per layer
 *    - Methods suffixed with _block
 * 
 * Acceleration strategy (Block Distribution):
 * - Case 1 (diagonal): Execute on Layer 0 only (pivot dependencies)
 * - Case 2 (horizontal U): Different U[k,j] blocks on different layers
 * - Case 3 (vertical L): Different L[i,k] blocks on different layers
 * - Case 4 (trailing): Different A[i,j] blocks on different layers
 */
class BlockLUSimulator3D {
public:
    // Configuration
    const SimConfig3D& config;
    
    // 3D PE Array
    PEArray3D pe_array;
    
    // Memory subsystem (same as 2D - off-chip)
    Memory memory;
    
    // Statistics
    SimStats3D stats;
    
    // Original matrix A for verification
    std::vector<float> A_original;
    
    // Constructor
    BlockLUSimulator3D(const SimConfig3D& cfg);
    
    // Run the complete LU decomposition
    void run();
    
    // Verify result (L * U should equal original A)
    bool verify(float tolerance = 1e-4f);
    
    // Get statistics
    const SimStats3D& getStats() const { return stats; }
    
    // Print final results
    void printResults();
    
private:
    // =========================================================================
    // BLOCK DISTRIBUTION SCHEME (Active)
    // =========================================================================
    
    // Execute Case 1: Diagonal block LU decomposition
    // Executes on Layer 0 only due to sequential pivot dependencies
    void executeCase1_block(uint32_t block_k);
    
    // Execute Case 2: All horizontal U blocks in parallel
    // Each layer computes a different U[k,j] block
    void executeCase2_block(uint32_t block_k, uint32_t num_blocks);
    
    // Execute Case 3: All vertical L blocks in parallel
    // Each layer computes a different L[i,k] block
    void executeCase3_block(uint32_t block_k, uint32_t num_blocks);
    
    // Execute Case 4: All trailing blocks in parallel
    // Each layer computes different A[i,j] blocks
    void executeCase4_block(uint32_t block_k, uint32_t num_blocks);
    
    // Helper: Execute single horizontal U block on a specific layer
    void executeSingleCase2OnLayer(uint32_t layer, uint32_t block_k, uint32_t block_j,
                                   const std::vector<std::vector<float>>& L_diag);
    
    // Helper: Execute single vertical L block on a specific layer
    void executeSingleCase3OnLayer(uint32_t layer, uint32_t block_k, uint32_t block_i,
                                   const std::vector<std::vector<float>>& U_diag);
    
    // Helper: Execute single trailing block update on a specific layer
    void executeSingleCase4OnLayer(uint32_t layer, uint32_t block_k, 
                                   uint32_t block_i, uint32_t block_j);
    
    // =========================================================================
    // ROW DISTRIBUTION SCHEME (Legacy - kept but unused)
    // =========================================================================
    
    // Execute Case 1: Diagonal block LU decomposition (row distribution)
    void executeCase1_row(uint32_t block_k);
    
    // Execute Case 2: Horizontal U block update (row distribution)
    void executeCase2_row(uint32_t block_k, uint32_t block_j);
    
    // Execute Case 3: Vertical L block update (row distribution)
    void executeCase3_row(uint32_t block_k, uint32_t block_i);
    
    // Execute Case 4: Trailing matrix update (row distribution)
    void executeCase4_row(uint32_t block_k, uint32_t block_i, uint32_t block_j);
    
    // =========================================================================
    // HELPER METHODS
    // =========================================================================
    
    // Broadcast U row to all layers via TSV
    uint64_t broadcastURowToAllLayers(const std::vector<float>& u_row);
    
    // Broadcast data block to all layers
    uint64_t broadcastBlockToAllLayers(const std::vector<std::vector<float>>& block);
    
    // Broadcast block to specific layer
    uint64_t broadcastBlockToLayer(uint32_t dst_layer, 
                                   const std::vector<std::vector<float>>& block);
    
    // Get row indices assigned to a specific layer (for row distribution)
    std::vector<uint32_t> getLayerRows(uint32_t layer, uint32_t total_rows);
    
    // Load block data to specific layer's PEs
    void loadBlockToLayer(uint32_t layer, const std::vector<std::vector<float>>& block,
                          uint32_t start_row, uint32_t num_rows);
    
    // Load full block to specific layer's PEs
    void loadFullBlockToLayer(uint32_t layer, const std::vector<std::vector<float>>& block);
    
    // Write results from layer's PEs to memory
    void writeBlockFromLayer(uint32_t layer, std::vector<std::vector<float>>& block,
                             uint32_t start_row, uint32_t num_rows);
    
    // Write full block from layer's PEs
    void writeFullBlockFromLayer(uint32_t layer, std::vector<std::vector<float>>& block);
    
    // Get blocks assigned to a specific layer for block distribution
    std::vector<uint32_t> getLayerBlockIndices(uint32_t layer, uint32_t total_blocks);
    
    // Get 2D block indices assigned to a specific layer for Case 4
    std::vector<std::pair<uint32_t, uint32_t>> getLayerTrailingBlocks(
        uint32_t layer, uint32_t block_k, uint32_t num_blocks);
};

#endif // BLOCK_LU_SIMULATOR_3D_HPP
