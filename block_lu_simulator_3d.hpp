#ifndef BLOCK_LU_SIMULATOR_3D_HPP
#define BLOCK_LU_SIMULATOR_3D_HPP

#include "config_3d.hpp"
#include "pe_array_3d.hpp"
#include "memory.hpp"

/**
 * 3D Block LU Decomposition Simulator
 * 
 * Implements LU decomposition using a 3D systolic array architecture
 * with row distribution across Z layers for parallelization.
 * 
 * Acceleration strategy:
 * - Case 1 (diagonal): Limited parallelization due to pivot dependencies
 * - Case 2 (horizontal U): Row distribution across layers
 * - Case 3 (vertical L): Row distribution across layers (best speedup)
 * - Case 4 (trailing): Row distribution across layers
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
    // Execute Case 1: Diagonal block LU decomposition
    // Limited 3D benefit due to sequential pivot dependencies
    void executeCase1(uint32_t block_k);
    
    // Execute Case 2: Horizontal U block update
    // U[(k-1)*b:k*b, (j-1)*b:j*b] = solve using L from Case 1
    // Row distribution: different rows processed on different layers
    void executeCase2(uint32_t block_k, uint32_t block_j);
    
    // Execute Case 3: Vertical L block update
    // L[(i-1)*b:i*b, (k-1)*b:k*b] = solve using U from Case 1
    // Row distribution: different rows of L processed on different layers
    void executeCase3(uint32_t block_k, uint32_t block_i);
    
    // Execute Case 4: Trailing matrix Schur complement update
    // A = A - L * U (outer product subtraction)
    // Row distribution: different rows of A processed on different layers
    void executeCase4(uint32_t block_k, uint32_t block_i, uint32_t block_j);
    
    // Helper: Broadcast U row to all layers via TSV
    uint64_t broadcastURowToAllLayers(const std::vector<float>& u_row);
    
    // Helper: Broadcast data block to all layers
    uint64_t broadcastBlockToAllLayers(const std::vector<std::vector<float>>& block);
    
    // Helper: Get row indices assigned to a specific layer
    std::vector<uint32_t> getLayerRows(uint32_t layer, uint32_t total_rows);
    
    // Helper: Load block data to specific layer's PEs
    void loadBlockToLayer(uint32_t layer, const std::vector<std::vector<float>>& block,
                          uint32_t start_row, uint32_t num_rows);
    
    // Helper: Write results from layer's PEs to memory
    void writeBlockFromLayer(uint32_t layer, std::vector<std::vector<float>>& block,
                             uint32_t start_row, uint32_t num_rows);
};

#endif // BLOCK_LU_SIMULATOR_3D_HPP
