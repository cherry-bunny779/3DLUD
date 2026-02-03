#ifndef BLOCK_LU_SIMULATOR_HPP
#define BLOCK_LU_SIMULATOR_HPP

#include "config.hpp"
#include "pe_array.hpp"
#include "memory.hpp"

/**
 * Block LU Decomposition Simulator
 * 
 * Implements the block LU algorithm with four cases:
 * - Case 1: Diagonal block LU decomposition (lu_nopivot_basic)
 * - Case 2: Horizontal block update (U = solve L\A for blocks right of diagonal)
 * - Case 3: Vertical block update (L = solve A/U for blocks below diagonal)
 * - Case 4: Trailing matrix update (A = A - L*U for remaining blocks)
 */
class BlockLUSimulator {
public:
    // Configuration
    SimConfig config;
    
    // Hardware components
    PEArray pe_array;
    Memory memory;
    
    // Statistics
    SimStats stats;
    
    // Original matrix for verification
    std::vector<float> A_original;
    
    // Constructor
    BlockLUSimulator(const SimConfig& cfg);
    
    // Initialize with random matrix
    void initializeRandom(uint32_t seed = 42);
    
    // Initialize with specific matrix
    void initializeFromArray(const std::vector<float>& data);
    
    // Run the full block LU decomposition simulation
    void run();
    
    // Verify result
    bool verify(float tolerance = 1e-4) const;
    
    // Get statistics
    const SimStats& getStats() const { return stats; }
    
    // Print results
    void printResults() const;
    
private:
    // Case 1: LU decomposition on diagonal block (block_k, block_k)
    // This is the sequential pivot-based decomposition within a single block
    uint64_t executeCase1(uint32_t block_k);
    
    // Case 2: Update U block at (block_k, block_j) where block_j > block_k
    // Given L from diagonal block, compute U = L^{-1} * A
    uint64_t executeCase2(uint32_t block_k, uint32_t block_j);
    
    // Case 3: Update L block at (block_i, block_k) where block_i > block_k
    // Given U from diagonal block, compute L = A * U^{-1}
    uint64_t executeCase3(uint32_t block_i, uint32_t block_k);
    
    // Case 4: Update trailing block at (block_i, block_j) where block_i > block_k, block_j > block_k
    // Compute A = A - L * U (Schur complement update)
    uint64_t executeCase4(uint32_t block_i, uint32_t block_j, uint32_t block_k);
    
    // Helper: Load a block from memory to PE array
    uint64_t loadBlockToPEs(const std::vector<float>& block);
    
    // Helper: Write results from PE array to memory
    uint64_t writePEsToBlock(std::vector<float>& block);
    
    // Verbose logging
    void log(const std::string& message) const;
    void logCycleState() const;
};

#endif // BLOCK_LU_SIMULATOR_HPP
