#ifndef MEMORY_HPP
#define MEMORY_HPP

#include "config.hpp"
#include <vector>
#include <cstdint>

/**
 * Memory class
 * Represents off-chip memory storing matrices A, L, U in 2D layout
 */
class Memory {
public:
    // Configuration
    const SimConfig& config;
    
    // Matrix storage (row-major order)
    std::vector<float> A;   // Input matrix (modified in place)
    std::vector<float> L;   // Lower triangular result
    std::vector<float> U;   // Upper triangular result
    
    // Dimensions
    uint32_t n;  // Matrix dimension
    
    // Constructor
    Memory(const SimConfig& cfg);
    
    // Initialize matrix A with random values
    void initializeRandom(uint32_t seed = 42);
    
    // Initialize matrix A with specific values (for testing)
    void initializeFromArray(const std::vector<float>& data);
    
    // Initialize L as identity, U as zeros
    void initializeLU();
    
    // Access functions (2D indexing)
    float getA(uint32_t row, uint32_t col) const;
    void setA(uint32_t row, uint32_t col, float value);
    
    float getL(uint32_t row, uint32_t col) const;
    void setL(uint32_t row, uint32_t col, float value);
    
    float getU(uint32_t row, uint32_t col) const;
    void setU(uint32_t row, uint32_t col, float value);
    
    // Block access functions
    // Extract a block from matrix A at block position (block_row, block_col)
    std::vector<float> getBlockA(uint32_t block_row, uint32_t block_col) const;
    void setBlockA(uint32_t block_row, uint32_t block_col, const std::vector<float>& block);
    
    std::vector<float> getBlockL(uint32_t block_row, uint32_t block_col) const;
    void setBlockL(uint32_t block_row, uint32_t block_col, const std::vector<float>& block);
    
    std::vector<float> getBlockU(uint32_t block_row, uint32_t block_col) const;
    void setBlockU(uint32_t block_row, uint32_t block_col, const std::vector<float>& block);
    
    // Print matrices (for debugging)
    void printA() const;
    void printL() const;
    void printU() const;
    
    // Verify LU decomposition: check if L*U ≈ A_original
    bool verify(const std::vector<float>& A_original, float tolerance = 1e-4) const;
    
private:
    // Convert 2D index to 1D
    uint32_t idx(uint32_t row, uint32_t col) const {
        return row * n + col;
    }
};

#endif // MEMORY_HPP
