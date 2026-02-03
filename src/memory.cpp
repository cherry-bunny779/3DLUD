#include "memory.hpp"
#include <iostream>
#include <iomanip>
#include <random>
#include <cmath>

Memory::Memory(const SimConfig& cfg)
    : config(cfg)
    , n(cfg.matrix_size)
{
    A.resize(n * n, 0.0f);
    L.resize(n * n, 0.0f);
    U.resize(n * n, 0.0f);
}

void Memory::initializeRandom(uint32_t seed) {
    std::mt19937 gen(seed);
    std::normal_distribution<float> dist(0.0f, std::sqrt(10.0f));
    
    // Generate random matrix
    for (uint32_t i = 0; i < n * n; i++) {
        A[i] = dist(gen);
    }
    
    // Ensure matrix is not singular by making it diagonally dominant
    for (uint32_t i = 0; i < n; i++) {
        float row_sum = 0.0f;
        for (uint32_t j = 0; j < n; j++) {
            if (i != j) {
                row_sum += std::abs(getA(i, j));
            }
        }
        // Make diagonal element larger than sum of other elements in row
        setA(i, i, row_sum + 1.0f + std::abs(getA(i, i)));
    }
    
    initializeLU();
}

void Memory::initializeFromArray(const std::vector<float>& data) {
    if (data.size() != n * n) {
        std::cerr << "Error: Input data size mismatch\n";
        return;
    }
    A = data;
    initializeLU();
}

void Memory::initializeLU() {
    // Initialize L as identity matrix
    for (uint32_t i = 0; i < n; i++) {
        for (uint32_t j = 0; j < n; j++) {
            setL(i, j, (i == j) ? 1.0f : 0.0f);
            setU(i, j, 0.0f);
        }
    }
}

float Memory::getA(uint32_t row, uint32_t col) const {
    return A[idx(row, col)];
}

void Memory::setA(uint32_t row, uint32_t col, float value) {
    A[idx(row, col)] = value;
}

float Memory::getL(uint32_t row, uint32_t col) const {
    return L[idx(row, col)];
}

void Memory::setL(uint32_t row, uint32_t col, float value) {
    L[idx(row, col)] = value;
}

float Memory::getU(uint32_t row, uint32_t col) const {
    return U[idx(row, col)];
}

void Memory::setU(uint32_t row, uint32_t col, float value) {
    U[idx(row, col)] = value;
}

std::vector<float> Memory::getBlockA(uint32_t block_row, uint32_t block_col) const {
    uint32_t b = config.block_size;
    std::vector<float> block(b * b);
    
    uint32_t start_row = block_row * b;
    uint32_t start_col = block_col * b;
    
    for (uint32_t i = 0; i < b; i++) {
        for (uint32_t j = 0; j < b; j++) {
            block[i * b + j] = getA(start_row + i, start_col + j);
        }
    }
    return block;
}

void Memory::setBlockA(uint32_t block_row, uint32_t block_col, const std::vector<float>& block) {
    uint32_t b = config.block_size;
    uint32_t start_row = block_row * b;
    uint32_t start_col = block_col * b;
    
    for (uint32_t i = 0; i < b; i++) {
        for (uint32_t j = 0; j < b; j++) {
            setA(start_row + i, start_col + j, block[i * b + j]);
        }
    }
}

std::vector<float> Memory::getBlockL(uint32_t block_row, uint32_t block_col) const {
    uint32_t b = config.block_size;
    std::vector<float> block(b * b);
    
    uint32_t start_row = block_row * b;
    uint32_t start_col = block_col * b;
    
    for (uint32_t i = 0; i < b; i++) {
        for (uint32_t j = 0; j < b; j++) {
            block[i * b + j] = getL(start_row + i, start_col + j);
        }
    }
    return block;
}

void Memory::setBlockL(uint32_t block_row, uint32_t block_col, const std::vector<float>& block) {
    uint32_t b = config.block_size;
    uint32_t start_row = block_row * b;
    uint32_t start_col = block_col * b;
    
    for (uint32_t i = 0; i < b; i++) {
        for (uint32_t j = 0; j < b; j++) {
            setL(start_row + i, start_col + j, block[i * b + j]);
        }
    }
}

std::vector<float> Memory::getBlockU(uint32_t block_row, uint32_t block_col) const {
    uint32_t b = config.block_size;
    std::vector<float> block(b * b);
    
    uint32_t start_row = block_row * b;
    uint32_t start_col = block_col * b;
    
    for (uint32_t i = 0; i < b; i++) {
        for (uint32_t j = 0; j < b; j++) {
            block[i * b + j] = getU(start_row + i, start_col + j);
        }
    }
    return block;
}

void Memory::setBlockU(uint32_t block_row, uint32_t block_col, const std::vector<float>& block) {
    uint32_t b = config.block_size;
    uint32_t start_row = block_row * b;
    uint32_t start_col = block_col * b;
    
    for (uint32_t i = 0; i < b; i++) {
        for (uint32_t j = 0; j < b; j++) {
            setU(start_row + i, start_col + j, block[i * b + j]);
        }
    }
}

void Memory::printA() const {
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "Matrix A:\n";
    for (uint32_t i = 0; i < n; i++) {
        for (uint32_t j = 0; j < n; j++) {
            std::cout << std::setw(10) << getA(i, j);
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}

void Memory::printL() const {
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "Matrix L:\n";
    for (uint32_t i = 0; i < n; i++) {
        for (uint32_t j = 0; j < n; j++) {
            std::cout << std::setw(10) << getL(i, j);
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}

void Memory::printU() const {
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "Matrix U:\n";
    for (uint32_t i = 0; i < n; i++) {
        for (uint32_t j = 0; j < n; j++) {
            std::cout << std::setw(10) << getU(i, j);
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}

bool Memory::verify(const std::vector<float>& A_original, float tolerance) const {
    // Compute L * U and compare with original A
    float max_error = 0.0f;
    
    for (uint32_t i = 0; i < n; i++) {
        for (uint32_t j = 0; j < n; j++) {
            float lu_ij = 0.0f;
            for (uint32_t k = 0; k < n; k++) {
                lu_ij += getL(i, k) * getU(k, j);
            }
            float error = std::abs(lu_ij - A_original[i * n + j]);
            if (error > max_error) {
                max_error = error;
            }
        }
    }
    
    std::cout << "Verification: max |L*U - A| = " << max_error << "\n";
    return max_error < tolerance;
}
