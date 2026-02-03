#include "block_lu_simulator_3d.hpp"
#include <iostream>
#include <cmath>
#include <algorithm>

// Helper: Convert flat vector to 2D block
static std::vector<std::vector<float>> flatTo2D(const std::vector<float>& flat, uint32_t b) {
    std::vector<std::vector<float>> block(b, std::vector<float>(b, 0.0f));
    for (uint32_t i = 0; i < b; i++) {
        for (uint32_t j = 0; j < b; j++) {
            block[i][j] = flat[i * b + j];
        }
    }
    return block;
}

// Helper: Convert 2D block to flat vector
static std::vector<float> toFlat(const std::vector<std::vector<float>>& block, uint32_t b) {
    std::vector<float> flat(b * b, 0.0f);
    for (uint32_t i = 0; i < b; i++) {
        for (uint32_t j = 0; j < b; j++) {
            flat[i * b + j] = block[i][j];
        }
    }
    return flat;
}

BlockLUSimulator3D::BlockLUSimulator3D(const SimConfig3D& cfg)
    : config(cfg)
    , pe_array(cfg)
    , memory(cfg)  // Pass config, not just matrix_size
    , stats(cfg.num_layers)
{
}

std::vector<uint32_t> BlockLUSimulator3D::getLayerRows(uint32_t layer, uint32_t total_rows) {
    std::vector<uint32_t> rows;
    
    // Row distribution: row j goes to layer (j % num_layers)
    for (uint32_t j = 0; j < total_rows; j++) {
        if ((j % config.num_layers) == layer) {
            rows.push_back(j);
        }
    }
    
    return rows;
}

uint64_t BlockLUSimulator3D::broadcastURowToAllLayers(const std::vector<float>& u_row) {
    // Broadcast U row data to all layers
    // Source is layer 0, broadcast to all others
    uint64_t max_latency = 0;
    
    for (uint32_t c = 0; c < u_row.size() && c < config.pe_array_size; c++) {
        uint64_t latency = pe_array.broadcastToAllLayers(0, 0, c, u_row[c]);
        max_latency = std::max(max_latency, latency);
    }
    
    stats.tsv_transfers += config.num_layers - 1;
    
    return max_latency;
}

uint64_t BlockLUSimulator3D::broadcastBlockToAllLayers(const std::vector<std::vector<float>>& block) {
    uint64_t total_cycles = 0;
    uint32_t b = block.size();
    
    // Broadcast each element of the block
    for (uint32_t r = 0; r < b && r < config.pe_array_size; r++) {
        for (uint32_t c = 0; c < b && c < config.pe_array_size; c++) {
            uint64_t latency = pe_array.broadcastToAllLayers(0, r, c, block[r][c]);
            total_cycles = std::max(total_cycles, latency);
        }
    }
    
    stats.tsv_transfers += (config.num_layers - 1) * b * b;
    
    return total_cycles;
}

void BlockLUSimulator3D::loadBlockToLayer(uint32_t layer, 
                                           const std::vector<std::vector<float>>& block,
                                           uint32_t start_row, uint32_t num_rows) {
    uint32_t b = config.block_size;
    
    for (uint32_t local_r = 0; local_r < num_rows && (start_row + local_r) < b; local_r++) {
        uint32_t global_r = start_row + local_r;
        for (uint32_t c = 0; c < b && c < config.pe_array_size; c++) {
            auto& pe = pe_array.getPE(layer, local_r, c);
            pe.reg_a = block[global_r][c];
            pe.startLoad(pe_array.getMemoryLoadDelay(layer, local_r, c));
        }
    }
}

void BlockLUSimulator3D::writeBlockFromLayer(uint32_t layer,
                                              std::vector<std::vector<float>>& block,
                                              uint32_t start_row, uint32_t num_rows) {
    uint32_t b = config.block_size;
    
    for (uint32_t local_r = 0; local_r < num_rows && (start_row + local_r) < b; local_r++) {
        uint32_t global_r = start_row + local_r;
        for (uint32_t c = 0; c < b && c < config.pe_array_size; c++) {
            auto& pe = pe_array.getPE(layer, local_r, c);
            block[global_r][c] = pe.reg_result;
            pe.startWrite(pe_array.getMemoryWriteDelay(layer, local_r, c));
        }
    }
}

void BlockLUSimulator3D::executeCase1(uint32_t block_k) {
    // Case 1: Diagonal block LU decomposition
    // Limited 3D parallelism due to sequential pivot dependencies
    // We process on layer 0 only (same as 2D)
    
    uint64_t start_cycle = pe_array.current_cycle;
    uint32_t b = config.block_size;
    
    // Get diagonal block
    auto A_flat = memory.getBlockA(block_k, block_k);
    auto A_block = flatTo2D(A_flat, b);
    
    // Initialize L and U blocks
    std::vector<std::vector<float>> L_block(b, std::vector<float>(b, 0.0f));
    std::vector<std::vector<float>> U_block(b, std::vector<float>(b, 0.0f));
    
    // Initialize L diagonal to 1
    for (uint32_t i = 0; i < b; i++) {
        L_block[i][i] = 1.0f;
    }
    
    // Load A block to layer 0
    for (uint32_t r = 0; r < b && r < config.pe_array_size; r++) {
        for (uint32_t c = 0; c < b && c < config.pe_array_size; c++) {
            auto& pe = pe_array.getPE(0, r, c);
            pe.reg_a = A_block[r][c];
            pe.startLoad(config.mem_load_delay);
        }
    }
    pe_array.waitUntilIdle();
    
    // Doolittle algorithm with sequential pivot dependencies
    for (uint32_t k = 0; k < b; k++) {
        float pivot = A_block[k][k];
        U_block[k][k] = pivot;
        
        // Column k of L: L[j,k] = A[j,k] / pivot for j > k
        for (uint32_t j = k + 1; j < b; j++) {
            auto& pe = pe_array.getPE(0, j, k);
            pe.reg_a = A_block[j][k];
            pe.reg_b = pivot;
            pe.startDIV(A_block[j][k], pivot, config.div_latency);
            stats.div_operations++;
        }
        pe_array.waitUntilIdle();
        
        // Store L column results
        for (uint32_t j = k + 1; j < b; j++) {
            L_block[j][k] = pe_array.getPE(0, j, k).reg_result;
        }
        
        // Row k of U: U[k,l] = A[k,l] for l >= k
        for (uint32_t l = k + 1; l < b; l++) {
            U_block[k][l] = A_block[k][l];
        }
        
        // Update trailing submatrix: A[j,l] -= L[j,k] * U[k,l]
        for (uint32_t j = k + 1; j < b; j++) {
            for (uint32_t l = k + 1; l < b; l++) {
                auto& pe = pe_array.getPE(0, j, l);
                pe.startMAC(A_block[j][l], L_block[j][k], U_block[k][l], config.mac_latency);
                stats.mac_operations++;
            }
        }
        pe_array.waitUntilIdle();
        
        // Update A_block with results
        for (uint32_t j = k + 1; j < b; j++) {
            for (uint32_t l = k + 1; l < b; l++) {
                A_block[j][l] = pe_array.getPE(0, j, l).reg_result;
            }
        }
    }
    
    // Write results back
    for (uint32_t r = 0; r < b && r < config.pe_array_size; r++) {
        for (uint32_t c = 0; c < b && c < config.pe_array_size; c++) {
            pe_array.getPE(0, r, c).startWrite(config.mem_write_delay);
        }
    }
    pe_array.waitUntilIdle();
    
    // Store L and U blocks
    memory.setBlockL(block_k, block_k, toFlat(L_block, b));
    memory.setBlockU(block_k, block_k, toFlat(U_block, b));
    
    stats.case1_cycles += (pe_array.current_cycle - start_cycle);
    
    if (config.verbose) {
        std::cout << "Case 1 [" << block_k << "," << block_k << "]: "
                  << (pe_array.current_cycle - start_cycle) << " cycles\n";
    }
}

void BlockLUSimulator3D::executeCase2(uint32_t block_k, uint32_t block_j) {
    // Case 2: Horizontal U block update (forward substitution)
    // U = L^(-1) * A, where L is from diagonal block
    // Row distribution: different rows processed on different layers
    
    uint64_t start_cycle = pe_array.current_cycle;
    uint32_t b = config.block_size;
    
    auto A_flat = memory.getBlockA(block_k, block_j);
    auto L_flat = memory.getBlockL(block_k, block_k);
    auto A_block = flatTo2D(A_flat, b);
    auto L_block = flatTo2D(L_flat, b);
    std::vector<std::vector<float>> U_block(b, std::vector<float>(b, 0.0f));
    
    // Broadcast L block to all layers
    uint64_t broadcast_cycles = broadcastBlockToAllLayers(L_block);
    stats.tsv_transfer_cycles += broadcast_cycles;
    
    // Load A block - distribute rows across layers
    for (uint32_t z = 0; z < config.num_layers; z++) {
        auto layer_rows = getLayerRows(z, b);
        for (uint32_t local_r = 0; local_r < layer_rows.size(); local_r++) {
            uint32_t global_r = layer_rows[local_r];
            for (uint32_t c = 0; c < b && c < config.pe_array_size; c++) {
                auto& pe = pe_array.getPE(z, local_r, c);
                pe.reg_a = A_block[global_r][c];
                pe.startLoad(config.mem_load_delay);
            }
        }
    }
    pe_array.waitUntilIdle();
    
    // Forward substitution: process each row
    // For each k (column of L), update remaining columns
    for (uint32_t k = 0; k < b; k++) {
        // Row k of U = Row k of A (after previous updates)
        // This must be done sequentially due to dependencies
        
        // Each layer processes its assigned rows in parallel
        for (uint32_t z = 0; z < config.num_layers; z++) {
            auto layer_rows = getLayerRows(z, b);
            for (uint32_t local_r = 0; local_r < layer_rows.size(); local_r++) {
                uint32_t j = layer_rows[local_r];
                if (j <= k) continue;  // Only update rows below pivot
                
                for (uint32_t l = 0; l < b; l++) {
                    auto& pe = pe_array.getPE(z, local_r, l);
                    // A[j,l] = A[j,l] - L[j,k] * A[k,l]
                    // Note: A[k,l] is U[k,l] for previous iterations
                    float a_jl = A_block[j][l];
                    float l_jk = L_block[j][k];
                    float a_kl = A_block[k][l];
                    pe.startMAC(a_jl, l_jk, a_kl, config.mac_latency);
                    stats.mac_operations++;
                }
            }
        }
        pe_array.waitUntilIdle();
        
        // Collect results from all layers
        for (uint32_t z = 0; z < config.num_layers; z++) {
            auto layer_rows = getLayerRows(z, b);
            for (uint32_t local_r = 0; local_r < layer_rows.size(); local_r++) {
                uint32_t j = layer_rows[local_r];
                if (j <= k) continue;
                for (uint32_t l = 0; l < b; l++) {
                    A_block[j][l] = pe_array.getPE(z, local_r, l).reg_result;
                }
            }
        }
    }
    
    // Result is in A_block, copy to U_block
    U_block = A_block;
    
    // Write results
    for (uint32_t z = 0; z < config.num_layers; z++) {
        auto layer_rows = getLayerRows(z, b);
        for (uint32_t local_r = 0; local_r < layer_rows.size(); local_r++) {
            for (uint32_t c = 0; c < b && c < config.pe_array_size; c++) {
                pe_array.getPE(z, local_r, c).startWrite(config.mem_write_delay);
            }
        }
    }
    pe_array.waitUntilIdle();
    
    memory.setBlockU(block_k, block_j, toFlat(U_block, b));
    
    stats.case2_cycles += (pe_array.current_cycle - start_cycle);
    
    if (config.verbose) {
        std::cout << "Case 2 [" << block_k << "," << block_j << "]: "
                  << (pe_array.current_cycle - start_cycle) << " cycles\n";
    }
}

void BlockLUSimulator3D::executeCase3(uint32_t block_k, uint32_t block_i) {
    // Case 3: Vertical L block update (backward substitution variant)
    // L = A * U^(-1), where U is from diagonal block
    // ROW DISTRIBUTION: This is where 3D shines - rows are independent!
    
    uint64_t start_cycle = pe_array.current_cycle;
    uint32_t b = config.block_size;
    
    auto A_flat = memory.getBlockA(block_i, block_k);
    auto U_flat = memory.getBlockU(block_k, block_k);
    auto A_block = flatTo2D(A_flat, b);
    auto U_block = flatTo2D(U_flat, b);
    std::vector<std::vector<float>> L_block(b, std::vector<float>(b, 0.0f));
    
    // Broadcast U block to all layers via TSV
    uint64_t broadcast_cycles = broadcastBlockToAllLayers(U_block);
    stats.tsv_transfer_cycles += broadcast_cycles;
    
    // Simulate broadcast delay
    for (uint64_t i = 0; i < broadcast_cycles; i++) {
        pe_array.tick();
    }
    
    // Load A block - distribute rows across layers
    for (uint32_t z = 0; z < config.num_layers; z++) {
        auto layer_rows = getLayerRows(z, b);
        for (uint32_t local_r = 0; local_r < layer_rows.size(); local_r++) {
            uint32_t global_r = layer_rows[local_r];
            for (uint32_t c = 0; c < b && c < config.pe_array_size; c++) {
                auto& pe = pe_array.getPE(z, local_r, c);
                pe.reg_a = A_block[global_r][c];
                pe.startLoad(config.mem_load_delay);
            }
        }
    }
    pe_array.waitUntilIdle();
    
    // For each column k of L (and U), compute L column
    // Key insight: ALL ROWS ARE INDEPENDENT - perfect for 3D parallelism
    for (uint32_t k = 0; k < b; k++) {
        float pivot = U_block[k][k];
        
        // Step 1: Division - L[j,k] = A[j,k] / U[k,k]
        // All rows j can be processed IN PARALLEL across layers
        for (uint32_t z = 0; z < config.num_layers; z++) {
            auto layer_rows = getLayerRows(z, b);
            for (uint32_t local_r = 0; local_r < layer_rows.size(); local_r++) {
                uint32_t j = layer_rows[local_r];
                auto& pe = pe_array.getPE(z, local_r, k);
                pe.startDIV(A_block[j][k], pivot, config.div_latency);
                stats.div_operations++;
            }
        }
        pe_array.waitUntilIdle();
        
        // Collect L column results from all layers
        for (uint32_t z = 0; z < config.num_layers; z++) {
            auto layer_rows = getLayerRows(z, b);
            for (uint32_t local_r = 0; local_r < layer_rows.size(); local_r++) {
                uint32_t j = layer_rows[local_r];
                L_block[j][k] = pe_array.getPE(z, local_r, k).reg_result;
            }
        }
        
        // Step 2: Update remaining columns - A[j,l] -= L[j,k] * U[k,l]
        // Again, all rows j are independent - process in parallel
        for (uint32_t z = 0; z < config.num_layers; z++) {
            auto layer_rows = getLayerRows(z, b);
            for (uint32_t local_r = 0; local_r < layer_rows.size(); local_r++) {
                uint32_t j = layer_rows[local_r];
                for (uint32_t l = k + 1; l < b; l++) {
                    auto& pe = pe_array.getPE(z, local_r, l);
                    pe.startMAC(A_block[j][l], L_block[j][k], U_block[k][l], config.mac_latency);
                    stats.mac_operations++;
                }
            }
        }
        pe_array.waitUntilIdle();
        
        // Collect updated A values from all layers
        for (uint32_t z = 0; z < config.num_layers; z++) {
            auto layer_rows = getLayerRows(z, b);
            for (uint32_t local_r = 0; local_r < layer_rows.size(); local_r++) {
                uint32_t j = layer_rows[local_r];
                for (uint32_t l = k + 1; l < b; l++) {
                    A_block[j][l] = pe_array.getPE(z, local_r, l).reg_result;
                }
            }
        }
    }
    
    // Write results back - each layer writes its rows
    for (uint32_t z = 0; z < config.num_layers; z++) {
        auto layer_rows = getLayerRows(z, b);
        for (uint32_t local_r = 0; local_r < layer_rows.size(); local_r++) {
            for (uint32_t c = 0; c < b && c < config.pe_array_size; c++) {
                pe_array.getPE(z, local_r, c).startWrite(config.mem_write_delay);
            }
        }
    }
    pe_array.waitUntilIdle();
    
    memory.setBlockL(block_i, block_k, toFlat(L_block, b));
    
    stats.case3_cycles += (pe_array.current_cycle - start_cycle);
    
    if (config.verbose) {
        std::cout << "Case 3 [" << block_i << "," << block_k << "]: "
                  << (pe_array.current_cycle - start_cycle) << " cycles"
                  << " (3D parallel across " << config.num_layers << " layers)\n";
    }
}

void BlockLUSimulator3D::executeCase4(uint32_t block_k, uint32_t block_i, uint32_t block_j) {
    // Case 4: Trailing matrix Schur complement update
    // A_trail = A_trail - L_sub * U_sub
    // ROW DISTRIBUTION: Rows of A are distributed across layers
    
    uint64_t start_cycle = pe_array.current_cycle;
    uint32_t b = config.block_size;
    
    auto A_flat = memory.getBlockA(block_i, block_j);
    auto L_flat = memory.getBlockL(block_i, block_k);
    auto U_flat = memory.getBlockU(block_k, block_j);
    auto A_block = flatTo2D(A_flat, b);
    auto L_block = flatTo2D(L_flat, b);
    auto U_block = flatTo2D(U_flat, b);
    
    // Broadcast U block to all layers (L rows are distributed)
    uint64_t broadcast_cycles = broadcastBlockToAllLayers(U_block);
    stats.tsv_transfer_cycles += broadcast_cycles;
    
    // Simulate broadcast delay
    for (uint64_t i = 0; i < broadcast_cycles; i++) {
        pe_array.tick();
    }
    
    // Load A and L blocks - distribute rows across layers
    for (uint32_t z = 0; z < config.num_layers; z++) {
        auto layer_rows = getLayerRows(z, b);
        for (uint32_t local_r = 0; local_r < layer_rows.size(); local_r++) {
            uint32_t global_r = layer_rows[local_r];
            for (uint32_t c = 0; c < b && c < config.pe_array_size; c++) {
                auto& pe = pe_array.getPE(z, local_r, c);
                pe.reg_a = A_block[global_r][c];
                pe.reg_b = L_block[global_r][c];  // Store L row for this PE
                pe.startLoad(config.mem_load_delay);
            }
        }
    }
    pe_array.waitUntilIdle();
    
    // Compute A = A - L * U (outer product)
    // For each k: A[j,l] -= L[j,k] * U[k,l]
    // Rows j are distributed across layers - each layer works independently
    for (uint32_t k = 0; k < b; k++) {
        for (uint32_t z = 0; z < config.num_layers; z++) {
            auto layer_rows = getLayerRows(z, b);
            for (uint32_t local_r = 0; local_r < layer_rows.size(); local_r++) {
                uint32_t j = layer_rows[local_r];
                for (uint32_t l = 0; l < b; l++) {
                    auto& pe = pe_array.getPE(z, local_r, l);
                    pe.startMAC(A_block[j][l], L_block[j][k], U_block[k][l], config.mac_latency);
                    stats.mac_operations++;
                }
            }
        }
        pe_array.waitUntilIdle();
        
        // Collect results from all layers
        for (uint32_t z = 0; z < config.num_layers; z++) {
            auto layer_rows = getLayerRows(z, b);
            for (uint32_t local_r = 0; local_r < layer_rows.size(); local_r++) {
                uint32_t j = layer_rows[local_r];
                for (uint32_t l = 0; l < b; l++) {
                    A_block[j][l] = pe_array.getPE(z, local_r, l).reg_result;
                }
            }
        }
    }
    
    // Write results back
    for (uint32_t z = 0; z < config.num_layers; z++) {
        auto layer_rows = getLayerRows(z, b);
        for (uint32_t local_r = 0; local_r < layer_rows.size(); local_r++) {
            for (uint32_t c = 0; c < b && c < config.pe_array_size; c++) {
                pe_array.getPE(z, local_r, c).startWrite(config.mem_write_delay);
            }
        }
    }
    pe_array.waitUntilIdle();
    
    memory.setBlockA(block_i, block_j, toFlat(A_block, b));
    
    stats.case4_cycles += (pe_array.current_cycle - start_cycle);
    
    if (config.verbose) {
        std::cout << "Case 4 [" << block_i << "," << block_j << "]: "
                  << (pe_array.current_cycle - start_cycle) << " cycles"
                  << " (3D parallel across " << config.num_layers << " layers)\n";
    }
}

void BlockLUSimulator3D::run() {
    if (config.verbose) {
        std::cout << "\n=== Starting 3D Block LU Decomposition ===\n";
        std::cout << "Layers: " << config.num_layers << "\n";
        std::cout << "Row distribution: row j -> layer (j % " << config.num_layers << ")\n\n";
    }
    
    // Reset state
    pe_array.reset();
    memory.initializeRandom(42);  // Fixed seed for reproducibility
    memory.initializeLU();
    
    // Store original A for verification
    A_original = memory.A;
    
    stats = SimStats3D(config.num_layers);
    
    uint32_t num_blocks = config.matrix_size / config.block_size;
    
    // Block LU decomposition main loop
    for (uint32_t k = 0; k < num_blocks; k++) {
        if (config.verbose) {
            std::cout << "\n--- Block iteration k = " << k << " ---\n";
        }
        
        // Case 1: Diagonal block LU
        executeCase1(k);
        
        // Case 2: Horizontal U blocks
        for (uint32_t j = k + 1; j < num_blocks; j++) {
            executeCase2(k, j);
        }
        
        // Case 3: Vertical L blocks
        for (uint32_t i = k + 1; i < num_blocks; i++) {
            executeCase3(k, i);
        }
        
        // Case 4: Trailing matrix updates
        for (uint32_t i = k + 1; i < num_blocks; i++) {
            for (uint32_t j = k + 1; j < num_blocks; j++) {
                executeCase4(k, i, j);
            }
        }
    }
    
    // Collect final statistics
    stats.total_cycles = pe_array.current_cycle;
    pe_array.collectStats(stats);
    
    if (config.verbose) {
        std::cout << "\n=== 3D Simulation Complete ===\n";
    }
}

bool BlockLUSimulator3D::verify(float tolerance) {
    return memory.verify(A_original, tolerance);
}

void BlockLUSimulator3D::printResults() {
    config.print();
    stats.print();
}
