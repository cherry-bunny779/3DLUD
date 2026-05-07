#include "block_lu_simulator_3d.hpp"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <map>

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
    , memory(cfg)
    , stats(cfg.num_layers)
{
}

// =============================================================================
// HELPER METHODS
// =============================================================================

std::vector<uint32_t> BlockLUSimulator3D::getLayerRows(uint32_t layer, uint32_t total_rows) {
    std::vector<uint32_t> rows;
    for (uint32_t j = 0; j < total_rows; j++) {
        if ((j % config.num_layers) == layer) {
            rows.push_back(j);
        }
    }
    return rows;
}

std::vector<uint32_t> BlockLUSimulator3D::getLayerBlockIndices(uint32_t layer, uint32_t total_blocks) {
    std::vector<uint32_t> blocks;
    for (uint32_t b = 0; b < total_blocks; b++) {
        if ((b % config.num_layers) == layer) {
            blocks.push_back(b);
        }
    }
    return blocks;
}

std::vector<std::pair<uint32_t, uint32_t>> BlockLUSimulator3D::getLayerTrailingBlocks(
    uint32_t layer, uint32_t block_k, uint32_t num_blocks) {
    
    std::vector<std::pair<uint32_t, uint32_t>> blocks;
    
    // Distribute trailing blocks in row-major order
    uint32_t block_idx = 0;
    for (uint32_t i = block_k + 1; i < num_blocks; i++) {
        for (uint32_t j = block_k + 1; j < num_blocks; j++) {
            if ((block_idx % config.num_layers) == layer) {
                blocks.push_back({i, j});
            }
            block_idx++;
        }
    }
    return blocks;
}

uint64_t BlockLUSimulator3D::broadcastURowToAllLayers(const std::vector<float>& u_row) {
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
    for (uint32_t r = 0; r < b && r < config.pe_array_size; r++) {
        for (uint32_t c = 0; c < b && c < config.pe_array_size; c++) {
            uint64_t latency = pe_array.broadcastToAllLayers(0, r, c, block[r][c]);
            total_cycles = std::max(total_cycles, latency);
        }
    }
    stats.tsv_transfers += (config.num_layers - 1) * b * b;
    return total_cycles;
}

uint64_t BlockLUSimulator3D::broadcastBlockToLayer(uint32_t dst_layer,
                                                    const std::vector<std::vector<float>>& block) {
    if (dst_layer == 0) return 0;  // No transfer needed
    
    uint64_t latency = config.getTSVLatency(0, dst_layer);
    uint32_t b = block.size();
    
    // Load data to destination layer's PEs
    for (uint32_t r = 0; r < b && r < config.pe_array_size; r++) {
        for (uint32_t c = 0; c < b && c < config.pe_array_size; c++) {
            auto& pe = pe_array.getPE(dst_layer, r, c);
            pe.reg_a = block[r][c];
        }
    }
    
    stats.tsv_transfers++;
    stats.tsv_transfer_cycles += latency;
    
    return latency;
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

void BlockLUSimulator3D::loadFullBlockToLayer(uint32_t layer,
                                               const std::vector<std::vector<float>>& block) {
    uint32_t b = block.size();
    for (uint32_t r = 0; r < b && r < config.pe_array_size; r++) {
        for (uint32_t c = 0; c < b && c < config.pe_array_size; c++) {
            auto& pe = pe_array.getPE(layer, r, c);
            pe.reg_a = block[r][c];
            pe.startLoad(config.mem_load_delay);
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

void BlockLUSimulator3D::writeFullBlockFromLayer(uint32_t layer,
                                                  std::vector<std::vector<float>>& block) {
    uint32_t b = block.size();
    for (uint32_t r = 0; r < b && r < config.pe_array_size; r++) {
        for (uint32_t c = 0; c < b && c < config.pe_array_size; c++) {
            auto& pe = pe_array.getPE(layer, r, c);
            block[r][c] = pe.reg_result;
            pe.startWrite(config.mem_write_delay);
        }
    }
}

// =============================================================================
// BLOCK DISTRIBUTION SCHEME (Active)
// =============================================================================

void BlockLUSimulator3D::executeCase1_block(uint32_t block_k) {
    // Case 1: Diagonal block LU decomposition
    // Execute on Layer 0 only due to sequential pivot dependencies
    
    uint64_t start_cycle = pe_array.current_cycle;
    uint32_t b = config.block_size;
    
    auto A_flat = memory.getBlockA(block_k, block_k);
    auto A_block = flatTo2D(A_flat, b);
    
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
                  << (pe_array.current_cycle - start_cycle) << " cycles (Layer 0 only)\n";
    }
}

void BlockLUSimulator3D::executeSingleCase2OnLayer(uint32_t layer, uint32_t block_k, uint32_t block_j,
                                                    const std::vector<std::vector<float>>& L_diag) {
    // Execute forward substitution for U[block_k, block_j] on specified layer
    // This helper is kept for potential future use but computation is inlined in executeCase2_block
}

void BlockLUSimulator3D::executeCase2_block(uint32_t block_k, uint32_t num_blocks) {
    // Case 2: Compute all horizontal U blocks in parallel
    // U[k, k+1], U[k, k+2], ..., U[k, num_blocks-1]
    // Each layer processes different blocks
    
    uint64_t start_cycle = pe_array.current_cycle;
    uint32_t b = config.block_size;
    uint32_t num_horizontal = num_blocks - block_k - 1;
    
    if (num_horizontal == 0) return;
    
    // Get L[k,k] from diagonal (computed in Case 1)
    auto L_flat = memory.getBlockL(block_k, block_k);
    auto L_diag = flatTo2D(L_flat, b);
    
    // Broadcast L_diag to all active layers via TSV
    uint64_t broadcast_cycles = 0;
    uint32_t active_layers = std::min(num_horizontal, config.num_layers);
    
    for (uint32_t z = 1; z < active_layers; z++) {
        uint64_t latency = config.getTSVLatency(0, z);
        broadcast_cycles = std::max(broadcast_cycles, latency);
        stats.tsv_transfers++;
    }
    stats.tsv_transfer_cycles += broadcast_cycles;
    
    // Simulate broadcast delay
    for (uint64_t i = 0; i < broadcast_cycles; i++) {
        pe_array.tick();
    }
    
    // Assign blocks to layers
    // Layer z processes blocks {k+1+z, k+1+z+Z, k+1+z+2Z, ...}
    std::vector<std::vector<std::vector<std::vector<float>>>> A_blocks(config.num_layers);
    std::vector<std::vector<uint32_t>> layer_block_j(config.num_layers);
    
    for (uint32_t j = block_k + 1; j < num_blocks; j++) {
        uint32_t z = (j - block_k - 1) % config.num_layers;
        layer_block_j[z].push_back(j);
    }
    
    // Process blocks in rounds (each layer processes one block at a time)
    uint32_t max_blocks_per_layer = (num_horizontal + config.num_layers - 1) / config.num_layers;
    
    for (uint32_t round = 0; round < max_blocks_per_layer; round++) {
        // Load A blocks for this round
        for (uint32_t z = 0; z < config.num_layers; z++) {
            if (round < layer_block_j[z].size()) {
                uint32_t j = layer_block_j[z][round];
                auto A_flat = memory.getBlockA(block_k, j);
                auto A_block = flatTo2D(A_flat, b);
                
                for (uint32_t r = 0; r < b && r < config.pe_array_size; r++) {
                    for (uint32_t c = 0; c < b && c < config.pe_array_size; c++) {
                        auto& pe = pe_array.getPE(z, r, c);
                        pe.reg_a = A_block[r][c];
                        pe.startLoad(config.mem_load_delay);
                    }
                }
                
                // Store for later use
                if (A_blocks[z].size() <= round) {
                    A_blocks[z].resize(round + 1);
                }
                A_blocks[z][round] = A_block;
            }
        }
        pe_array.waitUntilIdle();
        
        // Forward substitution on all layers in parallel
        for (uint32_t k = 0; k < b; k++) {
            // All layers process their blocks in parallel
            for (uint32_t z = 0; z < config.num_layers; z++) {
                if (round < layer_block_j[z].size()) {
                    auto& A_block = A_blocks[z][round];
                    
                    for (uint32_t row = k + 1; row < b; row++) {
                        for (uint32_t col = 0; col < b; col++) {
                            auto& pe = pe_array.getPE(z, row, col);
                            float a_jl = A_block[row][col];
                            float l_jk = L_diag[row][k];
                            float a_kl = A_block[k][col];
                            pe.startMAC(a_jl, l_jk, a_kl, config.mac_latency);
                            stats.mac_operations++;
                        }
                    }
                }
            }
            pe_array.waitUntilIdle();
            
            // Collect results from all layers
            for (uint32_t z = 0; z < config.num_layers; z++) {
                if (round < layer_block_j[z].size()) {
                    auto& A_block = A_blocks[z][round];
                    for (uint32_t row = k + 1; row < b; row++) {
                        for (uint32_t col = 0; col < b; col++) {
                            A_block[row][col] = pe_array.getPE(z, row, col).reg_result;
                        }
                    }
                }
            }
        }
        
        // Write results back
        for (uint32_t z = 0; z < config.num_layers; z++) {
            if (round < layer_block_j[z].size()) {
                uint32_t j = layer_block_j[z][round];
                auto& U_block = A_blocks[z][round];
                
                for (uint32_t r = 0; r < b && r < config.pe_array_size; r++) {
                    for (uint32_t c = 0; c < b && c < config.pe_array_size; c++) {
                        pe_array.getPE(z, r, c).startWrite(config.mem_write_delay);
                    }
                }
                
                memory.setBlockU(block_k, j, toFlat(U_block, b));
            }
        }
        pe_array.waitUntilIdle();
    }
    
    stats.case2_cycles += (pe_array.current_cycle - start_cycle);
    
    if (config.verbose) {
        std::cout << "Case 2 [k=" << block_k << "]: "
                  << (pe_array.current_cycle - start_cycle) << " cycles"
                  << " (" << num_horizontal << " blocks across " 
                  << std::min(num_horizontal, config.num_layers) << " layers)\n";
    }
}

void BlockLUSimulator3D::executeSingleCase3OnLayer(uint32_t layer, uint32_t block_k, uint32_t block_i,
                                                    const std::vector<std::vector<float>>& U_diag) {
    // Execute column-wise L computation for L[block_i, block_k] on specified layer
    // This helper is kept for potential future use but computation is inlined in executeCase3_block
}

void BlockLUSimulator3D::executeCase3_block(uint32_t block_k, uint32_t num_blocks) {
    // Case 3: Compute all vertical L blocks in parallel
    // L[k+1, k], L[k+2, k], ..., L[num_blocks-1, k]
    // Each layer processes different blocks
    
    uint64_t start_cycle = pe_array.current_cycle;
    uint32_t b = config.block_size;
    uint32_t num_vertical = num_blocks - block_k - 1;
    
    if (num_vertical == 0) return;
    
    // Get U[k,k] from diagonal (computed in Case 1)
    auto U_flat = memory.getBlockU(block_k, block_k);
    auto U_diag = flatTo2D(U_flat, b);
    
    // Broadcast U_diag to all active layers via TSV
    uint64_t broadcast_cycles = 0;
    uint32_t active_layers = std::min(num_vertical, config.num_layers);
    
    for (uint32_t z = 1; z < active_layers; z++) {
        uint64_t latency = config.getTSVLatency(0, z);
        broadcast_cycles = std::max(broadcast_cycles, latency);
        stats.tsv_transfers++;
    }
    stats.tsv_transfer_cycles += broadcast_cycles;
    
    // Simulate broadcast delay
    for (uint64_t i = 0; i < broadcast_cycles; i++) {
        pe_array.tick();
    }
    
    // Assign blocks to layers
    std::vector<std::vector<uint32_t>> layer_block_i(config.num_layers);
    for (uint32_t i = block_k + 1; i < num_blocks; i++) {
        uint32_t z = (i - block_k - 1) % config.num_layers;
        layer_block_i[z].push_back(i);
    }
    
    // Storage for A and L blocks per layer
    std::vector<std::vector<std::vector<std::vector<float>>>> A_blocks(config.num_layers);
    std::vector<std::vector<std::vector<std::vector<float>>>> L_blocks(config.num_layers);
    
    // Process blocks in rounds
    uint32_t max_blocks_per_layer = (num_vertical + config.num_layers - 1) / config.num_layers;
    
    for (uint32_t round = 0; round < max_blocks_per_layer; round++) {
        // Load A blocks for this round
        for (uint32_t z = 0; z < config.num_layers; z++) {
            if (round < layer_block_i[z].size()) {
                uint32_t i = layer_block_i[z][round];
                auto A_flat = memory.getBlockA(i, block_k);
                auto A_block = flatTo2D(A_flat, b);
                
                for (uint32_t r = 0; r < b && r < config.pe_array_size; r++) {
                    for (uint32_t c = 0; c < b && c < config.pe_array_size; c++) {
                        auto& pe = pe_array.getPE(z, r, c);
                        pe.reg_a = A_block[r][c];
                        pe.startLoad(config.mem_load_delay);
                    }
                }
                
                if (A_blocks[z].size() <= round) {
                    A_blocks[z].resize(round + 1);
                    L_blocks[z].resize(round + 1);
                }
                A_blocks[z][round] = A_block;
                L_blocks[z][round] = std::vector<std::vector<float>>(b, std::vector<float>(b, 0.0f));
            }
        }
        pe_array.waitUntilIdle();
        
        // Column-by-column computation (this is where 3D parallelism helps)
        for (uint32_t k = 0; k < b; k++) {
            float pivot = U_diag[k][k];
            
            // Division: L[j,k] = A[j,k] / U[k,k] - ALL ROWS PARALLEL
            for (uint32_t z = 0; z < config.num_layers; z++) {
                if (round < layer_block_i[z].size()) {
                    auto& A_block = A_blocks[z][round];
                    
                    for (uint32_t j = 0; j < b; j++) {
                        auto& pe = pe_array.getPE(z, j, k);
                        pe.startDIV(A_block[j][k], pivot, config.div_latency);
                        stats.div_operations++;
                    }
                }
            }
            pe_array.waitUntilIdle();
            
            // Collect L column results
            for (uint32_t z = 0; z < config.num_layers; z++) {
                if (round < layer_block_i[z].size()) {
                    auto& L_block = L_blocks[z][round];
                    for (uint32_t j = 0; j < b; j++) {
                        L_block[j][k] = pe_array.getPE(z, j, k).reg_result;
                    }
                }
            }
            
            // MAC update: A[j,l] -= L[j,k] * U[k,l] - ALL ROWS PARALLEL
            for (uint32_t z = 0; z < config.num_layers; z++) {
                if (round < layer_block_i[z].size()) {
                    auto& A_block = A_blocks[z][round];
                    auto& L_block = L_blocks[z][round];
                    
                    for (uint32_t j = 0; j < b; j++) {
                        for (uint32_t l = k + 1; l < b; l++) {
                            auto& pe = pe_array.getPE(z, j, l);
                            pe.startMAC(A_block[j][l], L_block[j][k], U_diag[k][l], config.mac_latency);
                            stats.mac_operations++;
                        }
                    }
                }
            }
            pe_array.waitUntilIdle();
            
            // Collect updated A values
            for (uint32_t z = 0; z < config.num_layers; z++) {
                if (round < layer_block_i[z].size()) {
                    auto& A_block = A_blocks[z][round];
                    for (uint32_t j = 0; j < b; j++) {
                        for (uint32_t l = k + 1; l < b; l++) {
                            A_block[j][l] = pe_array.getPE(z, j, l).reg_result;
                        }
                    }
                }
            }
        }
        
        // Write L results back
        for (uint32_t z = 0; z < config.num_layers; z++) {
            if (round < layer_block_i[z].size()) {
                uint32_t i = layer_block_i[z][round];
                auto& L_block = L_blocks[z][round];
                
                for (uint32_t r = 0; r < b && r < config.pe_array_size; r++) {
                    for (uint32_t c = 0; c < b && c < config.pe_array_size; c++) {
                        pe_array.getPE(z, r, c).startWrite(config.mem_write_delay);
                    }
                }
                
                memory.setBlockL(i, block_k, toFlat(L_block, b));
            }
        }
        pe_array.waitUntilIdle();
    }
    
    stats.case3_cycles += (pe_array.current_cycle - start_cycle);
    
    if (config.verbose) {
        std::cout << "Case 3 [k=" << block_k << "]: "
                  << (pe_array.current_cycle - start_cycle) << " cycles"
                  << " (" << num_vertical << " blocks across "
                  << std::min(num_vertical, config.num_layers) << " layers)\n";
    }
}

void BlockLUSimulator3D::executeSingleCase4OnLayer(uint32_t layer, uint32_t block_k,
                                                    uint32_t block_i, uint32_t block_j) {
    // Execute Schur complement for A[block_i, block_j] on specified layer
    // This helper is kept for potential future use but computation is inlined in executeCase4_block
}

void BlockLUSimulator3D::executeCase4_block(uint32_t block_k, uint32_t num_blocks) {
    // Case 4: Compute all trailing block updates in parallel
    // A[i,j] = A[i,j] - L[i,k] * U[k,j] for all i,j > k
    // Each layer processes different (i,j) blocks
    
    uint64_t start_cycle = pe_array.current_cycle;
    uint32_t b = config.block_size;
    uint32_t trailing_size = num_blocks - block_k - 1;
    uint32_t total_trailing = trailing_size * trailing_size;
    
    if (total_trailing == 0) return;
    
    // Get all L column blocks and U row blocks needed
    // L[k+1,k], L[k+2,k], ..., L[num_blocks-1,k]
    // U[k,k+1], U[k,k+2], ..., U[k,num_blocks-1]
    std::vector<std::vector<std::vector<float>>> L_col_blocks(trailing_size);
    std::vector<std::vector<std::vector<float>>> U_row_blocks(trailing_size);
    
    for (uint32_t idx = 0; idx < trailing_size; idx++) {
        uint32_t i = block_k + 1 + idx;
        auto L_flat = memory.getBlockL(i, block_k);
        L_col_blocks[idx] = flatTo2D(L_flat, b);
        
        auto U_flat = memory.getBlockU(block_k, i);
        U_row_blocks[idx] = flatTo2D(U_flat, b);
    }
    
    // Broadcast all L and U blocks to all layers
    // (Simpler approach: broadcast everything, layers ignore what they don't need)
    uint64_t broadcast_cycles = 0;
    for (uint32_t z = 1; z < config.num_layers; z++) {
        uint64_t latency = config.getTSVLatency(0, z);
        broadcast_cycles = std::max(broadcast_cycles, latency);
    }
    stats.tsv_transfers += (config.num_layers - 1) * trailing_size * 2;  // L and U blocks
    stats.tsv_transfer_cycles += broadcast_cycles;
    
    // Simulate broadcast delay
    for (uint64_t i = 0; i < broadcast_cycles; i++) {
        pe_array.tick();
    }
    
    // Assign trailing blocks to layers in row-major order
    std::vector<std::vector<std::pair<uint32_t, uint32_t>>> layer_blocks(config.num_layers);
    uint32_t block_idx = 0;
    for (uint32_t i = block_k + 1; i < num_blocks; i++) {
        for (uint32_t j = block_k + 1; j < num_blocks; j++) {
            uint32_t z = block_idx % config.num_layers;
            layer_blocks[z].push_back({i, j});
            block_idx++;
        }
    }
    
    // Storage for A blocks per layer
    std::vector<std::vector<std::vector<std::vector<float>>>> A_blocks(config.num_layers);
    
    // Process blocks in rounds
    uint32_t max_blocks_per_layer = (total_trailing + config.num_layers - 1) / config.num_layers;
    
    for (uint32_t round = 0; round < max_blocks_per_layer; round++) {
        // Load A blocks for this round
        for (uint32_t z = 0; z < config.num_layers; z++) {
            if (round < layer_blocks[z].size()) {
                auto [bi, bj] = layer_blocks[z][round];
                auto A_flat = memory.getBlockA(bi, bj);
                auto A_block = flatTo2D(A_flat, b);
                
                for (uint32_t r = 0; r < b && r < config.pe_array_size; r++) {
                    for (uint32_t c = 0; c < b && c < config.pe_array_size; c++) {
                        auto& pe = pe_array.getPE(z, r, c);
                        pe.reg_a = A_block[r][c];
                        pe.startLoad(config.mem_load_delay);
                    }
                }
                
                if (A_blocks[z].size() <= round) {
                    A_blocks[z].resize(round + 1);
                }
                A_blocks[z][round] = A_block;
            }
        }
        pe_array.waitUntilIdle();
        
        // Compute A = A - L * U (outer product) - ALL BLOCKS PARALLEL
        for (uint32_t k = 0; k < b; k++) {
            for (uint32_t z = 0; z < config.num_layers; z++) {
                if (round < layer_blocks[z].size()) {
                    auto [bi, bj] = layer_blocks[z][round];
                    auto& A_block = A_blocks[z][round];
                    
                    // Get the corresponding L and U blocks
                    uint32_t L_idx = bi - block_k - 1;
                    uint32_t U_idx = bj - block_k - 1;
                    auto& L_block = L_col_blocks[L_idx];
                    auto& U_block = U_row_blocks[U_idx];
                    
                    for (uint32_t row = 0; row < b; row++) {
                        for (uint32_t col = 0; col < b; col++) {
                            auto& pe = pe_array.getPE(z, row, col);
                            pe.startMAC(A_block[row][col], L_block[row][k], U_block[k][col], 
                                       config.mac_latency);
                            stats.mac_operations++;
                        }
                    }
                }
            }
            pe_array.waitUntilIdle();
            
            // Collect results
            for (uint32_t z = 0; z < config.num_layers; z++) {
                if (round < layer_blocks[z].size()) {
                    auto& A_block = A_blocks[z][round];
                    for (uint32_t row = 0; row < b; row++) {
                        for (uint32_t col = 0; col < b; col++) {
                            A_block[row][col] = pe_array.getPE(z, row, col).reg_result;
                        }
                    }
                }
            }
        }
        
        // Write results back
        for (uint32_t z = 0; z < config.num_layers; z++) {
            if (round < layer_blocks[z].size()) {
                auto [bi, bj] = layer_blocks[z][round];
                auto& A_block = A_blocks[z][round];
                
                for (uint32_t r = 0; r < b && r < config.pe_array_size; r++) {
                    for (uint32_t c = 0; c < b && c < config.pe_array_size; c++) {
                        pe_array.getPE(z, r, c).startWrite(config.mem_write_delay);
                    }
                }
                
                memory.setBlockA(bi, bj, toFlat(A_block, b));
            }
        }
        pe_array.waitUntilIdle();
    }
    
    stats.case4_cycles += (pe_array.current_cycle - start_cycle);
    
    if (config.verbose) {
        std::cout << "Case 4 [k=" << block_k << "]: "
                  << (pe_array.current_cycle - start_cycle) << " cycles"
                  << " (" << total_trailing << " blocks across "
                  << std::min(total_trailing, config.num_layers) << " layers)\n";
    }
}

void BlockLUSimulator3D::executeCase4_block_pipelined_v1(uint32_t block_k, uint32_t num_blocks) {
    // LEGACY V1: Arbitrary block-to-layer assignment with pipelining within rows
    // Kept for comparison purposes
    //
    // Pipelined Case 4: Use skewed systolic timing within each layer
    // Blocks assigned to the same layer with the same row index (i) are pipelined
    // 
    // Timing model (classic systolic skewing):
    // - PE[row,col] receives data at cycle (row + col)
    // - Single block latency: 3*b - 2 cycles
    // - Pipeline interval: b cycles between consecutive blocks in a row
    
    uint64_t start_cycle = pe_array.current_cycle;
    uint32_t b = config.block_size;
    uint32_t trailing_size = num_blocks - block_k - 1;
    uint32_t total_trailing = trailing_size * trailing_size;
    
    if (total_trailing == 0) return;
    
    // Get all L column blocks and U row blocks needed
    std::vector<std::vector<std::vector<float>>> L_col_blocks(trailing_size);
    std::vector<std::vector<std::vector<float>>> U_row_blocks(trailing_size);
    
    for (uint32_t idx = 0; idx < trailing_size; idx++) {
        uint32_t i = block_k + 1 + idx;
        auto L_flat = memory.getBlockL(i, block_k);
        L_col_blocks[idx] = flatTo2D(L_flat, b);
        
        auto U_flat = memory.getBlockU(block_k, i);
        U_row_blocks[idx] = flatTo2D(U_flat, b);
    }
    
    // Broadcast L and U blocks to all layers (TSV overhead)
    uint64_t broadcast_cycles = 0;
    for (uint32_t z = 1; z < config.num_layers; z++) {
        uint64_t latency = config.getTSVLatency(0, z);
        broadcast_cycles = std::max(broadcast_cycles, latency);
    }
    stats.tsv_transfers += (config.num_layers - 1) * trailing_size * 2;
    stats.tsv_transfer_cycles += broadcast_cycles;
    
    for (uint64_t i = 0; i < broadcast_cycles; i++) {
        pe_array.tick();
    }
    
    // Assign trailing blocks to layers and group by row within each layer
    // layer_row_blocks[z][i_idx] = list of j indices for row i on layer z
    std::vector<std::map<uint32_t, std::vector<uint32_t>>> layer_row_blocks(config.num_layers);
    uint32_t block_idx = 0;
    
    for (uint32_t i = block_k + 1; i < num_blocks; i++) {
        for (uint32_t j = block_k + 1; j < num_blocks; j++) {
            uint32_t z = block_idx % config.num_layers;
            uint32_t i_idx = i - block_k - 1;
            layer_row_blocks[z][i_idx].push_back(j);
            block_idx++;
        }
    }
    
    // Calculate pipelined cycle count for each layer
    // Each layer processes its assigned rows; rows are processed sequentially,
    // but blocks within a row are pipelined
    uint64_t max_layer_cycles = 0;
    
    for (uint32_t z = 0; z < config.num_layers; z++) {
        if (layer_row_blocks[z].empty()) continue;
        
        uint64_t layer_cycles = 0;
        
        for (auto& [i_idx, j_list] : layer_row_blocks[z]) {
            // Memory load for first block in row (overlapped with prev row's compute for later rows)
            layer_cycles += config.mem_load_delay;
            
            // Pipelined execution of all blocks in this row
            // N blocks take: (N-1)*b + single_block_latency cycles
            uint32_t n_blocks = j_list.size();
            uint64_t single_block_latency = 3 * b - 2;  // Skewed systolic timing
            uint64_t row_compute_cycles = (n_blocks > 1) ? 
                ((n_blocks - 1) * b + single_block_latency) : single_block_latency;
            layer_cycles += row_compute_cycles;
            
            // Memory write (can overlap with next row's load)
            layer_cycles += config.mem_write_delay;
        }
        
        max_layer_cycles = std::max(max_layer_cycles, layer_cycles);
    }
    
    // All layers execute in parallel, so total time is max across layers
    // Add the cycles to the PE array
    for (uint64_t i = 0; i < max_layer_cycles; i++) {
        pe_array.tick();
    }
    
    // Actually compute the block values for correctness
    // (The timing is already accounted for above)
    for (uint32_t z = 0; z < config.num_layers; z++) {
        for (auto& [i_idx, j_list] : layer_row_blocks[z]) {
            uint32_t i = block_k + 1 + i_idx;
            auto& L_block = L_col_blocks[i_idx];
            
            for (uint32_t j : j_list) {
                uint32_t j_idx = j - block_k - 1;
                auto& U_block = U_row_blocks[j_idx];
                
                // Load A block
                auto A_flat = memory.getBlockA(i, j);
                auto A_block = flatTo2D(A_flat, b);
                
                // Compute A = A - L * U (outer product)
                for (uint32_t k = 0; k < b; k++) {
                    for (uint32_t row = 0; row < b; row++) {
                        for (uint32_t col = 0; col < b; col++) {
                            A_block[row][col] -= L_block[row][k] * U_block[k][col];
                            stats.mac_operations++;
                        }
                    }
                }
                
                // Store result
                memory.setBlockA(i, j, toFlat(A_block, b));
            }
        }
    }
    
    stats.case4_cycles += (pe_array.current_cycle - start_cycle);
    
    if (config.verbose) {
        std::cout << "Case 4 [k=" << block_k << "] (pipelined-v1): "
                  << (pe_array.current_cycle - start_cycle) << " cycles"
                  << " (" << total_trailing << " blocks across "
                  << std::min(total_trailing, config.num_layers) << " layers)\n";
    }
}

void BlockLUSimulator3D::executeCase4_block_pipelined(uint32_t block_k, uint32_t num_blocks) {
    // NEW DESIGN: Spatial + Temporal Pipelining
    //
    // Key insight: Each trailing update A^(i,j) -= L^(i,k) × U^(k,j) is independent.
    //
    // SPATIAL PARALLELISM (across layers):
    //   - Assign different ROW indices (i) to different layers
    //   - Layer z processes rows: i = k+1+z, k+1+z+Z, k+1+z+2Z, ...
    //   - Each layer has a DIFFERENT L block stationary
    //
    // TEMPORAL PIPELINING (within each layer):
    //   - Stream the SAME sequence of U blocks to ALL layers simultaneously
    //   - U^(k, k+1), U^(k, k+2), ..., U^(k, P-1) streamed with skewed timing
    //   - Each layer uses its stationary L with the streaming U blocks
    //
    // Result: All layers compute in parallel, each producing a full row of results
    //
    // Timing model:
    //   - L loading: One-time load per row assignment (b×b elements)
    //   - U streaming: Skewed entry, columns enter at offset j
    //   - Single block latency: 3b-2 cycles
    //   - Pipeline interval: b cycles between consecutive U blocks
    //   - N U blocks pipelined: (N-1)×b + (3b-2) = (N+2)b - 2 cycles
    
    uint64_t start_cycle = pe_array.current_cycle;
    uint32_t b = config.block_size;
    uint32_t trailing_size = num_blocks - block_k - 1;  // Number of rows/cols in trailing matrix
    
    if (trailing_size == 0) return;
    
    // Preload all L column blocks and U row blocks
    std::vector<std::vector<std::vector<float>>> L_col_blocks(trailing_size);  // L^(k+1+idx, k)
    std::vector<std::vector<std::vector<float>>> U_row_blocks(trailing_size);  // U^(k, k+1+idx)
    
    for (uint32_t idx = 0; idx < trailing_size; idx++) {
        uint32_t i = block_k + 1 + idx;
        auto L_flat = memory.getBlockL(i, block_k);
        L_col_blocks[idx] = flatTo2D(L_flat, b);
        
        auto U_flat = memory.getBlockU(block_k, i);
        U_row_blocks[idx] = flatTo2D(U_flat, b);
    }
    
    // Assign rows to layers in round-robin fashion
    // layer_rows[z] = list of row indices (0-based into trailing matrix) assigned to layer z
    std::vector<std::vector<uint32_t>> layer_rows(config.num_layers);
    for (uint32_t row_idx = 0; row_idx < trailing_size; row_idx++) {
        uint32_t z = row_idx % config.num_layers;
        layer_rows[z].push_back(row_idx);
    }
    
    // Count active layers
    uint32_t active_layers = std::min(trailing_size, config.num_layers);
    
    // TSV broadcast overhead for U blocks to all layers
    // All layers receive the same U block sequence
    uint64_t broadcast_cycles = 0;
    for (uint32_t z = 1; z < active_layers; z++) {
        uint64_t latency = config.getTSVLatency(0, z);
        broadcast_cycles = std::max(broadcast_cycles, latency);
    }
    stats.tsv_transfers += active_layers * trailing_size;  // U blocks broadcast
    stats.tsv_transfer_cycles += broadcast_cycles;
    
    for (uint64_t i = 0; i < broadcast_cycles; i++) {
        pe_array.tick();
    }
    
    // Calculate timing for each layer
    // All layers stream the SAME U block sequence, but may have different numbers of L rows
    uint64_t max_layer_cycles = 0;
    
    // U streaming parameters (same for all layers)
    uint64_t single_block_latency = 3 * b - 2;  // Skewed systolic timing for one L×U
    uint64_t u_stream_cycles = (trailing_size > 1) ? 
        ((trailing_size - 1) * b + single_block_latency) : single_block_latency;

    // Bandwidth contention slowdown (applies to streaming portion of each layer's work)
    // See SimConfig3D::getBandwidthSlowdownFactor for the demand/supply derivation.
    const double bw_slowdown = config.getBandwidthSlowdownFactor(active_layers);

    for (uint32_t z = 0; z < config.num_layers; z++) {
        if (layer_rows[z].empty()) continue;
        
        uint32_t num_row_batches = layer_rows[z].size();
        uint64_t layer_cycles = 0;
        
        // For each row assigned to this layer:
        // 1. Load L block for this row (mem_load_delay)
        // 2. Stream ALL U blocks through (u_stream_cycles, scaled by bandwidth slowdown)
        // 3. Store results for this row (mem_write_delay)
        // 
        // Rows within a layer are processed SEQUENTIALLY (different L blocks needed)
        for (uint32_t batch = 0; batch < num_row_batches; batch++) {
            layer_cycles += config.mem_load_delay;    // Load L and initial A values
            layer_cycles += static_cast<uint64_t>(
                static_cast<double>(u_stream_cycles) * bw_slowdown);  // Streaming, with contention
            layer_cycles += config.mem_write_delay;    // Store results
        }
        
        max_layer_cycles = std::max(max_layer_cycles, layer_cycles);
    }
    
    // All layers execute in parallel, so total time is max across layers
    for (uint64_t i = 0; i < max_layer_cycles; i++) {
        pe_array.tick();
    }
    
    // Actually compute the block values for correctness
    for (uint32_t z = 0; z < config.num_layers; z++) {
        for (uint32_t row_idx : layer_rows[z]) {
            uint32_t i = block_k + 1 + row_idx;
            auto& L_block = L_col_blocks[row_idx];
            
            // This layer processes row i with ALL U blocks
            for (uint32_t col_idx = 0; col_idx < trailing_size; col_idx++) {
                uint32_t j = block_k + 1 + col_idx;
                auto& U_block = U_row_blocks[col_idx];
                
                // Load A block
                auto A_flat = memory.getBlockA(i, j);
                auto A_block = flatTo2D(A_flat, b);
                
                // Compute A = A - L * U
                for (uint32_t k = 0; k < b; k++) {
                    for (uint32_t row = 0; row < b; row++) {
                        for (uint32_t col = 0; col < b; col++) {
                            A_block[row][col] -= L_block[row][k] * U_block[k][col];
                            stats.mac_operations++;
                        }
                    }
                }
                
                // Store result
                memory.setBlockA(i, j, toFlat(A_block, b));
            }
        }
    }
    
    stats.case4_cycles += (pe_array.current_cycle - start_cycle);
    
    if (config.verbose) {
        std::cout << "Case 4 [k=" << block_k << "] (pipelined-v2): "
                  << (pe_array.current_cycle - start_cycle) << " cycles"
                  << " (" << trailing_size << " rows × " << trailing_size << " cols, "
                  << active_layers << " layers active, "
                  << "U stream: " << u_stream_cycles << " cycles"
                  << ", BW slowdown: " << bw_slowdown << "x)\n";
    }
}

// =============================================================================
// ROW DISTRIBUTION SCHEME (Legacy - kept but unused)
// =============================================================================

void BlockLUSimulator3D::executeCase1_row(uint32_t block_k) {
    // Case 1: Diagonal block LU decomposition (same as block distribution)
    // Limited 3D parallelism due to sequential pivot dependencies
    // We process on layer 0 only (same as 2D)
    
    uint64_t start_cycle = pe_array.current_cycle;
    uint32_t b = config.block_size;
    
    auto A_flat = memory.getBlockA(block_k, block_k);
    auto A_block = flatTo2D(A_flat, b);
    
    std::vector<std::vector<float>> L_block(b, std::vector<float>(b, 0.0f));
    std::vector<std::vector<float>> U_block(b, std::vector<float>(b, 0.0f));
    
    for (uint32_t i = 0; i < b; i++) {
        L_block[i][i] = 1.0f;
    }
    
    for (uint32_t r = 0; r < b && r < config.pe_array_size; r++) {
        for (uint32_t c = 0; c < b && c < config.pe_array_size; c++) {
            auto& pe = pe_array.getPE(0, r, c);
            pe.reg_a = A_block[r][c];
            pe.startLoad(config.mem_load_delay);
        }
    }
    pe_array.waitUntilIdle();
    
    for (uint32_t k = 0; k < b; k++) {
        float pivot = A_block[k][k];
        U_block[k][k] = pivot;
        
        for (uint32_t j = k + 1; j < b; j++) {
            auto& pe = pe_array.getPE(0, j, k);
            pe.reg_a = A_block[j][k];
            pe.reg_b = pivot;
            pe.startDIV(A_block[j][k], pivot, config.div_latency);
            stats.div_operations++;
        }
        pe_array.waitUntilIdle();
        
        for (uint32_t j = k + 1; j < b; j++) {
            L_block[j][k] = pe_array.getPE(0, j, k).reg_result;
        }
        
        for (uint32_t l = k + 1; l < b; l++) {
            U_block[k][l] = A_block[k][l];
        }
        
        for (uint32_t j = k + 1; j < b; j++) {
            for (uint32_t l = k + 1; l < b; l++) {
                auto& pe = pe_array.getPE(0, j, l);
                pe.startMAC(A_block[j][l], L_block[j][k], U_block[k][l], config.mac_latency);
                stats.mac_operations++;
            }
        }
        pe_array.waitUntilIdle();
        
        for (uint32_t j = k + 1; j < b; j++) {
            for (uint32_t l = k + 1; l < b; l++) {
                A_block[j][l] = pe_array.getPE(0, j, l).reg_result;
            }
        }
    }
    
    for (uint32_t r = 0; r < b && r < config.pe_array_size; r++) {
        for (uint32_t c = 0; c < b && c < config.pe_array_size; c++) {
            pe_array.getPE(0, r, c).startWrite(config.mem_write_delay);
        }
    }
    pe_array.waitUntilIdle();
    
    memory.setBlockL(block_k, block_k, toFlat(L_block, b));
    memory.setBlockU(block_k, block_k, toFlat(U_block, b));
    
    stats.case1_cycles += (pe_array.current_cycle - start_cycle);
}

void BlockLUSimulator3D::executeCase2_row(uint32_t block_k, uint32_t block_j) {
    // Case 2: Horizontal U block update (forward substitution)
    // Row distribution: different rows processed on different layers
    
    uint64_t start_cycle = pe_array.current_cycle;
    uint32_t b = config.block_size;
    
    auto A_flat = memory.getBlockA(block_k, block_j);
    auto L_flat = memory.getBlockL(block_k, block_k);
    auto A_block = flatTo2D(A_flat, b);
    auto L_block = flatTo2D(L_flat, b);
    std::vector<std::vector<float>> U_block(b, std::vector<float>(b, 0.0f));
    
    uint64_t broadcast_cycles = broadcastBlockToAllLayers(L_block);
    stats.tsv_transfer_cycles += broadcast_cycles;
    
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
    
    for (uint32_t k = 0; k < b; k++) {
        for (uint32_t z = 0; z < config.num_layers; z++) {
            auto layer_rows = getLayerRows(z, b);
            for (uint32_t local_r = 0; local_r < layer_rows.size(); local_r++) {
                uint32_t j = layer_rows[local_r];
                if (j <= k) continue;
                
                for (uint32_t l = 0; l < b; l++) {
                    auto& pe = pe_array.getPE(z, local_r, l);
                    float a_jl = A_block[j][l];
                    float l_jk = L_block[j][k];
                    float a_kl = A_block[k][l];
                    pe.startMAC(a_jl, l_jk, a_kl, config.mac_latency);
                    stats.mac_operations++;
                }
            }
        }
        pe_array.waitUntilIdle();
        
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
    
    U_block = A_block;
    
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
}

void BlockLUSimulator3D::executeCase3_row(uint32_t block_k, uint32_t block_i) {
    // Case 3: Vertical L block update (backward substitution variant)
    // Row distribution: rows are independent - perfect for 3D parallelism
    
    uint64_t start_cycle = pe_array.current_cycle;
    uint32_t b = config.block_size;
    
    auto A_flat = memory.getBlockA(block_i, block_k);
    auto U_flat = memory.getBlockU(block_k, block_k);
    auto A_block = flatTo2D(A_flat, b);
    auto U_block = flatTo2D(U_flat, b);
    std::vector<std::vector<float>> L_block(b, std::vector<float>(b, 0.0f));
    
    uint64_t broadcast_cycles = broadcastBlockToAllLayers(U_block);
    stats.tsv_transfer_cycles += broadcast_cycles;
    
    for (uint64_t i = 0; i < broadcast_cycles; i++) {
        pe_array.tick();
    }
    
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
    
    for (uint32_t k = 0; k < b; k++) {
        float pivot = U_block[k][k];
        
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
        
        for (uint32_t z = 0; z < config.num_layers; z++) {
            auto layer_rows = getLayerRows(z, b);
            for (uint32_t local_r = 0; local_r < layer_rows.size(); local_r++) {
                uint32_t j = layer_rows[local_r];
                L_block[j][k] = pe_array.getPE(z, local_r, k).reg_result;
            }
        }
        
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
}

void BlockLUSimulator3D::executeCase4_row(uint32_t block_k, uint32_t block_i, uint32_t block_j) {
    // Case 4: Trailing matrix Schur complement update
    // Row distribution: rows of A are distributed across layers
    
    uint64_t start_cycle = pe_array.current_cycle;
    uint32_t b = config.block_size;
    
    auto A_flat = memory.getBlockA(block_i, block_j);
    auto L_flat = memory.getBlockL(block_i, block_k);
    auto U_flat = memory.getBlockU(block_k, block_j);
    auto A_block = flatTo2D(A_flat, b);
    auto L_block = flatTo2D(L_flat, b);
    auto U_block = flatTo2D(U_flat, b);
    
    uint64_t broadcast_cycles = broadcastBlockToAllLayers(U_block);
    stats.tsv_transfer_cycles += broadcast_cycles;
    
    for (uint64_t i = 0; i < broadcast_cycles; i++) {
        pe_array.tick();
    }
    
    for (uint32_t z = 0; z < config.num_layers; z++) {
        auto layer_rows = getLayerRows(z, b);
        for (uint32_t local_r = 0; local_r < layer_rows.size(); local_r++) {
            uint32_t global_r = layer_rows[local_r];
            for (uint32_t c = 0; c < b && c < config.pe_array_size; c++) {
                auto& pe = pe_array.getPE(z, local_r, c);
                pe.reg_a = A_block[global_r][c];
                pe.reg_b = L_block[global_r][c];
                pe.startLoad(config.mem_load_delay);
            }
        }
    }
    pe_array.waitUntilIdle();
    
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
}

// =============================================================================
// MAIN EXECUTION
// =============================================================================

void BlockLUSimulator3D::run() {
    if (config.verbose) {
        std::cout << "\n=== Starting 3D Block LU Decomposition ===\n";
        std::cout << "Layers: " << config.num_layers << "\n";
        std::cout << "Mapping: Block distribution (inter-block parallelism)\n\n";
    }
    
    // Reset state
    pe_array.reset();
    memory.initializeRandom(42);
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
        
        // Case 1: Diagonal block LU (Layer 0 only)
        executeCase1_block(k);
        
        // Case 2: All horizontal U blocks in parallel
        executeCase2_block(k, num_blocks);
        
        // Case 3: All vertical L blocks in parallel
        executeCase3_block(k, num_blocks);
        
        // Case 4: All trailing blocks in parallel
        if (config.pipeline_enabled) {
            if (config.pipeline_version == 1) {
                executeCase4_block_pipelined_v1(k, num_blocks);
            } else {
                executeCase4_block_pipelined(k, num_blocks);  // v2 (default)
            }
        } else {
            executeCase4_block(k, num_blocks);
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
