#include "block_lu_simulator.hpp"
#include <iostream>
#include <iomanip>

BlockLUSimulator::BlockLUSimulator(const SimConfig& cfg)
    : config(cfg)
    , pe_array(cfg)
    , memory(cfg)
    , stats()
{}

void BlockLUSimulator::initializeRandom(uint32_t seed) {
    memory.initializeRandom(seed);
    A_original = memory.A;  // Save original for verification
}

void BlockLUSimulator::initializeFromArray(const std::vector<float>& data) {
    memory.initializeFromArray(data);
    A_original = memory.A;  // Save original for verification
}

void BlockLUSimulator::log(const std::string& message) const {
    if (config.verbose) {
        std::cout << "[Cycle " << std::setw(6) << pe_array.current_cycle << "] " << message << "\n";
    }
}

void BlockLUSimulator::logCycleState() const {
    if (config.verbose) {
        pe_array.printState();
    }
}

uint64_t BlockLUSimulator::loadBlockToPEs(const std::vector<float>& block) {
    uint32_t b = config.block_size;
    uint64_t cycles = 0;
    
    // All PEs start loading simultaneously (parallel load)
    for (uint32_t i = 0; i < b; i++) {
        for (uint32_t j = 0; j < b; j++) {
            pe_array.getPE(i, j).reg_a = block[i * b + j];
            pe_array.getPE(i, j).startLoad(pe_array.getMemoryLoadDelay(i, j));
        }
    }
    
    // Wait for all loads to complete
    cycles = pe_array.waitUntilIdle();
    stats.memory_load_cycles += cycles;
    
    return cycles;
}

uint64_t BlockLUSimulator::writePEsToBlock(std::vector<float>& block) {
    uint32_t b = config.block_size;
    uint64_t cycles = 0;
    
    // Collect results from PEs
    for (uint32_t i = 0; i < b; i++) {
        for (uint32_t j = 0; j < b; j++) {
            block[i * b + j] = pe_array.getPE(i, j).reg_result;
            pe_array.getPE(i, j).startWrite(pe_array.getMemoryWriteDelay(i, j));
        }
    }
    
    // Wait for all writes to complete
    cycles = pe_array.waitUntilIdle();
    stats.memory_write_cycles += cycles;
    
    return cycles;
}

/**
 * Case 1: Diagonal block LU decomposition (lu_nopivot_basic)
 * 
 * Algorithm (from MATLAB):
 * for k = 1:b-1
 *     for j = k+1:b
 *         L(j,k) = A(j,k) / A(k,k)           <- Division
 *         for l = k:b
 *             A(j,l) = A(j,l) - L(j,k) * A(k,l)  <- MAC
 *         end
 *     end
 * end
 * U = A
 * 
 * This has sequential dependencies on the pivot element A(k,k)
 */
uint64_t BlockLUSimulator::executeCase1(uint32_t block_k) {
    log("=== Case 1: Diagonal LU on block (" + std::to_string(block_k) + "," + std::to_string(block_k) + ") ===");
    
    uint32_t b = config.block_size;
    uint64_t case_cycles = 0;
    
    // Load diagonal block A into PE array
    std::vector<float> A_block = memory.getBlockA(block_k, block_k);
    case_cycles += loadBlockToPEs(A_block);
    
    // Local copies for computation
    std::vector<float> L_block(b * b, 0.0f);
    std::vector<float> U_block(b * b, 0.0f);
    
    // Initialize L as identity for this block
    for (uint32_t i = 0; i < b; i++) {
        L_block[i * b + i] = 1.0f;
    }
    
    // Execute LU decomposition with cycle counting
    for (uint32_t k = 0; k < b - 1; k++) {
        log("  Pivot step k=" + std::to_string(k));
        
        for (uint32_t j = k + 1; j < b; j++) {
            // Division: L(j,k) = A(j,k) / A(k,k)
            float a_jk = A_block[j * b + k];
            float a_kk = A_block[k * b + k];
            
            pe_array.getPE(j, k).startDIV(a_jk, a_kk, config.div_latency);
            stats.div_operations++;
        }
        
        // Wait for all divisions to complete
        case_cycles += pe_array.waitUntilIdle();
        
        // Store L values
        for (uint32_t j = k + 1; j < b; j++) {
            L_block[j * b + k] = pe_array.getPE(j, k).getResult();
        }
        
        // MAC operations: A(j,l) = A(j,l) - L(j,k) * A(k,l) for all j > k, l >= k
        for (uint32_t j = k + 1; j < b; j++) {
            float l_jk = L_block[j * b + k];
            
            for (uint32_t l = k; l < b; l++) {
                float a_jl = A_block[j * b + l];
                float a_kl = A_block[k * b + l];
                
                pe_array.getPE(j, l).startMAC(a_jl, l_jk, a_kl, config.mac_latency);
                stats.mac_operations++;
            }
        }
        
        // Wait for all MACs to complete
        case_cycles += pe_array.waitUntilIdle();
        
        // Update A_block with results
        for (uint32_t j = k + 1; j < b; j++) {
            for (uint32_t l = k; l < b; l++) {
                A_block[j * b + l] = pe_array.getPE(j, l).getResult();
            }
        }
        
        if (config.verbose) {
            logCycleState();
        }
    }
    
    // U = final A_block (upper triangular part)
    U_block = A_block;
    
    // Write L and U blocks back to memory
    std::vector<float> L_result(b * b), U_result(b * b);
    for (uint32_t i = 0; i < b; i++) {
        for (uint32_t j = 0; j < b; j++) {
            pe_array.getPE(i, j).reg_result = L_block[i * b + j];
        }
    }
    case_cycles += writePEsToBlock(L_result);
    memory.setBlockL(block_k, block_k, L_block);
    
    for (uint32_t i = 0; i < b; i++) {
        for (uint32_t j = 0; j < b; j++) {
            pe_array.getPE(i, j).reg_result = U_block[i * b + j];
        }
    }
    case_cycles += writePEsToBlock(U_result);
    memory.setBlockU(block_k, block_k, U_block);
    
    stats.case1_cycles += case_cycles;
    log("  Case 1 complete: " + std::to_string(case_cycles) + " cycles");
    
    return case_cycles;
}

/**
 * Case 2: Horizontal block update (compute U blocks right of diagonal)
 * 
 * Given L from diagonal block, solve for U: L * U = A
 * This is forward substitution applied block-wise
 * 
 * Algorithm (from MATLAB case2 function):
 * for k = 1:b
 *     for j = k+1:b
 *         for l = 1:b
 *             A(j,l) = A(j,l) - L(j,k) * A(k,l)
 *         end
 *     end
 * end
 * U = A
 */
uint64_t BlockLUSimulator::executeCase2(uint32_t block_k, uint32_t block_j) {
    log("=== Case 2: Horizontal U update on block (" + std::to_string(block_k) + "," + std::to_string(block_j) + ") ===");
    
    uint32_t b = config.block_size;
    uint64_t case_cycles = 0;
    
    // Load A block and L block (from diagonal)
    std::vector<float> A_block = memory.getBlockA(block_k, block_j);
    std::vector<float> L_block = memory.getBlockL(block_k, block_k);
    
    case_cycles += loadBlockToPEs(A_block);
    
    // Execute forward substitution
    for (uint32_t k = 0; k < b; k++) {
        for (uint32_t j = k + 1; j < b; j++) {
            float l_jk = L_block[j * b + k];
            
            // All columns can be updated in parallel
            for (uint32_t l = 0; l < b; l++) {
                float a_jl = A_block[j * b + l];
                float a_kl = A_block[k * b + l];
                
                pe_array.getPE(j, l).startMAC(a_jl, l_jk, a_kl, config.mac_latency);
                stats.mac_operations++;
            }
        }
        
        // Wait for this k-iteration to complete
        case_cycles += pe_array.waitUntilIdle();
        
        // Update A_block
        for (uint32_t j = k + 1; j < b; j++) {
            for (uint32_t l = 0; l < b; l++) {
                A_block[j * b + l] = pe_array.getPE(j, l).getResult();
            }
        }
    }
    
    // Write U block to memory
    for (uint32_t i = 0; i < b; i++) {
        for (uint32_t j = 0; j < b; j++) {
            pe_array.getPE(i, j).reg_result = A_block[i * b + j];
        }
    }
    std::vector<float> U_result(b * b);
    case_cycles += writePEsToBlock(U_result);
    memory.setBlockU(block_k, block_j, A_block);
    
    stats.case2_cycles += case_cycles;
    log("  Case 2 complete: " + std::to_string(case_cycles) + " cycles");
    
    return case_cycles;
}

/**
 * Case 3: Vertical block update (compute L blocks below diagonal)
 * 
 * Given U from diagonal block, solve for L: L * U = A
 * 
 * Algorithm (from MATLAB case3 function):
 * for k = 1:b
 *     for j = 1:b
 *         L(j,k) = A(j,k) / U(k,k)           <- Division (scalar / scalar)
 *         for l = k+1:b
 *             A(j,l) = A(j,l) - L(j,k) * U(k,l)  <- Scalar broadcasts to row
 *         end
 *     end
 * end
 */
uint64_t BlockLUSimulator::executeCase3(uint32_t block_i, uint32_t block_k) {
    log("=== Case 3: Vertical L update on block (" + std::to_string(block_i) + "," + std::to_string(block_k) + ") ===");
    
    uint32_t b = config.block_size;
    uint64_t case_cycles = 0;
    
    // Load A block and U block (from diagonal)
    std::vector<float> A_block = memory.getBlockA(block_i, block_k);
    std::vector<float> U_block = memory.getBlockU(block_k, block_k);
    std::vector<float> L_block(b * b, 0.0f);
    
    case_cycles += loadBlockToPEs(A_block);
    
    // Execute backward substitution variant
    for (uint32_t k = 0; k < b; k++) {
        // Division: L(j,k) = A(j,k) / U(k,k) for all j
        // All rows can compute division in parallel
        float u_kk = U_block[k * b + k];
        
        for (uint32_t j = 0; j < b; j++) {
            float a_jk = A_block[j * b + k];
            pe_array.getPE(j, k).startDIV(a_jk, u_kk, config.div_latency);
            stats.div_operations++;
        }
        
        case_cycles += pe_array.waitUntilIdle();
        
        // Store L values
        for (uint32_t j = 0; j < b; j++) {
            L_block[j * b + k] = pe_array.getPE(j, k).getResult();
        }
        
        // MAC: A(j,l) = A(j,l) - L(j,k) * U(k,l) for all j, l > k
        // Scalar L(j,k) broadcasts to multiply row U(k,:)
        for (uint32_t j = 0; j < b; j++) {
            float l_jk = L_block[j * b + k];
            
            for (uint32_t l = k + 1; l < b; l++) {
                float a_jl = A_block[j * b + l];
                float u_kl = U_block[k * b + l];
                
                pe_array.getPE(j, l).startMAC(a_jl, l_jk, u_kl, config.mac_latency);
                stats.mac_operations++;
            }
        }
        
        case_cycles += pe_array.waitUntilIdle();
        
        // Update A_block
        for (uint32_t j = 0; j < b; j++) {
            for (uint32_t l = k + 1; l < b; l++) {
                A_block[j * b + l] = pe_array.getPE(j, l).getResult();
            }
        }
    }
    
    // Write L block to memory
    for (uint32_t i = 0; i < b; i++) {
        for (uint32_t j = 0; j < b; j++) {
            pe_array.getPE(i, j).reg_result = L_block[i * b + j];
        }
    }
    std::vector<float> L_result(b * b);
    case_cycles += writePEsToBlock(L_result);
    memory.setBlockL(block_i, block_k, L_block);
    
    stats.case3_cycles += case_cycles;
    log("  Case 3 complete: " + std::to_string(case_cycles) + " cycles");
    
    return case_cycles;
}

/**
 * Case 4: Trailing matrix update (Schur complement)
 * 
 * Algorithm (from MATLAB case4 function):
 * A_trail = A_trail - L_sub * U_sub
 * 
 * for k = 1:b
 *     for j = 1:b
 *         for l = 1:b
 *             A_trail(j,l) = A_trail(j,l) - L_sub(j,k) * U_sub(k,l)
 *         end
 *     end
 * end
 * 
 * This is essentially a matrix multiplication (outer product accumulation)
 */
uint64_t BlockLUSimulator::executeCase4(uint32_t block_i, uint32_t block_j, uint32_t block_k) {
    log("=== Case 4: Trailing update on block (" + std::to_string(block_i) + "," + std::to_string(block_j) + ") with k=" + std::to_string(block_k) + " ===");
    
    uint32_t b = config.block_size;
    uint64_t case_cycles = 0;
    
    // Load blocks
    std::vector<float> A_block = memory.getBlockA(block_i, block_j);
    std::vector<float> L_block = memory.getBlockL(block_i, block_k);
    std::vector<float> U_block = memory.getBlockU(block_k, block_j);
    
    case_cycles += loadBlockToPEs(A_block);
    
    // Execute outer product accumulation: A = A - L * U
    for (uint32_t k = 0; k < b; k++) {
        // All (j,l) pairs can be computed in parallel for a given k
        for (uint32_t j = 0; j < b; j++) {
            float l_jk = L_block[j * b + k];
            
            for (uint32_t l = 0; l < b; l++) {
                float a_jl = A_block[j * b + l];
                float u_kl = U_block[k * b + l];
                
                pe_array.getPE(j, l).startMAC(a_jl, l_jk, u_kl, config.mac_latency);
                stats.mac_operations++;
            }
        }
        
        case_cycles += pe_array.waitUntilIdle();
        
        // Update A_block
        for (uint32_t j = 0; j < b; j++) {
            for (uint32_t l = 0; l < b; l++) {
                A_block[j * b + l] = pe_array.getPE(j, l).getResult();
            }
        }
    }
    
    // Write updated A block back to memory
    for (uint32_t i = 0; i < b; i++) {
        for (uint32_t j = 0; j < b; j++) {
            pe_array.getPE(i, j).reg_result = A_block[i * b + j];
        }
    }
    std::vector<float> A_result(b * b);
    case_cycles += writePEsToBlock(A_result);
    memory.setBlockA(block_i, block_j, A_block);
    
    stats.case4_cycles += case_cycles;
    log("  Case 4 complete: " + std::to_string(case_cycles) + " cycles");
    
    return case_cycles;
}

/**
 * Main simulation loop implementing block LU decomposition
 * 
 * From MATLAB algorithm:
 * for k = 1:P
 *     for i = k:P
 *         for j = k:P
 *             if (i == k && j == k)  -> Case 1: diagonal
 *             if (i == k && j > k)   -> Case 2: horizontal
 *             if (i > k && j == k)   -> Case 3: vertical
 *             if (i > k && j > k)    -> Case 4: trailing
 *         end
 *     end
 * end
 */
void BlockLUSimulator::run() {
    uint32_t P = config.getNumBlocks();  // Number of blocks along each dimension
    
    std::cout << "Starting Block LU Decomposition Simulation\n";
    std::cout << "Matrix size: " << config.matrix_size << "x" << config.matrix_size << "\n";
    std::cout << "Block size: " << config.block_size << "\n";
    std::cout << "Number of blocks: " << P << "x" << P << "\n";
    std::cout << "PE array size: " << config.pe_array_size << "x" << config.pe_array_size << "\n\n";
    
    pe_array.reset();
    stats = SimStats();
    
    for (uint32_t k = 0; k < P; k++) {
        log("======== Block iteration k=" + std::to_string(k) + " ========");
        
        for (uint32_t i = k; i < P; i++) {
            for (uint32_t j = k; j < P; j++) {
                if (i == k && j == k) {
                    // Case 1: Diagonal block LU decomposition
                    stats.total_cycles += executeCase1(k);
                }
                else if (i == k && j > k) {
                    // Case 2: Horizontal U block update
                    stats.total_cycles += executeCase2(k, j);
                }
                else if (i > k && j == k) {
                    // Case 3: Vertical L block update
                    stats.total_cycles += executeCase3(i, k);
                }
                else if (i > k && j > k) {
                    // Case 4: Trailing matrix update
                    stats.total_cycles += executeCase4(i, j, k);
                }
            }
        }
    }
    
    // Collect final PE utilization statistics
    pe_array.collectStats(stats);
    
    std::cout << "Simulation complete.\n\n";
}

bool BlockLUSimulator::verify(float tolerance) const {
    return memory.verify(A_original, tolerance);
}

void BlockLUSimulator::printResults() const {
    config.print();
    stats.print();
    
    if (config.matrix_size <= 8) {
        std::cout << "\nMatrix Results:\n";
        memory.printL();
        memory.printU();
    }
}

// Implementation of print functions
void SimConfig::print() const {
    std::cout << "=== Simulation Configuration ===\n";
    std::cout << "Matrix size:       " << matrix_size << "x" << matrix_size << "\n";
    std::cout << "PE array size:     " << pe_array_size << "x" << pe_array_size << "\n";
    std::cout << "Block size:        " << block_size << "\n";
    std::cout << "Number of blocks:  " << getNumBlocks() << "x" << getNumBlocks() << "\n";
    std::cout << "MAC latency:       " << mac_latency << " cycles\n";
    std::cout << "DIV latency:       " << div_latency << " cycles\n";
    std::cout << "Memory load delay: " << mem_load_delay << " cycles\n";
    std::cout << "Memory write delay:" << mem_write_delay << " cycles\n";
    std::cout << "Verbose mode:      " << (verbose ? "ON" : "OFF") << "\n\n";
}

void SimStats::print() const {
    std::cout << "=== Simulation Statistics ===\n";
    std::cout << "Total cycles:           " << total_cycles << "\n\n";
    
    std::cout << "Cycles by case:\n";
    std::cout << "  Case 1 (Diagonal):    " << case1_cycles << " (" 
              << std::fixed << std::setprecision(1) 
              << (100.0 * case1_cycles / total_cycles) << "%)\n";
    std::cout << "  Case 2 (Horizontal):  " << case2_cycles << " (" 
              << (100.0 * case2_cycles / total_cycles) << "%)\n";
    std::cout << "  Case 3 (Vertical):    " << case3_cycles << " (" 
              << (100.0 * case3_cycles / total_cycles) << "%)\n";
    std::cout << "  Case 4 (Trailing):    " << case4_cycles << " (" 
              << (100.0 * case4_cycles / total_cycles) << "%)\n\n";
    
    std::cout << "Memory cycles:\n";
    std::cout << "  Load cycles:          " << memory_load_cycles << "\n";
    std::cout << "  Write cycles:         " << memory_write_cycles << "\n\n";
    
    std::cout << "Operations:\n";
    std::cout << "  MAC operations:       " << mac_operations << "\n";
    std::cout << "  DIV operations:       " << div_operations << "\n\n";
    
    std::cout << "PE Utilization:\n";
    std::cout << "  Active PE-cycles:     " << total_pe_active_cycles << "\n";
    std::cout << "  Total PE-cycles:      " << total_pe_possible_cycles << "\n";
    std::cout << "  Utilization:          " << std::fixed << std::setprecision(2) 
              << getPEUtilization() << "%\n\n";
}
