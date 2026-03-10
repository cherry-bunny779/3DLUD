#include "pipelined_case4.hpp"
#include <iostream>
#include <algorithm>
#include <cassert>
#include <iomanip>

PipelinedCase4Processor::PipelinedCase4Processor(const Config& config)
    : config_(config), b_(config.block_size), active_compute_buffer_(0), 
      active_store_buffer_(0) {
    
    initializePEStates();
    
    // Initialize double buffers
    uint32_t block_elements = b_ * b_;
    for (int i = 0; i < 2; i++) {
        U_buffer_[i].resize(block_elements, 0.0f);
        A_buffer_[i].resize(block_elements, 0.0f);
        A_result_[i].resize(block_elements, 0.0f);
    }
    
    stats_.reset();
}

void PipelinedCase4Processor::initializePEStates() {
    pe_states_.resize(b_ * b_);
    for (uint32_t r = 0; r < b_; r++) {
        for (uint32_t c = 0; c < b_; c++) {
            PEStatePipeline& pe = getPE(r, c);
            pe.row = r;
            pe.col = c;
            pe.reset();
        }
    }
}

void PipelinedCase4Processor::resetPEStates() {
    for (auto& pe : pe_states_) {
        pe.reset();
    }
}

void PipelinedCase4Processor::scheduleMemoryOps(uint32_t num_blocks, bool shared_L) {
    scheduled_ops_.clear();
    
    // Load shared block (L or U) at cycle 0
    MemoryOp shared_load;
    shared_load.type = shared_L ? MemoryOp::LOAD_L : MemoryOp::LOAD_U;
    shared_load.block_idx = 0;  // Shared block
    shared_load.start_cycle = 0;
    shared_load.end_cycle = config_.mem_load_delay;
    shared_load.completed = false;
    scheduled_ops_.push_back(shared_load);
    
    if (config_.enable_overlap) {
        // Overlapped scheduling: start loading next block while computing current
        // Block 0: Load U/L and A at cycle 0, ready at mem_load_delay
        // Block j: Start load when block j-1 is ~halfway through compute
        
        for (uint32_t blk = 0; blk < num_blocks; blk++) {
            uint32_t block_compute_start = blk * getPipelineInterval();
            
            // Varying block (U if shared_L, else L)
            MemoryOp varying_load;
            varying_load.type = shared_L ? MemoryOp::LOAD_U : MemoryOp::LOAD_L;
            varying_load.block_idx = blk;
            
            if (blk == 0) {
                // First block loads immediately after shared block
                varying_load.start_cycle = 0;
            } else {
                // Start loading when previous block is partway through
                // Ideally ready by the time block j starts computing
                uint32_t desired_ready = block_compute_start + config_.mem_load_delay;
                varying_load.start_cycle = (desired_ready > config_.mem_load_delay) 
                    ? desired_ready - config_.mem_load_delay : 0;
            }
            varying_load.end_cycle = varying_load.start_cycle + config_.mem_load_delay;
            varying_load.completed = false;
            scheduled_ops_.push_back(varying_load);
            
            // Load A block
            MemoryOp a_load;
            a_load.type = MemoryOp::LOAD_A;
            a_load.block_idx = blk;
            a_load.start_cycle = varying_load.start_cycle;  // Load together
            a_load.end_cycle = a_load.start_cycle + config_.mem_load_delay;
            a_load.completed = false;
            scheduled_ops_.push_back(a_load);
            
            // Store A result
            // Block j completes at getBlockCompletionCycle(block_compute_start)
            // Can start storing immediately after completion
            uint32_t block_complete = getBlockCompletionCycle(block_compute_start);
            MemoryOp a_store;
            a_store.type = MemoryOp::STORE_A;
            a_store.block_idx = blk;
            a_store.start_cycle = block_complete;
            a_store.end_cycle = a_store.start_cycle + config_.mem_write_delay;
            a_store.completed = false;
            scheduled_ops_.push_back(a_store);
        }
    } else {
        // Non-overlapped: sequential load-compute-store for each block
        uint32_t current_cycle = config_.mem_load_delay;  // After shared block loads
        
        for (uint32_t blk = 0; blk < num_blocks; blk++) {
            // Load varying block
            MemoryOp varying_load;
            varying_load.type = shared_L ? MemoryOp::LOAD_U : MemoryOp::LOAD_L;
            varying_load.block_idx = blk;
            varying_load.start_cycle = current_cycle;
            varying_load.end_cycle = current_cycle + config_.mem_load_delay;
            scheduled_ops_.push_back(varying_load);
            
            // Load A
            MemoryOp a_load;
            a_load.type = MemoryOp::LOAD_A;
            a_load.block_idx = blk;
            a_load.start_cycle = current_cycle;
            a_load.end_cycle = current_cycle + config_.mem_load_delay;
            scheduled_ops_.push_back(a_load);
            
            current_cycle += config_.mem_load_delay;
            
            // Compute (use skewed timing)
            uint32_t compute_end = current_cycle + 3 * (b_ - 1) + 1;
            
            // Store
            MemoryOp a_store;
            a_store.type = MemoryOp::STORE_A;
            a_store.block_idx = blk;
            a_store.start_cycle = compute_end;
            a_store.end_cycle = compute_end + config_.mem_write_delay;
            scheduled_ops_.push_back(a_store);
            
            current_cycle = a_store.end_cycle;
        }
    }
}

uint32_t PipelinedCase4Processor::getMemoryReadyCycle(uint32_t block_idx) const {
    // Find when both varying block and A block are loaded
    uint32_t varying_ready = 0;
    uint32_t a_ready = 0;
    
    for (const auto& op : scheduled_ops_) {
        if (op.block_idx == block_idx) {
            if (op.type == MemoryOp::LOAD_U || op.type == MemoryOp::LOAD_L) {
                varying_ready = op.end_cycle;
            } else if (op.type == MemoryOp::LOAD_A) {
                a_ready = op.end_cycle;
            }
        }
    }
    
    return std::max(varying_ready, a_ready);
}

bool PipelinedCase4Processor::isBlockDataReady(uint32_t block_idx, uint64_t cycle) const {
    return cycle >= getMemoryReadyCycle(block_idx);
}

void PipelinedCase4Processor::computeBlockMAC(
    const std::vector<float>& L_block,
    const std::vector<float>& U_block,
    std::vector<float>& A_block) {
    
    // A -= L * U (standard matrix multiply-subtract)
    for (uint32_t i = 0; i < b_; i++) {
        for (uint32_t j = 0; j < b_; j++) {
            float sum = 0.0f;
            for (uint32_t k = 0; k < b_; k++) {
                sum += L_block[i * b_ + k] * U_block[k * b_ + j];
            }
            A_block[i * b_ + j] -= sum;
        }
    }
}

uint32_t PipelinedCase4Processor::countActivePEs(
    uint64_t cycle, 
    const std::vector<uint32_t>& block_start_cycles) const {
    
    uint32_t active_count = 0;
    
    for (uint32_t r = 0; r < b_; r++) {
        for (uint32_t c = 0; c < b_; c++) {
            uint32_t skew = getSkewDelay(r, c);
            
            // Check each block
            for (size_t blk = 0; blk < block_start_cycles.size(); blk++) {
                uint32_t pe_start = block_start_cycles[blk] + skew;
                uint32_t pe_end = pe_start + b_ * config_.mac_latency;
                
                if (cycle >= pe_start && cycle < pe_end) {
                    // Check if data is ready
                    if (isBlockDataReady(static_cast<uint32_t>(blk), cycle)) {
                        active_count++;
                        break;  // PE can only work on one block at a time
                    }
                }
            }
        }
    }
    
    return active_count;
}

bool PipelinedCase4Processor::allBlocksComplete(uint32_t num_blocks, uint64_t cycle) const {
    if (num_blocks == 0) return true;
    
    // Last block starts at (num_blocks - 1) * pipeline_interval
    // Last PE [b-1, b-1] completes that block at:
    uint32_t last_block_start = (num_blocks - 1) * getPipelineInterval();
    uint32_t last_completion = getBlockCompletionCycle(last_block_start);
    
    return cycle > last_completion;
}

void PipelinedCase4Processor::simulatePipelinedExecution(
    const std::vector<float>& shared_block,
    const std::vector<std::vector<float>>& varying_blocks,
    std::vector<std::vector<float>>& A_blocks,
    bool shared_L) {
    
    uint32_t num_blocks = static_cast<uint32_t>(varying_blocks.size());
    
    // Calculate block start cycles (pipelined)
    std::vector<uint32_t> block_start_cycles(num_blocks);
    for (uint32_t blk = 0; blk < num_blocks; blk++) {
        // In pipelined mode, blocks start every pipeline_interval cycles
        // But must wait for memory to be ready
        uint32_t ideal_start = blk * getPipelineInterval();
        uint32_t mem_ready = getMemoryReadyCycle(blk);
        block_start_cycles[blk] = std::max(ideal_start, mem_ready);
    }
    
    // Find first compute start (after initial memory loads)
    uint32_t first_compute_cycle = block_start_cycles[0];
    
    // Find last completion cycle
    uint32_t last_block_start = block_start_cycles[num_blocks - 1];
    uint32_t last_pe_completion = getBlockCompletionCycle(last_block_start);
    
    // Find when all stores complete
    uint32_t last_store_complete = 0;
    for (const auto& op : scheduled_ops_) {
        if (op.type == MemoryOp::STORE_A) {
            last_store_complete = std::max(last_store_complete, op.end_cycle);
        }
    }
    
    uint64_t total_cycles = std::max(last_pe_completion, last_store_complete);
    
    // Perform actual computation (for correctness)
    for (uint32_t blk = 0; blk < num_blocks; blk++) {
        if (shared_L) {
            computeBlockMAC(shared_block, varying_blocks[blk], A_blocks[blk]);
        } else {
            computeBlockMAC(varying_blocks[blk], shared_block, A_blocks[blk]);
        }
    }
    
    // Count MAC operations
    stats_.mac_operations += static_cast<uint64_t>(num_blocks) * b_ * b_ * b_;
    stats_.blocks_processed += num_blocks;
    
    // Cycle-by-cycle PE activity tracking
    uint64_t pe_active_total = 0;
    uint64_t compute_cycles = 0;
    
    for (uint64_t cycle = first_compute_cycle; cycle <= last_pe_completion; cycle++) {
        uint32_t active = countActivePEs(cycle, block_start_cycles);
        pe_active_total += active;
        if (active > 0) {
            compute_cycles++;
        }
    }
    
    // Memory stalls: cycles where PEs could compute but data not ready
    uint64_t memory_stalls = 0;
    for (uint32_t blk = 0; blk < num_blocks; blk++) {
        uint32_t ideal_start = blk * getPipelineInterval();
        uint32_t actual_start = block_start_cycles[blk];
        if (actual_start > ideal_start) {
            memory_stalls += (actual_start - ideal_start);
        }
    }
    
    // Pipeline fill/drain
    uint32_t first_pe_start = block_start_cycles[0];  // PE[0,0]
    uint32_t last_pe_start_block0 = first_pe_start + getSkewDelay(b_ - 1, b_ - 1);
    stats_.pipeline_fill_cycles = last_pe_start_block0 - first_pe_start;
    
    uint32_t first_pe_done_last = block_start_cycles[num_blocks - 1] + b_ * config_.mac_latency;
    uint32_t last_pe_done_last = last_pe_completion;
    stats_.pipeline_drain_cycles = last_pe_done_last - first_pe_done_last;
    
    stats_.total_cycles += total_cycles;
    stats_.compute_cycles += compute_cycles;
    stats_.memory_stall_cycles += memory_stalls;
    stats_.pe_active_cycles += pe_active_total;
    stats_.pe_possible_cycles += total_cycles * b_ * b_;
    
    if (config_.verbose) {
        std::cout << "Pipelined execution: " << num_blocks << " blocks\n";
        std::cout << "  First compute cycle: " << first_compute_cycle << "\n";
        std::cout << "  Last PE completion: " << last_pe_completion << "\n";
        std::cout << "  Total cycles: " << total_cycles << "\n";
        std::cout << "  PE utilization: " << std::fixed << std::setprecision(1)
                  << (100.0 * pe_active_total / (total_cycles * b_ * b_)) << "%\n";
    }
}

uint64_t PipelinedCase4Processor::processBlockRow(
    const std::vector<float>& L_block,
    const std::vector<std::vector<float>>& U_blocks,
    std::vector<std::vector<float>>& A_blocks) {
    
    uint32_t num_blocks = static_cast<uint32_t>(U_blocks.size());
    if (num_blocks == 0) return 0;
    
    assert(U_blocks.size() == A_blocks.size());
    assert(L_block.size() == b_ * b_);
    for (const auto& blk : U_blocks) {
        assert(blk.size() == b_ * b_);
    }
    for (const auto& blk : A_blocks) {
        assert(blk.size() == b_ * b_);
    }
    
    resetPEStates();
    scheduleMemoryOps(num_blocks, true /* shared_L */);
    simulatePipelinedExecution(L_block, U_blocks, A_blocks, true);
    
    return stats_.total_cycles;
}

uint64_t PipelinedCase4Processor::processBlockCol(
    const std::vector<std::vector<float>>& L_blocks,
    const std::vector<float>& U_block,
    std::vector<std::vector<float>>& A_blocks) {
    
    uint32_t num_blocks = static_cast<uint32_t>(L_blocks.size());
    if (num_blocks == 0) return 0;
    
    assert(L_blocks.size() == A_blocks.size());
    assert(U_block.size() == b_ * b_);
    for (const auto& blk : L_blocks) {
        assert(blk.size() == b_ * b_);
    }
    for (const auto& blk : A_blocks) {
        assert(blk.size() == b_ * b_);
    }
    
    resetPEStates();
    scheduleMemoryOps(num_blocks, false /* shared_L = false means shared_U */);
    simulatePipelinedExecution(U_block, L_blocks, A_blocks, false);
    
    return stats_.total_cycles;
}

uint64_t PipelinedCase4Processor::processSequential(
    const std::vector<std::vector<float>>& L_blocks,
    const std::vector<std::vector<float>>& U_blocks,
    std::vector<std::vector<float>>& A_blocks) {
    
    uint32_t num_blocks = static_cast<uint32_t>(L_blocks.size());
    if (num_blocks == 0) return 0;
    
    assert(L_blocks.size() == U_blocks.size());
    assert(L_blocks.size() == A_blocks.size());
    
    stats_.reset();
    
    uint64_t total_cycles = 0;
    
    for (uint32_t blk = 0; blk < num_blocks; blk++) {
        // Load L, U, A
        uint64_t load_cycles = config_.mem_load_delay;  // Assume parallel load
        
        // Compute with skewing
        // Single block: PE[0,0] starts at 0, PE[b-1,b-1] completes at 3(b-1)
        uint64_t compute_cycles = 3 * (b_ - 1) + 1;
        
        // Store
        uint64_t store_cycles = config_.mem_write_delay;
        
        total_cycles += load_cycles + compute_cycles + store_cycles;
        
        // Perform actual computation
        computeBlockMAC(L_blocks[blk], U_blocks[blk], A_blocks[blk]);
        stats_.mac_operations += b_ * b_ * b_;
    }
    
    stats_.total_cycles = total_cycles;
    stats_.compute_cycles = num_blocks * (3 * (b_ - 1) + 1);
    stats_.blocks_processed = num_blocks;
    stats_.pe_active_cycles = num_blocks * b_ * b_ * b_;  // All PEs active for b cycles each
    stats_.pe_possible_cycles = total_cycles * b_ * b_;
    
    return total_cycles;
}
