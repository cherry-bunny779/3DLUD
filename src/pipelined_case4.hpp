#ifndef PIPELINED_CASE4_HPP
#define PIPELINED_CASE4_HPP

#include <vector>
#include <cstdint>
#include <queue>
#include <functional>

/**
 * PE state for pipelined systolic execution.
 * Tracks per-PE progress through block computations with skewed dataflow.
 */
struct PEStatePipeline {
    uint32_t row, col;              // Position in PE array
    
    // Current work tracking
    int32_t current_block_idx;      // Which block in batch (-1 = idle)
    uint32_t k_iteration;           // 0 to b-1 within current block
    float accumulator;              // Running sum for current output
    
    // Timing (relative to batch start)
    uint32_t block_start_cycle;     // When this PE starts current block
    uint32_t block_end_cycle;       // When this PE completes current block
    
    // Per-cycle state
    bool active_this_cycle;         // Did MAC this cycle
    float result;                   // Completed value
    bool result_ready;              // Output available
    
    PEStatePipeline() : row(0), col(0), current_block_idx(-1), k_iteration(0),
                        accumulator(0.0f), block_start_cycle(0), block_end_cycle(0),
                        active_this_cycle(false), result(0.0f), result_ready(false) {}
    
    void reset() {
        current_block_idx = -1;
        k_iteration = 0;
        accumulator = 0.0f;
        block_start_cycle = 0;
        block_end_cycle = 0;
        active_this_cycle = false;
        result = 0.0f;
        result_ready = false;
    }
};

/**
 * Memory operation for scheduling.
 */
struct MemoryOp {
    enum Type { LOAD_L, LOAD_U, LOAD_A, STORE_A };
    Type type;
    uint32_t block_idx;             // Which block this op is for
    uint32_t start_cycle;           // When op starts
    uint32_t end_cycle;             // When op completes (data available)
    bool completed;
    
    MemoryOp() : type(LOAD_U), block_idx(0), start_cycle(0), 
                 end_cycle(0), completed(false) {}
    
    MemoryOp(Type t, uint32_t idx, uint32_t start, uint32_t end)
        : type(t), block_idx(idx), start_cycle(start), 
          end_cycle(end), completed(false) {}
};

/**
 * Statistics for pipelined execution.
 */
struct PipelineStats {
    uint64_t total_cycles;
    uint64_t compute_cycles;        // Cycles with at least one PE active
    uint64_t memory_stall_cycles;   // Cycles waiting for memory
    uint64_t pipeline_fill_cycles;  // Initial fill latency (first PE start to last PE start)
    uint64_t pipeline_drain_cycles; // Drain latency (first PE done to last PE done)
    uint64_t mac_operations;
    uint64_t blocks_processed;
    
    // Utilization tracking
    uint64_t pe_active_cycles;      // Sum of (active cycles per PE)
    uint64_t pe_possible_cycles;    // Total PE-cycles available
    
    PipelineStats() : total_cycles(0), compute_cycles(0), memory_stall_cycles(0),
                      pipeline_fill_cycles(0), pipeline_drain_cycles(0),
                      mac_operations(0), blocks_processed(0),
                      pe_active_cycles(0), pe_possible_cycles(0) {}
    
    void reset() {
        total_cycles = 0;
        compute_cycles = 0;
        memory_stall_cycles = 0;
        pipeline_fill_cycles = 0;
        pipeline_drain_cycles = 0;
        mac_operations = 0;
        blocks_processed = 0;
        pe_active_cycles = 0;
        pe_possible_cycles = 0;
    }
    
    double getUtilization() const {
        if (pe_possible_cycles == 0) return 0.0;
        return 100.0 * pe_active_cycles / pe_possible_cycles;
    }
};

/**
 * Pipelined Case 4 (trailing update) processor.
 * 
 * Implements classic systolic dataflow with skewed inputs:
 * - PE[i,j] receives data at cycle (i + j) due to row/col skewing
 * - Blocks sharing L (same row) or U (same column) can be pipelined
 * - New block starts every b cycles (after previous block's data fully streamed)
 * - Memory operations overlapped with computation via double-buffering
 * 
 * Timing model (b=4 example, single block):
 *   PE[0,0]: starts @0, completes @3   (k=0,1,2,3)
 *   PE[0,3]: starts @3, completes @6   (skew = 3)
 *   PE[3,0]: starts @3, completes @6
 *   PE[3,3]: starts @6, completes @9   (last PE)
 *   Single block latency = 2*(b-1) + b = 3*b - 2 = 10 cycles
 * 
 * Pipelined blocks (N blocks):
 *   Block j starts at cycle j * b
 *   Total = (N-1)*b + (3*b - 2) = N*b + 2*b - 2
 */
class PipelinedCase4Processor {
public:
    struct Config {
        uint32_t block_size;        // b
        uint32_t mem_load_delay;    // Cycles to load one block
        uint32_t mem_write_delay;   // Cycles to write one block
        uint32_t mac_latency;       // Cycles per MAC (typically 1)
        bool enable_overlap;        // Double-buffer memory ops with compute
        bool verbose;               // Debug output
        
        Config() : block_size(4), mem_load_delay(5), mem_write_delay(5),
                   mac_latency(1), enable_overlap(true), verbose(false) {}
    };
    
    explicit PipelinedCase4Processor(const Config& config);
    
    /**
     * Process a row of trailing blocks (all share the same L block).
     * L is loaded once, U blocks streamed in pipeline fashion.
     * 
     * @param L_block     Shared L block (b x b, row-major)
     * @param U_blocks    Vector of U blocks, one per trailing block
     * @param A_blocks    Vector of A blocks, updated in place (A -= L*U)
     * @return            Total cycle count for this batch
     */
    uint64_t processBlockRow(
        const std::vector<float>& L_block,
        const std::vector<std::vector<float>>& U_blocks,
        std::vector<std::vector<float>>& A_blocks
    );
    
    /**
     * Process a column of trailing blocks (all share the same U block).
     * U is loaded once, L blocks streamed in pipeline fashion.
     */
    uint64_t processBlockCol(
        const std::vector<std::vector<float>>& L_blocks,
        const std::vector<float>& U_block,
        std::vector<std::vector<float>>& A_blocks
    );
    
    /**
     * Process blocks sequentially without pipelining.
     * Baseline for comparison.
     */
    uint64_t processSequential(
        const std::vector<std::vector<float>>& L_blocks,
        const std::vector<std::vector<float>>& U_blocks,
        std::vector<std::vector<float>>& A_blocks
    );
    
    const PipelineStats& getStats() const { return stats_; }
    void resetStats() { stats_.reset(); }
    const Config& getConfig() const { return config_; }

private:
    Config config_;
    uint32_t b_;                    // Block size shorthand
    std::vector<PEStatePipeline> pe_states_;  // b*b PE states (flattened)
    PipelineStats stats_;
    
    // Memory operation scheduling
    std::vector<MemoryOp> scheduled_ops_;
    
    // Double buffers for overlapped memory access
    std::vector<float> U_buffer_[2];
    std::vector<float> A_buffer_[2];
    std::vector<float> A_result_[2];  // Results waiting to be stored
    int active_compute_buffer_;        // Buffer being computed
    int active_store_buffer_;          // Buffer being stored
    
    // Access helpers
    PEStatePipeline& getPE(uint32_t row, uint32_t col) {
        return pe_states_[row * b_ + col];
    }
    const PEStatePipeline& getPE(uint32_t row, uint32_t col) const {
        return pe_states_[row * b_ + col];
    }
    
    /**
     * Skew delay for PE[row, col].
     * Classic systolic: data arrives at cycle (row + col).
     */
    uint32_t getSkewDelay(uint32_t row, uint32_t col) const {
        return row + col;
    }
    
    /**
     * Cycle when PE[row,col] completes a block that started at base_cycle.
     * = base_cycle + skew_delay + (b-1)*mac_latency + mac_latency
     * = base_cycle + row + col + b*mac_latency
     */
    uint32_t getCompletionCycle(uint32_t row, uint32_t col, uint32_t base_cycle) const {
        return base_cycle + getSkewDelay(row, col) + b_ * config_.mac_latency;
    }
    
    /**
     * Cycle when last PE completes for a block starting at base_cycle.
     */
    uint32_t getBlockCompletionCycle(uint32_t base_cycle) const {
        // Last PE is [b-1, b-1]
        return getCompletionCycle(b_ - 1, b_ - 1, base_cycle);
    }
    
    /**
     * Pipeline interval: cycles between starting consecutive blocks.
     */
    uint32_t getPipelineInterval() const {
        return b_ * config_.mac_latency;
    }
    
    // Initialization
    void initializePEStates();
    void resetPEStates();
    
    // Memory scheduling
    void scheduleMemoryOps(uint32_t num_blocks, bool shared_L);
    uint32_t getMemoryReadyCycle(uint32_t block_idx) const;
    bool isBlockDataReady(uint32_t block_idx, uint64_t cycle) const;
    
    // Cycle-accurate simulation
    void simulatePipelinedExecution(
        const std::vector<float>& shared_block,     // L (if shared_L) or U
        const std::vector<std::vector<float>>& varying_blocks,  // U or L blocks
        std::vector<std::vector<float>>& A_blocks,
        bool shared_L
    );
    
    void simulateCycle(
        uint64_t cycle,
        const std::vector<float>& shared_block,
        const std::vector<std::vector<float>>& varying_blocks,
        std::vector<std::vector<float>>& A_blocks,
        bool shared_L,
        const std::vector<uint32_t>& block_start_cycles
    );
    
    // Block computation (actual math for verification)
    void computeBlockMAC(
        const std::vector<float>& L_block,
        const std::vector<float>& U_block,
        std::vector<float>& A_block
    );
    
    // State queries
    bool allBlocksComplete(uint32_t num_blocks, uint64_t cycle) const;
    uint32_t countActivePEs(uint64_t cycle, const std::vector<uint32_t>& block_start_cycles) const;
};

#endif // PIPELINED_CASE4_HPP
