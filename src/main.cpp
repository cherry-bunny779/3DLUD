#include "block_lu_simulator.hpp"
#include "block_lu_simulator_3d.hpp"
#include "config_3d.hpp"
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include <cmath>

// Helper function to check if a number is a perfect square
bool isPerfectSquare(uint32_t n) {
    if (n == 0) return false;
    uint32_t root = static_cast<uint32_t>(std::sqrt(n));
    return root * root == n;
}

// Helper function to get integer square root
uint32_t intSqrt(uint32_t n) {
    return static_cast<uint32_t>(std::sqrt(n));
}

void printUsage(const char* program_name) {
    std::cout << "Usage: " << program_name << " [options]\n\n";
    std::cout << "Options:\n";
    std::cout << "  -n <size>      Matrix size (power of 2, default: 16)\n";
    std::cout << "  -p <size>      PE array size (power of 2, default: 4)\n";
    std::cout << "  -b <size>      Block size (power of 2, default: equals PE array size)\n";
    std::cout << "  -mac <cycles>  MAC latency in cycles (default: 1)\n";
    std::cout << "  -div <cycles>  Division latency in cycles (default: 10)\n";
    std::cout << "  -load <cycles> Memory load delay in cycles (default: 5)\n";
    std::cout << "  -write <cycles> Memory write delay in cycles (default: 5)\n";
    std::cout << "  -seed <value>  Random seed for matrix generation (default: 42)\n";
    std::cout << "  -v             Enable verbose mode (cycle-by-cycle trace)\n";
    std::cout << "  -pipeline      Enable pipelined Case 4 execution (row-wise skewed systolic)\n";
    std::cout << "  -h, --help     Show this help message\n\n";
    std::cout << "3D Options:\n";
    std::cout << "  -3d            Enable 3D simulation mode\n";
    std::cout << "  -z <layers>    Number of Z layers (default: 4)\n";
    std::cout << "  -tsv <cycles>  TSV latency per hop in cycles (default: 1)\n\n";
    std::cout << "Bandwidth Contention Model (Case 4 only, used in -compare-part):\n";
    std::cout << "  -contention <m>   Enable bandwidth contention model with mode m:\n";
    std::cout << "                      0 = ideal      (no slowdown, legacy behavior)\n";
    std::cout << "                      1 = realistic  (U broadcast free, L/A serialized)\n";
    std::cout << "                      2 = pessimistic (full per-region serialization)\n";
    std::cout << "                    If omitted, model is disabled (slowdown = 1.0x).\n\n";
    std::cout << "Comparison Modes:\n";
    std::cout << "  -compare          Run both 2D and 3D with same PE array, print comparison\n";
    std::cout << "  -compare-pipeline Compare pipelined vs non-pipelined 2D execution\n";
    std::cout << "  -compare-part     Compare 2D vs 2D-Partitioned vs 3D (same total PEs)\n";
    std::cout << "                    2D-Partitioned: p×p array split into 4 regions of (p/2)×(p/2)\n";
    std::cout << "                    Shows benefit of spatial parallelism and TSV overhead cost\n";
    std::cout << "  -equal-pe         Equal PE budget: 2D (p x p) vs 3D ((p/sqrt(Z)) x (p/sqrt(Z)) x Z)\n";
    std::cout << "                    Requires: Z is a perfect square and sqrt(Z) divides p\n\n";
    std::cout << "Examples:\n";
    std::cout << "  " << program_name << " -n 32 -p 8 -mac 1 -div 10 -v\n";
    std::cout << "  " << program_name << " -n 64 -p 8 -pipeline          # Enable pipelined trailing update\n";
    std::cout << "  " << program_name << " -n 64 -p 8 -compare-pipeline  # Compare pipelined vs non-pipelined\n";
    std::cout << "  " << program_name << " -n 64 -p 8 -compare-part      # 2D vs 2D-Partitioned vs 3D\n";
    std::cout << "  " << program_name << " -n 32 -p 8 -3d -z 4\n";
    std::cout << "  " << program_name << " -n 32 -p 8 -compare -z 4\n";
    std::cout << "  " << program_name << " -n 64 -p 8 -equal-pe -z 4     # 2D: 8x8=64 PEs, 3D: 4x4x4=64 PEs\n\n";
}

void printComparison(const SimStats& stats_2d, const SimStats3D& stats_3d, 
                     uint32_t num_layers) {
    std::cout << "\n";
    std::cout << "========================================================================\n";
    std::cout << "                    2D vs 3D Comparison Results                         \n";
    std::cout << "========================================================================\n";
    std::cout << "\n";
    
    // Total cycles
    double speedup = static_cast<double>(stats_2d.total_cycles) / 
                     std::max(stats_3d.total_cycles, (uint64_t)1);
    std::cout << "Total Cycles:\n";
    std::cout << "    2D:          " << std::setw(12) << stats_2d.total_cycles << "\n";
    std::cout << "    3D (" << num_layers << "L):     " << std::setw(12) << stats_3d.total_cycles << "\n";
    std::cout << "    Speedup:     " << std::fixed << std::setprecision(2) << std::setw(12) 
              << speedup << "x\n\n";
    
    // Case breakdown
    std::cout << "Case Breakdown (cycles):\n";
    std::cout << "                          2D            3D        Speedup\n";
    
    auto calc_speedup = [](uint64_t a, uint64_t b) {
        return static_cast<double>(a) / std::max(b, (uint64_t)1);
    };
    
    std::cout << "    Case 1 (diag):   " << std::setw(10) << stats_2d.case1_cycles 
              << "    " << std::setw(10) << stats_3d.case1_cycles 
              << "    " << std::setw(6) << calc_speedup(stats_2d.case1_cycles, stats_3d.case1_cycles) << "x\n";
    std::cout << "    Case 2 (horiz):  " << std::setw(10) << stats_2d.case2_cycles 
              << "    " << std::setw(10) << stats_3d.case2_cycles 
              << "    " << std::setw(6) << calc_speedup(stats_2d.case2_cycles, stats_3d.case2_cycles) << "x\n";
    std::cout << "    Case 3 (vert):   " << std::setw(10) << stats_2d.case3_cycles 
              << "    " << std::setw(10) << stats_3d.case3_cycles 
              << "    " << std::setw(6) << calc_speedup(stats_2d.case3_cycles, stats_3d.case3_cycles) << "x\n";
    std::cout << "    Case 4 (trail):  " << std::setw(10) << stats_2d.case4_cycles 
              << "    " << std::setw(10) << stats_3d.case4_cycles 
              << "    " << std::setw(6) << calc_speedup(stats_2d.case4_cycles, stats_3d.case4_cycles) << "x\n\n";
    
    // Utilization
    double util_2d = (stats_2d.total_pe_possible_cycles > 0) ?
        (100.0 * stats_2d.total_pe_active_cycles / stats_2d.total_pe_possible_cycles) : 0.0;
    double util_3d = (stats_3d.total_pe_possible_cycles > 0) ?
        (100.0 * stats_3d.total_pe_active_cycles / stats_3d.total_pe_possible_cycles) : 0.0;
    
    std::cout << "PE Utilization:\n";
    std::cout << "    2D:          " << std::setw(12) << std::setprecision(2) << util_2d << "%\n";
    std::cout << "    3D:          " << std::setw(12) << util_3d << "%\n\n";
    
    // Operations
    std::cout << "Operations:\n";
    std::cout << "    MACs:        " << std::setw(12) << stats_2d.mac_operations 
              << " (2D)  " << std::setw(12) << stats_3d.mac_operations << " (3D)\n";
    std::cout << "    DIVs:        " << std::setw(12) << stats_2d.div_operations 
              << " (2D)  " << std::setw(12) << stats_3d.div_operations << " (3D)\n\n";
    
    // TSV overhead
    std::cout << "3D-Specific Overhead:\n";
    std::cout << "    TSV Transfer Cycles: " << std::setw(10) << stats_3d.tsv_transfer_cycles << "\n";
    double tsv_overhead = (stats_3d.total_cycles > 0) ?
        (100.0 * stats_3d.tsv_transfer_cycles / stats_3d.total_cycles) : 0.0;
    std::cout << "    TSV Overhead:        " << std::setw(10) << std::setprecision(1) 
              << tsv_overhead << "%\n";
    
    std::cout << "\n========================================================================\n";
}

void printEqualPEComparison(const SimStats& stats_2d, const SimStats3D& stats_3d,
                            uint32_t pe_2d, uint32_t pe_3d, uint32_t block_2d, uint32_t block_3d,
                            uint32_t num_layers, uint32_t total_pes) {
    std::cout << "\n";
    std::cout << "========================================================================\n";
    std::cout << "           Equal PE Budget: 2D vs 3D Comparison Results                \n";
    std::cout << "========================================================================\n";
    std::cout << "\n";
    
    // PE configuration summary
    std::cout << "Silicon Budget (Total PEs): " << total_pes << "\n\n";
    std::cout << "PE Array Configuration:\n";
    std::cout << "    2D: " << pe_2d << " x " << pe_2d << " = " << (pe_2d * pe_2d) << " PEs\n";
    std::cout << "    3D: " << pe_3d << " x " << pe_3d << " x " << num_layers 
              << " = " << (pe_3d * pe_3d * num_layers) << " PEs\n\n";
    std::cout << "Block Size:\n";
    std::cout << "    2D: " << block_2d << " x " << block_2d << "\n";
    std::cout << "    3D: " << block_3d << " x " << block_3d << "\n\n";
    
    // Total cycles
    double speedup = static_cast<double>(stats_2d.total_cycles) / 
                     std::max(stats_3d.total_cycles, (uint64_t)1);
    std::cout << "Total Cycles:\n";
    std::cout << "    2D:          " << std::setw(12) << stats_2d.total_cycles << "\n";
    std::cout << "    3D (" << num_layers << "L):     " << std::setw(12) << stats_3d.total_cycles << "\n";
    std::cout << "    Speedup:     " << std::fixed << std::setprecision(2) << std::setw(12) 
              << speedup << "x";
    if (speedup > 1.0) {
        std::cout << " (3D is faster)";
    } else if (speedup < 1.0) {
        std::cout << " (2D is faster)";
    } else {
        std::cout << " (equal)";
    }
    std::cout << "\n\n";
    
    // Case breakdown
    std::cout << "Case Breakdown (cycles):\n";
    std::cout << "                          2D            3D        Speedup\n";
    
    auto calc_speedup = [](uint64_t a, uint64_t b) {
        return static_cast<double>(a) / std::max(b, (uint64_t)1);
    };
    
    std::cout << "    Case 1 (diag):   " << std::setw(10) << stats_2d.case1_cycles 
              << "    " << std::setw(10) << stats_3d.case1_cycles 
              << "    " << std::setw(6) << calc_speedup(stats_2d.case1_cycles, stats_3d.case1_cycles) << "x\n";
    std::cout << "    Case 2 (horiz):  " << std::setw(10) << stats_2d.case2_cycles 
              << "    " << std::setw(10) << stats_3d.case2_cycles 
              << "    " << std::setw(6) << calc_speedup(stats_2d.case2_cycles, stats_3d.case2_cycles) << "x\n";
    std::cout << "    Case 3 (vert):   " << std::setw(10) << stats_2d.case3_cycles 
              << "    " << std::setw(10) << stats_3d.case3_cycles 
              << "    " << std::setw(6) << calc_speedup(stats_2d.case3_cycles, stats_3d.case3_cycles) << "x\n";
    std::cout << "    Case 4 (trail):  " << std::setw(10) << stats_2d.case4_cycles 
              << "    " << std::setw(10) << stats_3d.case4_cycles 
              << "    " << std::setw(6) << calc_speedup(stats_2d.case4_cycles, stats_3d.case4_cycles) << "x\n\n";
    
    // Utilization
    double util_2d = (stats_2d.total_pe_possible_cycles > 0) ?
        (100.0 * stats_2d.total_pe_active_cycles / stats_2d.total_pe_possible_cycles) : 0.0;
    double util_3d = (stats_3d.total_pe_possible_cycles > 0) ?
        (100.0 * stats_3d.total_pe_active_cycles / stats_3d.total_pe_possible_cycles) : 0.0;
    
    std::cout << "PE Utilization:\n";
    std::cout << "    2D:          " << std::setw(12) << std::setprecision(2) << util_2d << "%\n";
    std::cout << "    3D:          " << std::setw(12) << util_3d << "%\n\n";
    
    // Operations
    std::cout << "Operations:\n";
    std::cout << "    MACs:        " << std::setw(12) << stats_2d.mac_operations 
              << " (2D)  " << std::setw(12) << stats_3d.mac_operations << " (3D)\n";
    std::cout << "    DIVs:        " << std::setw(12) << stats_2d.div_operations 
              << " (2D)  " << std::setw(12) << stats_3d.div_operations << " (3D)\n\n";
    
    // Number of blocks (affects parallelism potential)
    std::cout << "Algorithm Granularity:\n";
    std::cout << "    2D blocks: n/b = n/" << block_2d << " blocks per dimension\n";
    std::cout << "    3D blocks: n/b = n/" << block_3d << " blocks per dimension\n";
    std::cout << "    (Smaller blocks = more blocks = more parallelism potential)\n\n";
    
    // TSV overhead
    std::cout << "3D-Specific Overhead:\n";
    std::cout << "    TSV Transfer Cycles: " << std::setw(10) << stats_3d.tsv_transfer_cycles << "\n";
    double tsv_overhead = (stats_3d.total_cycles > 0) ?
        (100.0 * stats_3d.tsv_transfer_cycles / stats_3d.total_cycles) : 0.0;
    std::cout << "    TSV Overhead:        " << std::setw(10) << std::setprecision(1) 
              << tsv_overhead << "%\n";
    
    // Summary
    std::cout << "\n------------------------------------------------------------------------\n";
    std::cout << "CONCLUSION: With equal silicon budget (" << total_pes << " PEs),\n";
    if (speedup > 1.05) {
        std::cout << "    3D stacking is " << std::setprecision(2) << speedup << "x FASTER than 2D.\n";
        std::cout << "    The layer parallelism outweighs the smaller per-layer PE array.\n";
    } else if (speedup < 0.95) {
        std::cout << "    2D is " << std::setprecision(2) << (1.0/speedup) << "x FASTER than 3D stacking.\n";
        std::cout << "    The larger PE array outweighs the layer parallelism benefit.\n";
    } else {
        std::cout << "    2D and 3D have similar performance.\n";
    }
    std::cout << "========================================================================\n";
}

void printPipelineComparison(const SimStats& stats_seq, const SimStats& stats_pipe,
                              uint32_t matrix_size, uint32_t pe_size, uint32_t block_size) {
    std::cout << "\n";
    std::cout << "========================================================================\n";
    std::cout << "          Pipelined vs Non-Pipelined Comparison Results                \n";
    std::cout << "========================================================================\n";
    std::cout << "\n";
    
    std::cout << "Configuration:\n";
    std::cout << "    Matrix size:   " << matrix_size << " x " << matrix_size << "\n";
    std::cout << "    PE array:      " << pe_size << " x " << pe_size << "\n";
    std::cout << "    Block size:    " << block_size << " x " << block_size << "\n\n";
    
    // Total cycles
    double speedup = static_cast<double>(stats_seq.total_cycles) / 
                     std::max(stats_pipe.total_cycles, (uint64_t)1);
    std::cout << "Total Cycles:\n";
    std::cout << "    Sequential:  " << std::setw(12) << stats_seq.total_cycles << "\n";
    std::cout << "    Pipelined:   " << std::setw(12) << stats_pipe.total_cycles << "\n";
    std::cout << "    Speedup:     " << std::fixed << std::setprecision(2) << std::setw(12) 
              << speedup << "x\n\n";
    
    // Case breakdown
    std::cout << "Case Breakdown (cycles):\n";
    std::cout << "                       Sequential    Pipelined     Speedup\n";
    
    auto calc_speedup = [](uint64_t a, uint64_t b) {
        return static_cast<double>(a) / std::max(b, (uint64_t)1);
    };
    
    std::cout << "    Case 1 (diag):   " << std::setw(10) << stats_seq.case1_cycles 
              << "    " << std::setw(10) << stats_pipe.case1_cycles 
              << "    " << std::setw(6) << calc_speedup(stats_seq.case1_cycles, stats_pipe.case1_cycles) << "x\n";
    std::cout << "    Case 2 (horiz):  " << std::setw(10) << stats_seq.case2_cycles 
              << "    " << std::setw(10) << stats_pipe.case2_cycles 
              << "    " << std::setw(6) << calc_speedup(stats_seq.case2_cycles, stats_pipe.case2_cycles) << "x\n";
    std::cout << "    Case 3 (vert):   " << std::setw(10) << stats_seq.case3_cycles 
              << "    " << std::setw(10) << stats_pipe.case3_cycles 
              << "    " << std::setw(6) << calc_speedup(stats_seq.case3_cycles, stats_pipe.case3_cycles) << "x\n";
    std::cout << "    Case 4 (trail):  " << std::setw(10) << stats_seq.case4_cycles 
              << "    " << std::setw(10) << stats_pipe.case4_cycles 
              << "    " << std::setw(6) << calc_speedup(stats_seq.case4_cycles, stats_pipe.case4_cycles) << "x\n\n";
    
    // Case 4 improvement breakdown
    double case4_reduction = 100.0 * (1.0 - static_cast<double>(stats_pipe.case4_cycles) / 
                                      std::max(stats_seq.case4_cycles, (uint64_t)1));
    std::cout << "Case 4 Improvement:\n";
    std::cout << "    Cycle reduction: " << std::setw(10) << std::setprecision(1) 
              << case4_reduction << "%\n";
    std::cout << "    (Pipelining overlaps memory loads with computation)\n\n";
    
    // Memory cycles
    std::cout << "Memory Cycles:\n";
    std::cout << "    Sequential Load:  " << std::setw(10) << stats_seq.memory_load_cycles << "\n";
    std::cout << "    Pipelined Load:   " << std::setw(10) << stats_pipe.memory_load_cycles << "\n";
    std::cout << "    Sequential Write: " << std::setw(10) << stats_seq.memory_write_cycles << "\n";
    std::cout << "    Pipelined Write:  " << std::setw(10) << stats_pipe.memory_write_cycles << "\n\n";
    
    // Utilization
    double util_seq = (stats_seq.total_pe_possible_cycles > 0) ?
        (100.0 * stats_seq.total_pe_active_cycles / stats_seq.total_pe_possible_cycles) : 0.0;
    double util_pipe = (stats_pipe.total_pe_possible_cycles > 0) ?
        (100.0 * stats_pipe.total_pe_active_cycles / stats_pipe.total_pe_possible_cycles) : 0.0;
    
    std::cout << "PE Utilization:\n";
    std::cout << "    Sequential:  " << std::setw(12) << std::setprecision(2) << util_seq << "%\n";
    std::cout << "    Pipelined:   " << std::setw(12) << util_pipe << "%\n\n";
    
    // Operations (should be identical)
    std::cout << "Operations (should be identical):\n";
    std::cout << "    MACs:        " << std::setw(12) << stats_seq.mac_operations 
              << " (seq)  " << std::setw(12) << stats_pipe.mac_operations << " (pipe)\n";
    std::cout << "    DIVs:        " << std::setw(12) << stats_seq.div_operations 
              << " (seq)  " << std::setw(12) << stats_pipe.div_operations << " (pipe)\n";
    
    // Summary
    std::cout << "\n------------------------------------------------------------------------\n";
    std::cout << "CONCLUSION:\n";
    if (speedup > 1.05) {
        std::cout << "    Pipelining provides " << std::setprecision(2) << speedup 
                  << "x speedup by overlapping memory and computation.\n";
        std::cout << "    Case 4 (trailing update) cycles reduced by " 
                  << std::setprecision(1) << case4_reduction << "%.\n";
    } else {
        std::cout << "    Pipelining provides minimal benefit for this configuration.\n";
        std::cout << "    Consider larger matrices with more trailing blocks per row.\n";
    }
    std::cout << "========================================================================\n";
}

void printPartitionedComparison(const SimStats& stats_2d,
                                 const SimStats& stats_2d_pipe,
                                 const SimStats3D& stats_2d_part,
                                 const SimStats3D& stats_3d_old,
                                 const SimStats3D& stats_3d_new,
                                 uint32_t matrix_size,
                                 uint32_t pe_2d, uint32_t pe_part,
                                 uint32_t block_2d, uint32_t block_part,
                                 uint32_t tsv_latency) {
    std::cout << "\n";
    std::cout << "=============================================================================================================================\n";
    std::cout << "          2D vs 2D-Pipelined vs 2D-Partitioned vs 3D-Old vs 3D-New Comparison (Same Total PE Budget)                        \n";
    std::cout << "=============================================================================================================================\n";
    std::cout << "\n";
    
    uint32_t total_pes = pe_2d * pe_2d;
    
    std::cout << "Configuration:\n";
    std::cout << "    Matrix size:       " << matrix_size << " x " << matrix_size << "\n";
    std::cout << "    Total PE budget:   " << total_pes << " PEs\n\n";
    
    std::cout << "Architecture Breakdown:\n";
    std::cout << "    ┌─────────────────┬──────────────┬──────────────┬──────────────┬──────────────┬──────────────┐\n";
    std::cout << "    │                 │ 2D Standard  │ 2D Pipelined │ 2D Partition │ 3D Old (V1)  │ 3D New (V2)  │\n";
    std::cout << "    ├─────────────────┼──────────────┼──────────────┼──────────────┼──────────────┼──────────────┤\n";
    std::cout << "    │ PE Array        │ " << std::setw(4) << pe_2d << " x " << std::setw(4) << pe_2d 
              << "   │ " << std::setw(4) << pe_2d << " x " << std::setw(4) << pe_2d
              << "   │ 4x(" << std::setw(2) << pe_part << "x" << std::setw(2) << pe_part << ")    │ "
              << std::setw(2) << pe_part << "x" << std::setw(2) << pe_part << "x4      │ "
              << std::setw(2) << pe_part << "x" << std::setw(2) << pe_part << "x4      │\n";
    std::cout << "    │ Block Size      │ " << std::setw(4) << block_2d << " x " << std::setw(4) << block_2d 
              << "   │ " << std::setw(4) << block_2d << " x " << std::setw(4) << block_2d
              << "   │ " << std::setw(4) << block_part << " x " << std::setw(4) << block_part 
              << "   │ " << std::setw(4) << block_part << " x " << std::setw(4) << block_part 
              << "   │ " << std::setw(4) << block_part << " x " << std::setw(4) << block_part << "   │\n";
    std::cout << "    │ Case 4 Mode     │ Sequential   │ Skewed Pipe  │ Row-Parallel │ Arb. Assign  │ Row-Parallel │\n";
    std::cout << "    │ Inter-region    │ N/A          │ N/A          │ 2D wire(0cy) │ TSV(" << tsv_latency << "cy/hop)│ TSV(" << tsv_latency << "cy/hop)│\n";
    std::cout << "    └─────────────────┴──────────────┴──────────────┴──────────────┴──────────────┴──────────────┘\n\n";
    
    // Total cycles comparison
    std::cout << "Performance Results:\n";
    std::cout << "    ┌─────────────────┬──────────────┬──────────────┬──────────────┬──────────────┬──────────────┐\n";
    std::cout << "    │ Metric          │ 2D Standard  │ 2D Pipelined │ 2D Partition │ 3D Old (V1)  │ 3D New (V2)  │\n";
    std::cout << "    ├─────────────────┼──────────────┼──────────────┼──────────────┼──────────────┼──────────────┤\n";
    
    std::cout << "    │ Total Cycles    │ " << std::setw(12) << stats_2d.total_cycles 
              << " │ " << std::setw(12) << stats_2d_pipe.total_cycles 
              << " │ " << std::setw(12) << stats_2d_part.total_cycles 
              << " │ " << std::setw(12) << stats_3d_old.total_cycles 
              << " │ " << std::setw(12) << stats_3d_new.total_cycles << " │\n";
    
    // Speedups relative to 2D standard
    double speedup_pipe = static_cast<double>(stats_2d.total_cycles) / 
                          std::max(stats_2d_pipe.total_cycles, (uint64_t)1);
    double speedup_part = static_cast<double>(stats_2d.total_cycles) / 
                          std::max(stats_2d_part.total_cycles, (uint64_t)1);
    double speedup_3d_old = static_cast<double>(stats_2d.total_cycles) / 
                            std::max(stats_3d_old.total_cycles, (uint64_t)1);
    double speedup_3d_new = static_cast<double>(stats_2d.total_cycles) / 
                            std::max(stats_3d_new.total_cycles, (uint64_t)1);
    
    std::cout << "    │ vs 2D Speedup   │ " << std::setw(12) << "1.00x" 
              << " │ " << std::setw(11) << std::fixed << std::setprecision(2) << speedup_pipe << "x"
              << " │ " << std::setw(11) << speedup_part << "x"
              << " │ " << std::setw(11) << speedup_3d_old << "x"
              << " │ " << std::setw(11) << speedup_3d_new << "x │\n";
    
    std::cout << "    ├─────────────────┼──────────────┼──────────────┼──────────────┼──────────────┼──────────────┤\n";
    
    // Case breakdown
    std::cout << "    │ Case 1 (diag)   │ " << std::setw(12) << stats_2d.case1_cycles 
              << " │ " << std::setw(12) << stats_2d_pipe.case1_cycles 
              << " │ " << std::setw(12) << stats_2d_part.case1_cycles 
              << " │ " << std::setw(12) << stats_3d_old.case1_cycles 
              << " │ " << std::setw(12) << stats_3d_new.case1_cycles << " │\n";
    std::cout << "    │ Case 2 (horiz)  │ " << std::setw(12) << stats_2d.case2_cycles 
              << " │ " << std::setw(12) << stats_2d_pipe.case2_cycles 
              << " │ " << std::setw(12) << stats_2d_part.case2_cycles 
              << " │ " << std::setw(12) << stats_3d_old.case2_cycles 
              << " │ " << std::setw(12) << stats_3d_new.case2_cycles << " │\n";
    std::cout << "    │ Case 3 (vert)   │ " << std::setw(12) << stats_2d.case3_cycles 
              << " │ " << std::setw(12) << stats_2d_pipe.case3_cycles 
              << " │ " << std::setw(12) << stats_2d_part.case3_cycles 
              << " │ " << std::setw(12) << stats_3d_old.case3_cycles 
              << " │ " << std::setw(12) << stats_3d_new.case3_cycles << " │\n";
    std::cout << "    │ Case 4 (trail)  │ " << std::setw(12) << stats_2d.case4_cycles 
              << " │ " << std::setw(12) << stats_2d_pipe.case4_cycles 
              << " │ " << std::setw(12) << stats_2d_part.case4_cycles 
              << " │ " << std::setw(12) << stats_3d_old.case4_cycles 
              << " │ " << std::setw(12) << stats_3d_new.case4_cycles << " │\n";
    
    std::cout << "    ├─────────────────┼──────────────┼──────────────┼──────────────┼──────────────┼──────────────┤\n";
    
    // Case 4 reduction percentage
    double case4_red_pipe = (stats_2d.case4_cycles > 0) ?
        (100.0 * (1.0 - static_cast<double>(stats_2d_pipe.case4_cycles) / stats_2d.case4_cycles)) : 0.0;
    double case4_red_part = (stats_2d.case4_cycles > 0) ?
        (100.0 * (1.0 - static_cast<double>(stats_2d_part.case4_cycles) / stats_2d.case4_cycles)) : 0.0;
    double case4_red_3d_old = (stats_2d.case4_cycles > 0) ?
        (100.0 * (1.0 - static_cast<double>(stats_3d_old.case4_cycles) / stats_2d.case4_cycles)) : 0.0;
    double case4_red_3d_new = (stats_2d.case4_cycles > 0) ?
        (100.0 * (1.0 - static_cast<double>(stats_3d_new.case4_cycles) / stats_2d.case4_cycles)) : 0.0;
    
    std::cout << "    │ Case4 Reduction │ " << std::setw(12) << "baseline"
              << " │ " << std::setw(10) << std::setprecision(1) << case4_red_pipe << "%"
              << " │ " << std::setw(10) << case4_red_part << "%"
              << " │ " << std::setw(10) << case4_red_3d_old << "%"
              << " │ " << std::setw(10) << case4_red_3d_new << "% │\n";
    
    std::cout << "    ├─────────────────┼──────────────┼──────────────┼──────────────┼──────────────┼──────────────┤\n";
    
    // Utilization
    double util_2d = (stats_2d.total_pe_possible_cycles > 0) ?
        (100.0 * stats_2d.total_pe_active_cycles / stats_2d.total_pe_possible_cycles) : 0.0;
    double util_pipe = (stats_2d_pipe.total_pe_possible_cycles > 0) ?
        (100.0 * stats_2d_pipe.total_pe_active_cycles / stats_2d_pipe.total_pe_possible_cycles) : 0.0;
    double util_part = (stats_2d_part.total_pe_possible_cycles > 0) ?
        (100.0 * stats_2d_part.total_pe_active_cycles / stats_2d_part.total_pe_possible_cycles) : 0.0;
    double util_3d_old = (stats_3d_old.total_pe_possible_cycles > 0) ?
        (100.0 * stats_3d_old.total_pe_active_cycles / stats_3d_old.total_pe_possible_cycles) : 0.0;
    double util_3d_new = (stats_3d_new.total_pe_possible_cycles > 0) ?
        (100.0 * stats_3d_new.total_pe_active_cycles / stats_3d_new.total_pe_possible_cycles) : 0.0;
    
    std::cout << "    │ PE Utilization  │ " << std::setw(10) << std::setprecision(1) << util_2d << "% "
              << "│ " << std::setw(10) << util_pipe << "% "
              << "│ " << std::setw(10) << util_part << "% "
              << "│ " << std::setw(10) << util_3d_old << "% "
              << "│ " << std::setw(10) << util_3d_new << "% │\n";
    
    std::cout << "    └─────────────────┴──────────────┴──────────────┴──────────────┴──────────────┴──────────────┘\n\n";
    
    // TSV overhead analysis
    std::cout << "TSV Overhead Analysis:\n";
    std::cout << "    2D/2D-Pipelined/2D-Partitioned TSV cycles: 0 (no vertical stacking)\n";
    std::cout << "    3D Old TSV transfer cycles: " << stats_3d_old.tsv_transfer_cycles << "\n";
    std::cout << "    3D New TSV transfer cycles: " << stats_3d_new.tsv_transfer_cycles << "\n";
    
    double tsv_overhead_old = (stats_3d_old.total_cycles > 0) ?
        (100.0 * stats_3d_old.tsv_transfer_cycles / stats_3d_old.total_cycles) : 0.0;
    double tsv_overhead_new = (stats_3d_new.total_cycles > 0) ?
        (100.0 * stats_3d_new.tsv_transfer_cycles / stats_3d_new.total_cycles) : 0.0;
    std::cout << "    TSV overhead (Old/New):     " << std::setprecision(1) << tsv_overhead_old << "% / " << tsv_overhead_new << "%\n\n";
    
    // V1 vs V2 comparison
    std::cout << "3D Case 4 V1 vs V2 Comparison:\n";
    double v2_vs_v1_speedup = static_cast<double>(stats_3d_old.case4_cycles) / 
                              std::max(stats_3d_new.case4_cycles, (uint64_t)1);
    double case4_v2_improvement = (stats_3d_old.case4_cycles > 0) ?
        (100.0 * (1.0 - static_cast<double>(stats_3d_new.case4_cycles) / stats_3d_old.case4_cycles)) : 0.0;
    std::cout << "    V1 (arbitrary block assign):  " << stats_3d_old.case4_cycles << " cycles\n";
    std::cout << "    V2 (row-parallel streaming):  " << stats_3d_new.case4_cycles << " cycles\n";
    std::cout << "    V2 speedup over V1:           " << std::setprecision(2) << v2_vs_v1_speedup << "x";
    if (case4_v2_improvement > 0) {
        std::cout << " (" << std::setprecision(1) << case4_v2_improvement << "% reduction)\n\n";
    } else {
        std::cout << " (" << std::setprecision(1) << -case4_v2_improvement << "% increase)\n\n";
    }
    
    // Find the best performer
    uint64_t min_cycles = std::min({stats_2d.total_cycles, stats_2d_pipe.total_cycles, 
                                    stats_2d_part.total_cycles, stats_3d_old.total_cycles,
                                    stats_3d_new.total_cycles});
    std::string best_arch;
    if (min_cycles == stats_2d.total_cycles) best_arch = "2D Standard";
    else if (min_cycles == stats_2d_pipe.total_cycles) best_arch = "2D Pipelined";
    else if (min_cycles == stats_2d_part.total_cycles) best_arch = "2D Partitioned";
    else if (min_cycles == stats_3d_old.total_cycles) best_arch = "3D Old (V1)";
    else best_arch = "3D New (V2)";
    
    // Summary
    /*
    std::cout << "-----------------------------------------------------------------------------------------------------------------------------\n";
    std::cout << "ANALYSIS:\n";
    
    std::cout << "    • BEST PERFORMER: " << best_arch << " with " << min_cycles << " cycles\n\n";
    
    if (speedup_pipe > 1.05) {
        std::cout << "    • 2D Pipelined: " << std::setprecision(2) << speedup_pipe 
                  << "x speedup via skewed systolic Case 4 (reduces Case 4 by " 
                  << std::setprecision(1) << case4_red_pipe << "%).\n";
    } else {
        std::cout << "    • 2D Pipelined: minimal speedup (" << std::setprecision(2) << speedup_pipe 
                  << "x) - pipelining overhead dominates for small matrices.\n";
    }
    
    if (speedup_part > 1.05) {
        std::cout << "    • 2D Partitioned: " << std::setprecision(2) << speedup_part 
                  << "x speedup via 4-way spatial parallelism.\n";
    } else if (speedup_part < 0.95) {
        std::cout << "    • 2D Partitioned: " << std::setprecision(2) << (1.0/speedup_part)
                  << "x SLOWER - smaller blocks increase block count overhead.\n";
    } else {
        std::cout << "    • 2D Partitioned: similar to 2D Standard (" << std::setprecision(2) << speedup_part << "x).\n";
    }
    
    std::cout << "    • 3D V1→V2 improvement: " << std::setprecision(2) << v2_vs_v1_speedup 
              << "x Case 4 speedup from row-parallel U streaming.\n";
    
    std::cout << "\n    Key Insight: V2's row-parallel design assigns different L blocks to layers\n";
    std::cout << "    and streams the SAME U sequence to all layers simultaneously, maximizing\n";
    std::cout << "    both spatial (across layers) and temporal (U streaming) parallelism.\n";
    std::cout << "=============================================================================================================================\n";
    */
}

int main(int argc, char* argv[]) {
    SimConfig config;
    uint32_t seed = 42;
    
    // 3D-specific options
    bool use_3d = false;
    bool compare_mode = false;
    bool equal_pe_mode = false;
    bool compare_pipeline_mode = false;
    bool compare_part_mode = false;
    uint32_t num_layers = 4;
    uint32_t tsv_latency = 1;
    uint32_t contention_mode = 1;       // 0=ideal, 1=realistic (default), 2=pessimistic
    bool     contention_mode_set = false;
    
    // Parse command line arguments
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-n") == 0 && i + 1 < argc) {
            config.matrix_size = std::atoi(argv[++i]);
        }
        else if (strcmp(argv[i], "-p") == 0 && i + 1 < argc) {
            config.pe_array_size = std::atoi(argv[++i]);
            // Default: block_size = pe_array_size unless explicitly set
            config.block_size = config.pe_array_size;
        }
        else if (strcmp(argv[i], "-b") == 0 && i + 1 < argc) {
            config.block_size = std::atoi(argv[++i]);
        }
        else if (strcmp(argv[i], "-mac") == 0 && i + 1 < argc) {
            config.mac_latency = std::atoi(argv[++i]);
        }
        else if (strcmp(argv[i], "-div") == 0 && i + 1 < argc) {
            config.div_latency = std::atoi(argv[++i]);
        }
        else if (strcmp(argv[i], "-load") == 0 && i + 1 < argc) {
            config.mem_load_delay = std::atoi(argv[++i]);
        }
        else if (strcmp(argv[i], "-write") == 0 && i + 1 < argc) {
            config.mem_write_delay = std::atoi(argv[++i]);
        }
        else if (strcmp(argv[i], "-seed") == 0 && i + 1 < argc) {
            seed = std::atoi(argv[++i]);
        }
        else if (strcmp(argv[i], "-v") == 0) {
            config.verbose = true;
        }
        else if (strcmp(argv[i], "-pipeline") == 0) {
            config.pipeline_enabled = true;
        }
        else if (strcmp(argv[i], "-3d") == 0) {
            use_3d = true;
        }
        else if (strcmp(argv[i], "-z") == 0 && i + 1 < argc) {
            num_layers = std::atoi(argv[++i]);
            use_3d = true;  // Implied
        }
        else if (strcmp(argv[i], "-tsv") == 0 && i + 1 < argc) {
            tsv_latency = std::atoi(argv[++i]);
            use_3d = true;  // Implied
        }
        else if (strcmp(argv[i], "-contention") == 0 && i + 1 < argc) {
            int mode = std::atoi(argv[++i]);
            if (mode < 0 || mode > 2) {
                std::cerr << "Error: -contention must be 0 (ideal), 1 (realistic), or 2 (pessimistic).\n";
                return 1;
            }
            contention_mode = static_cast<uint32_t>(mode);
            contention_mode_set = true;
        }
        else if (strcmp(argv[i], "-compare") == 0) {
            compare_mode = true;
        }
        else if (strcmp(argv[i], "-compare-pipeline") == 0) {
            compare_pipeline_mode = true;
        }
        else if (strcmp(argv[i], "-compare-part") == 0) {
            compare_part_mode = true;
        }
        else if (strcmp(argv[i], "-equal-pe") == 0) {
            equal_pe_mode = true;
        }
        else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            printUsage(argv[0]);
            return 0;
        }
        else {
            std::cerr << "Unknown option: " << argv[i] << "\n";
            printUsage(argv[0]);
            return 1;
        }
    }
    
    // Validate configuration
    if (!config.validate()) {
        std::cerr << "Error: Invalid configuration.\n";
        std::cerr << "Requirements:\n";
        std::cerr << "  - Matrix size must be a power of 2\n";
        std::cerr << "  - PE array size must be a power of 2\n";
        std::cerr << "  - Block size must be a power of 2\n";
        std::cerr << "  - Matrix size >= PE array size\n";
        std::cerr << "  - Block size <= PE array size\n";
        std::cerr << "  - Matrix size must be divisible by block size\n";
        return 1;
    }
    
    // Validate 3D-specific parameters
    if ((use_3d || compare_mode) && num_layers == 0) {
        std::cerr << "Error: Number of layers must be > 0\n";
        return 1;
    }
    
    // Pipeline comparison mode (2D only: sequential vs pipelined)
    if (compare_pipeline_mode) {
        std::cout << "========================================\n";
        std::cout << "Pipeline Comparison: Sequential vs Pipelined\n";
        std::cout << "========================================\n\n";
        
        // Run sequential (non-pipelined) simulation
        std::cout << "=== Running Sequential (Non-Pipelined) Simulation ===\n";
        SimConfig config_seq = config;
        config_seq.pipeline_enabled = false;
        BlockLUSimulator sim_seq(config_seq);
        sim_seq.initializeRandom(seed);
        sim_seq.run();
        bool verified_seq = sim_seq.verify();
        std::cout << "Verification: " << (verified_seq ? "PASSED" : "FAILED") << "\n\n";
        
        // Run pipelined simulation
        std::cout << "=== Running Pipelined Simulation ===\n";
        SimConfig config_pipe = config;
        config_pipe.pipeline_enabled = true;
        BlockLUSimulator sim_pipe(config_pipe);
        sim_pipe.initializeRandom(seed);
        sim_pipe.run();
        bool verified_pipe = sim_pipe.verify();
        std::cout << "Verification: " << (verified_pipe ? "PASSED" : "FAILED") << "\n";
        
        // Print comparison
        printPipelineComparison(sim_seq.getStats(), sim_pipe.getStats(),
                                config.matrix_size, config.pe_array_size, config.block_size);
        
        return (verified_seq && verified_pipe) ? 0 : 1;
    }
    
    // 2D vs 2D-Partitioned vs 3D comparison mode
    if (compare_part_mode) {
        // Validate: PE array must be divisible by 2 for 4-way partitioning
        if (config.pe_array_size < 4 || config.pe_array_size % 2 != 0) {
            std::cerr << "Error: For -compare-part, PE array size must be >= 4 and even.\n";
            std::cerr << "  Current: p = " << config.pe_array_size << "\n";
            return 1;
        }
        
        uint32_t pe_2d = config.pe_array_size;
        uint32_t block_2d = config.block_size;
        uint32_t pe_part = config.pe_array_size / 2;  // Each partition is (p/2) x (p/2)
        uint32_t block_part = pe_part;  // Block size matches partition size
        
        // Validate matrix size is compatible with partitioned block size
        if (config.matrix_size % block_part != 0) {
            std::cerr << "Error: Matrix size must be divisible by partitioned block size.\n";
            std::cerr << "  Matrix: " << config.matrix_size << ", Part block size: " << block_part << "\n";
            return 1;
        }
        
        std::cout << "========================================================================\n";
        std::cout << "2D vs 2D-Pipelined vs 2D-Partitioned vs 3D (V1/V2) Comparison\n";
        std::cout << "========================================================================\n\n";
        
        std::cout << "Configuration:\n";
        std::cout << "  Matrix size:          " << config.matrix_size << " x " << config.matrix_size << "\n";
        std::cout << "  Total PE budget:      " << (pe_2d * pe_2d) << " PEs\n";
        std::cout << "  2D Standard:          " << pe_2d << " x " << pe_2d << " PEs, block " << block_2d << " (sequential Case 4)\n";
        std::cout << "  2D Pipelined:         " << pe_2d << " x " << pe_2d << " PEs, block " << block_2d << " (skewed systolic)\n";
        std::cout << "  2D Partitioned:       4 x (" << pe_part << " x " << pe_part << ") PEs, block " << block_part << " (row-parallel)\n";
        std::cout << "  3D Old (V1):          " << pe_part << " x " << pe_part << " x 4 PEs, block " << block_part << " (arbitrary assign)\n";
        std::cout << "  3D New (V2):          " << pe_part << " x " << pe_part << " x 4 PEs, block " << block_part << " (row-parallel)\n";
        std::cout << "  TSV latency:          " << tsv_latency << " cycle(s)/hop\n\n";
        
        // 1. Run 2D Standard simulation (non-pipelined)
        std::cout << "=== Running 2D Standard Simulation ===\n";
        SimConfig config_2d = config;
        config_2d.pipeline_enabled = false;
        BlockLUSimulator sim_2d(config_2d);
        sim_2d.initializeRandom(seed);
        sim_2d.run();
        bool verified_2d = sim_2d.verify();
        std::cout << "Verification: " << (verified_2d ? "PASSED" : "FAILED") << "\n\n";
        
        // 2. Run 2D Pipelined simulation
        std::cout << "=== Running 2D Pipelined Simulation ===\n";
        SimConfig config_2d_pipe = config;
        config_2d_pipe.pipeline_enabled = true;
        BlockLUSimulator sim_2d_pipe(config_2d_pipe);
        sim_2d_pipe.initializeRandom(seed);
        sim_2d_pipe.run();
        bool verified_2d_pipe = sim_2d_pipe.verify();
        std::cout << "Verification: " << (verified_2d_pipe ? "PASSED" : "FAILED") << "\n\n";
        
        // 3. Run 2D Partitioned simulation (3D simulator with TSV=0, pipelined v2)
        std::cout << "=== Running 2D Partitioned Simulation (4 regions, pipelined) ===\n";
        SimConfig config_part_base;
        config_part_base.matrix_size = config.matrix_size;
        config_part_base.pe_array_size = pe_part;
        config_part_base.block_size = block_part;
        config_part_base.mac_latency = config.mac_latency;
        config_part_base.div_latency = config.div_latency;
        config_part_base.mem_load_delay = config.mem_load_delay;
        config_part_base.mem_write_delay = config.mem_write_delay;
        config_part_base.verbose = config.verbose;
        config_part_base.pipeline_enabled = true;  // Enable pipelining within each partition
        
        SimConfig3D config_part(config_part_base);
        config_part.num_layers = 4;
        config_part.tsv_latency = 0;  // No TSV overhead for 2D partitioned
        config_part.pipeline_version = 2;  // Use row-parallel (v2)
        // Bandwidth contention model: 2D Partitioned has a shared memory bus
        config_part.is_2d_partitioned = true;
        config_part.enable_bandwidth_model = contention_mode_set;
        config_part.contention_mode = contention_mode;
        config_part.memory_bandwidth_per_layer = config_part_base.block_size;
        BlockLUSimulator3D sim_part(config_part);
        sim_part.run();
        bool verified_part = sim_part.verify();
        std::cout << "Verification: " << (verified_part ? "PASSED" : "FAILED") << "\n\n";
        
        // 4. Run 3D simulation with V1 (old - arbitrary block assignment)
        std::cout << "=== Running 3D Old (V1) Simulation (Z=4, TSV=" << tsv_latency << ", arb. assign) ===\n";
        SimConfig3D config_3d_old(config_part_base);
        config_3d_old.num_layers = 4;
        config_3d_old.tsv_latency = tsv_latency;
        config_3d_old.pipeline_version = 1;  // V1: arbitrary block assignment
        // 3D: independent per-layer ports, contention model gives a flat 3x slowdown
        config_3d_old.is_2d_partitioned = false;
        config_3d_old.enable_bandwidth_model = contention_mode_set;
        config_3d_old.contention_mode = contention_mode;
        config_3d_old.memory_bandwidth_per_layer = config_part_base.block_size;
        BlockLUSimulator3D sim_3d_old(config_3d_old);
        sim_3d_old.run();
        bool verified_3d_old = sim_3d_old.verify();
        std::cout << "Verification: " << (verified_3d_old ? "PASSED" : "FAILED") << "\n\n";
        
        // 5. Run 3D simulation with V2 (new - row-parallel streaming)
        std::cout << "=== Running 3D New (V2) Simulation (Z=4, TSV=" << tsv_latency << ", row-parallel) ===\n";
        SimConfig3D config_3d_new(config_part_base);
        config_3d_new.num_layers = 4;
        config_3d_new.tsv_latency = tsv_latency;
        config_3d_new.pipeline_version = 2;  // V2: row-parallel streaming
        // 3D: independent per-layer ports, contention model gives a flat 3x slowdown
        config_3d_new.is_2d_partitioned = false;
        config_3d_new.enable_bandwidth_model = contention_mode_set;
        config_3d_new.contention_mode = contention_mode;
        config_3d_new.memory_bandwidth_per_layer = config_part_base.block_size;
        BlockLUSimulator3D sim_3d_new(config_3d_new);
        sim_3d_new.run();
        bool verified_3d_new = sim_3d_new.verify();
        std::cout << "Verification: " << (verified_3d_new ? "PASSED" : "FAILED") << "\n";
        
        // Print comparison
        printPartitionedComparison(sim_2d.getStats(), sim_2d_pipe.getStats(),
                                   sim_part.getStats(), sim_3d_old.getStats(), sim_3d_new.getStats(),
                                   config.matrix_size, pe_2d, pe_part, block_2d, block_part,
                                   tsv_latency);
        
        return (verified_2d && verified_2d_pipe && verified_part && verified_3d_old && verified_3d_new) ? 0 : 1;
    }
    
    // Equal PE comparison mode
    if (equal_pe_mode) {
        // Validate constraints for equal PE mode
        if (!isPerfectSquare(num_layers)) {
            std::cerr << "Error: For -equal-pe mode, Z (number of layers) must be a perfect square.\n";
            std::cerr << "  Valid values: 1, 4, 9, 16, 25, ...\n";
            std::cerr << "  Current value: Z = " << num_layers << "\n";
            return 1;
        }
        
        uint32_t sqrt_z = intSqrt(num_layers);
        if (config.pe_array_size % sqrt_z != 0) {
            std::cerr << "Error: For -equal-pe mode, sqrt(Z) must divide PE array size (p).\n";
            std::cerr << "  p = " << config.pe_array_size << ", Z = " << num_layers 
                      << ", sqrt(Z) = " << sqrt_z << "\n";
            std::cerr << "  " << config.pe_array_size << " % " << sqrt_z << " = " 
                      << (config.pe_array_size % sqrt_z) << " (must be 0)\n";
            return 1;
        }
        
        // Calculate 3D configuration
        uint32_t pe_2d = config.pe_array_size;  // 2D PE array: p x p
        uint32_t block_2d = config.block_size;  // 2D block size: b
        uint32_t pe_3d = config.pe_array_size / sqrt_z;  // 3D PE array: (p/√Z) x (p/√Z)
        uint32_t block_3d = pe_3d;  // 3D block size: same as 3D PE array size
        uint32_t total_pes = pe_2d * pe_2d;  // Total PE budget
        
        // Validate 3D block size is compatible with matrix
        if (config.matrix_size % block_3d != 0) {
            std::cerr << "Error: Matrix size must be divisible by 3D block size.\n";
            std::cerr << "  Matrix size = " << config.matrix_size << ", 3D block size = " << block_3d << "\n";
            return 1;
        }
        
        std::cout << "========================================\n";
        std::cout << "Equal PE Budget: 2D vs 3D Comparison\n";
        std::cout << "========================================\n\n";
        
        std::cout << "Configuration:\n";
        std::cout << "  Matrix size:       " << config.matrix_size << " x " << config.matrix_size << "\n";
        std::cout << "  Total PE budget:   " << total_pes << " PEs\n";
        std::cout << "  2D PE array:       " << pe_2d << " x " << pe_2d << " = " << (pe_2d * pe_2d) << " PEs\n";
        std::cout << "  2D block size:     " << block_2d << " x " << block_2d << "\n";
        std::cout << "  3D PE array:       " << pe_3d << " x " << pe_3d << " x " << num_layers 
                  << " = " << (pe_3d * pe_3d * num_layers) << " PEs\n";
        std::cout << "  3D block size:     " << block_3d << " x " << block_3d << "\n";
        std::cout << "  TSV latency:       " << tsv_latency << " cycle(s)/hop\n\n";
        
        // Run 2D simulation
        std::cout << "=== Running 2D Simulation ===\n";
        SimConfig config_2d = config;
        config_2d.pe_array_size = pe_2d;
        config_2d.block_size = block_2d;
        BlockLUSimulator sim_2d(config_2d);
        sim_2d.initializeRandom(seed);
        sim_2d.run();
        bool verified_2d = sim_2d.verify();
        std::cout << "Verification: " << (verified_2d ? "PASSED" : "FAILED") << "\n\n";
        
        // Run 3D simulation with reduced PE array
        std::cout << "=== Running 3D Simulation (" << pe_3d << "x" << pe_3d << "x" << num_layers << ") ===\n";
        SimConfig config_3d_base;
        config_3d_base.matrix_size = config.matrix_size;
        config_3d_base.pe_array_size = pe_3d;
        config_3d_base.block_size = block_3d;
        config_3d_base.mac_latency = config.mac_latency;
        config_3d_base.div_latency = config.div_latency;
        config_3d_base.mem_load_delay = config.mem_load_delay;
        config_3d_base.mem_write_delay = config.mem_write_delay;
        config_3d_base.verbose = config.verbose;
        
        SimConfig3D config_3d(config_3d_base);
        config_3d.num_layers = num_layers;
        config_3d.tsv_latency = tsv_latency;
        BlockLUSimulator3D sim_3d(config_3d);
        sim_3d.run();
        bool verified_3d = sim_3d.verify();
        std::cout << "Verification: " << (verified_3d ? "PASSED" : "FAILED") << "\n";
        
        // Print comparison
        printEqualPEComparison(sim_2d.getStats(), sim_3d.getStats(),
                               pe_2d, pe_3d, block_2d, block_3d, num_layers, total_pes);
        
        return (verified_2d && verified_3d) ? 0 : 1;
    }
    
    if (compare_mode) {
        // Run both 2D and 3D simulations for comparison
        std::cout << "========================================\n";
        std::cout << "2D vs 3D Block LU Decomposition Comparison\n";
        std::cout << "========================================\n\n";
        
        std::cout << "Configuration:\n";
        std::cout << "  Matrix size:   " << config.matrix_size << " x " << config.matrix_size << "\n";
        std::cout << "  PE array size: " << config.pe_array_size << " x " << config.pe_array_size << "\n";
        std::cout << "  Block size:    " << config.block_size << " x " << config.block_size << "\n";
        std::cout << "  Z layers:      " << num_layers << "\n";
        std::cout << "  TSV latency:   " << tsv_latency << " cycle(s)/hop\n\n";
        
        // Run 2D simulation
        std::cout << "=== Running 2D Simulation ===\n";
        BlockLUSimulator sim_2d(config);
        sim_2d.initializeRandom(seed);
        sim_2d.run();
        bool verified_2d = sim_2d.verify();
        std::cout << "Verification: " << (verified_2d ? "PASSED" : "FAILED") << "\n\n";
        
        // Run 3D simulation
        std::cout << "=== Running 3D Simulation (" << num_layers << " layers) ===\n";
        SimConfig3D config_3d(config);
        config_3d.num_layers = num_layers;
        config_3d.tsv_latency = tsv_latency;
        BlockLUSimulator3D sim_3d(config_3d);
        sim_3d.run();
        bool verified_3d = sim_3d.verify();
        std::cout << "Verification: " << (verified_3d ? "PASSED" : "FAILED") << "\n";
        
        // Print comparison
        printComparison(sim_2d.getStats(), sim_3d.getStats(), num_layers);
        
        return (verified_2d && verified_3d) ? 0 : 1;
    }
    else if (use_3d) {
        // Run 3D simulation only
        std::cout << "========================================\n";
        std::cout << "3D Block LU Decomposition Simulator\n";
        std::cout << "========================================\n\n";
        
        SimConfig3D config_3d(config);
        config_3d.num_layers = num_layers;
        config_3d.tsv_latency = tsv_latency;
        config_3d.print();
        
        BlockLUSimulator3D simulator(config_3d);
        simulator.run();
        
        // Verify result
        std::cout << "=== Verification ===\n";
        bool verified = simulator.verify();
        if (verified) {
            std::cout << "Result: PASSED\n\n";
        } else {
            std::cout << "Result: FAILED\n\n";
        }
        
        // Print results
        simulator.printResults();
        
        return verified ? 0 : 1;
    }
    else {
        // Run 2D simulation only (original behavior)
        std::cout << "========================================\n";
        std::cout << "2D Block LU Decomposition Simulator\n";
        std::cout << "========================================\n\n";
        
        BlockLUSimulator simulator(config);
        simulator.initializeRandom(seed);
        
        simulator.run();
        
        // Verify result
        std::cout << "=== Verification ===\n";
        bool verified = simulator.verify();
        if (verified) {
            std::cout << "Result: PASSED\n\n";
        } else {
            std::cout << "Result: FAILED\n\n";
        }
        
        // Print results
        simulator.printResults();
        
        return verified ? 0 : 1;
    }
}
