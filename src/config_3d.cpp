#include "config_3d.hpp"
#include <iostream>
#include <iomanip>

void SimConfig3D::print() const {
    std::cout << "=== 3D Block LU Decomposition Simulator Configuration ===\n";
    std::cout << "Matrix size:        " << matrix_size << " x " << matrix_size << "\n";
    std::cout << "PE array size:      " << pe_array_size << " x " << pe_array_size << "\n";
    std::cout << "Block size:         " << block_size << " x " << block_size << "\n";
    std::cout << "Number of layers:   " << num_layers << " (Z dimension)\n";
    std::cout << "TSV latency:        " << tsv_latency << " cycle(s) per hop\n";
    std::cout << "MAC latency:        " << mac_latency << " cycle(s)\n";
    std::cout << "DIV latency:        " << div_latency << " cycle(s)\n";
    std::cout << "Memory load delay:  " << mem_load_delay << " cycle(s)\n";
    std::cout << "Memory write delay: " << mem_write_delay << " cycle(s)\n";
    std::cout << "Load balance:       " << (isBalancedDistribution() ? "Balanced" : "Unbalanced") << "\n";
    if (!isBalancedDistribution()) {
        std::cout << "  (block_size " << block_size << " % num_layers " << num_layers 
                  << " = " << (block_size % num_layers) << " remainder rows)\n";
    }
    std::cout << "Verbose mode:       " << (verbose ? "ON" : "OFF") << "\n";
    std::cout << "=========================================================\n\n";
}

void SimStats3D::print() const {
    std::cout << "\n=== 3D Simulation Statistics ===\n";
    std::cout << "Total cycles:           " << total_cycles << "\n";
    std::cout << "  Case 1 (diagonal):    " << case1_cycles 
              << " (" << std::fixed << std::setprecision(1) 
              << (100.0 * case1_cycles / total_cycles) << "%)\n";
    std::cout << "  Case 2 (horizontal):  " << case2_cycles
              << " (" << (100.0 * case2_cycles / total_cycles) << "%)\n";
    std::cout << "  Case 3 (vertical):    " << case3_cycles
              << " (" << (100.0 * case3_cycles / total_cycles) << "%)\n";
    std::cout << "  Case 4 (trailing):    " << case4_cycles
              << " (" << (100.0 * case4_cycles / total_cycles) << "%)\n";
    
    std::cout << "\nTSV Statistics:\n";
    std::cout << "  TSV transfer cycles:  " << tsv_transfer_cycles << "\n";
    std::cout << "  TSV transfers:        " << tsv_transfers << "\n";
    
    std::cout << "\nPE Utilization:\n";
    double utilization = (total_pe_possible_cycles > 0) ?
        (100.0 * total_pe_active_cycles / total_pe_possible_cycles) : 0.0;
    std::cout << "  Overall:              " << std::setprecision(2) << utilization << "%\n";
    std::cout << "  Active PE-cycles:     " << total_pe_active_cycles << "\n";
    std::cout << "  Total PE-cycles:      " << total_pe_possible_cycles << "\n";
    
    std::cout << "\nPer-Layer Utilization:\n";
    for (size_t z = 0; z < layer_active_cycles.size(); z++) {
        uint64_t total = layer_active_cycles[z] + layer_idle_cycles[z];
        double layer_util = (total > 0) ? (100.0 * layer_active_cycles[z] / total) : 0.0;
        std::cout << "  Layer " << z << ":              " << std::setprecision(2) 
                  << layer_util << "%\n";
    }
    
    std::cout << "\nOperations:\n";
    std::cout << "  MAC operations:       " << mac_operations << "\n";
    std::cout << "  DIV operations:       " << div_operations << "\n";
    std::cout << "================================\n";
}
