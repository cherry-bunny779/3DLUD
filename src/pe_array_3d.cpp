#include "pe_array_3d.hpp"
#include <iostream>
#include <iomanip>
#include <algorithm>

PEArray3D::PEArray3D(const SimConfig3D& cfg)
    : config(cfg)
    , current_cycle(0)
{
    // Initialize 3D PE grid: [layer][row][col]
    pes.resize(config.num_layers);
    for (uint32_t z = 0; z < config.num_layers; z++) {
        pes[z].resize(config.pe_array_size);
        for (uint32_t r = 0; r < config.pe_array_size; r++) {
            pes[z][r].resize(config.pe_array_size);
            for (uint32_t c = 0; c < config.pe_array_size; c++) {
                pes[z][r][c] = ProcessingElement3D(r, c, z);
            }
        }
    }
}

void PEArray3D::reset() {
    current_cycle = 0;
    for (uint32_t z = 0; z < config.num_layers; z++) {
        for (uint32_t r = 0; r < config.pe_array_size; r++) {
            for (uint32_t c = 0; c < config.pe_array_size; c++) {
                pes[z][r][c].reset();
                pes[z][r][c].layer = z;
            }
        }
    }
}

ProcessingElement3D& PEArray3D::getPE(uint32_t layer, uint32_t row, uint32_t col) {
    return pes[layer][row][col];
}

const ProcessingElement3D& PEArray3D::getPE(uint32_t layer, uint32_t row, uint32_t col) const {
    return pes[layer][row][col];
}

void PEArray3D::tick() {
    current_cycle++;
    for (uint32_t z = 0; z < config.num_layers; z++) {
        for (uint32_t r = 0; r < config.pe_array_size; r++) {
            for (uint32_t c = 0; c < config.pe_array_size; c++) {
                pes[z][r][c].tick();
            }
        }
    }
}

bool PEArray3D::allIdle() const {
    for (uint32_t z = 0; z < config.num_layers; z++) {
        if (!layerIdle(z)) return false;
    }
    return true;
}

bool PEArray3D::layerIdle(uint32_t layer) const {
    for (uint32_t r = 0; r < config.pe_array_size; r++) {
        for (uint32_t c = 0; c < config.pe_array_size; c++) {
            if (pes[layer][r][c].isBusy()) {
                return false;
            }
        }
    }
    return true;
}

uint64_t PEArray3D::waitUntilIdle() {
    uint64_t cycles_waited = 0;
    while (!allIdle()) {
        tick();
        cycles_waited++;
    }
    return cycles_waited;
}

uint64_t PEArray3D::waitUntilLayerIdle(uint32_t layer) {
    uint64_t cycles_waited = 0;
    while (!layerIdle(layer)) {
        tick();
        cycles_waited++;
    }
    return cycles_waited;
}

uint32_t PEArray3D::getMemoryLoadDelay(uint32_t layer, uint32_t pe_row, uint32_t pe_col) const {
    // For 3D, memory access might vary by layer
    // For now, return constant delay (can be extended later)
    (void)layer;
    (void)pe_row;
    (void)pe_col;
    return config.mem_load_delay;
}

uint32_t PEArray3D::getMemoryWriteDelay(uint32_t layer, uint32_t pe_row, uint32_t pe_col) const {
    (void)layer;
    (void)pe_row;
    (void)pe_col;
    return config.mem_write_delay;
}

uint64_t PEArray3D::broadcastToAllLayers(uint32_t src_layer, uint32_t row, uint32_t col, float data) {
    // Broadcast from source layer to all other layers
    // Per-hop model: farthest layer determines total latency
    
    uint32_t max_hops = 0;
    
    for (uint32_t z = 0; z < config.num_layers; z++) {
        if (z != src_layer) {
            uint32_t hops = (z > src_layer) ? (z - src_layer) : (src_layer - z);
            max_hops = std::max(max_hops, hops);
            
            uint32_t latency = hops * config.tsv_latency;
            pes[z][row][col].startTSVReceive(latency);
            pes[z][row][col].setTSVReceivedData(data);
        }
    }
    
    // Source PE sends
    pes[src_layer][row][col].startTSVSend(data, max_hops * config.tsv_latency);
    
    return max_hops * config.tsv_latency;
}

uint64_t PEArray3D::sendToLayer(uint32_t src_layer, uint32_t dst_layer,
                                 uint32_t row, uint32_t col, float data) {
    uint32_t latency = config.getTSVLatency(src_layer, dst_layer);
    
    pes[src_layer][row][col].startTSVSend(data, latency);
    pes[dst_layer][row][col].startTSVReceive(latency);
    pes[dst_layer][row][col].setTSVReceivedData(data);
    
    return latency;
}

uint64_t PEArray3D::aggregateToLayer(uint32_t dst_layer, uint32_t row, uint32_t col) {
    // Aggregate (sum) values from all layers to destination layer
    // Uses pipelined approach: each layer sends to next, accumulating
    
    // For simplicity, we'll do sequential aggregation
    // Layer 0 -> Layer 1 -> ... -> dst_layer
    
    float sum = 0.0f;
    uint64_t total_cycles = 0;
    
    // Collect from layers below dst_layer
    for (uint32_t z = 0; z < dst_layer; z++) {
        sum += pes[z][row][col].reg_result;
        total_cycles += config.tsv_latency;
    }
    
    // Add dst_layer's own value
    sum += pes[dst_layer][row][col].reg_result;
    
    // Collect from layers above dst_layer
    for (uint32_t z = dst_layer + 1; z < config.num_layers; z++) {
        sum += pes[z][row][col].reg_result;
        total_cycles += config.tsv_latency;
    }
    
    pes[dst_layer][row][col].reg_result = sum;
    
    return total_cycles;
}

void PEArray3D::collectStats(SimStats3D& stats) const {
    stats.layer_active_cycles.resize(config.num_layers, 0);
    stats.layer_idle_cycles.resize(config.num_layers, 0);
    
    for (uint32_t z = 0; z < config.num_layers; z++) {
        for (uint32_t r = 0; r < config.pe_array_size; r++) {
            for (uint32_t c = 0; c < config.pe_array_size; c++) {
                const auto& pe = pes[z][r][c];
                
                stats.total_pe_active_cycles += pe.active_cycles;
                stats.total_pe_possible_cycles += pe.active_cycles + pe.idle_cycles;
                
                stats.layer_active_cycles[z] += pe.active_cycles + pe.tsv_active_cycles;
                stats.layer_idle_cycles[z] += pe.idle_cycles;
                
                stats.tsv_transfer_cycles += pe.tsv_active_cycles;
            }
        }
    }
}

void PEArray3D::printState() const {
    std::cout << "Cycle " << current_cycle << " - 3D PE Array State:\n";
    for (uint32_t z = 0; z < config.num_layers; z++) {
        printLayerState(z);
    }
    std::cout << "\n";
}

void PEArray3D::printLayerState(uint32_t layer) const {
    std::cout << "  Layer " << layer << ":\n";
    std::cout << "      ";
    for (uint32_t c = 0; c < config.pe_array_size; c++) {
        std::cout << std::setw(6) << "C" << c;
    }
    std::cout << "\n";
    
    for (uint32_t r = 0; r < config.pe_array_size; r++) {
        std::cout << "  R" << r << ": ";
        for (uint32_t c = 0; c < config.pe_array_size; c++) {
            std::cout << std::setw(6) << peStateToString(pes[layer][r][c].state);
        }
        std::cout << "\n";
    }
}

uint32_t PEArray3D::countActivePEs() const {
    uint32_t count = 0;
    for (uint32_t z = 0; z < config.num_layers; z++) {
        count += countActivePEsOnLayer(z);
    }
    return count;
}

uint32_t PEArray3D::countActivePEsOnLayer(uint32_t layer) const {
    uint32_t count = 0;
    for (uint32_t r = 0; r < config.pe_array_size; r++) {
        for (uint32_t c = 0; c < config.pe_array_size; c++) {
            const auto& pe = pes[layer][r][c];
            if (pe.state == PEState::COMPUTING || 
                pe.state == PEState::LOADING ||
                pe.state == PEState::WRITING ||
                pe.tsv_sending || pe.tsv_receiving) {
                count++;
            }
        }
    }
    return count;
}
