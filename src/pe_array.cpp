#include "pe_array.hpp"
#include <iostream>
#include <iomanip>

PEArray::PEArray(const SimConfig& cfg)
    : config(cfg)
    , current_cycle(0)
{
    // Initialize PE grid
    pes.resize(config.pe_array_size);
    for (uint32_t r = 0; r < config.pe_array_size; r++) {
        pes[r].resize(config.pe_array_size);
        for (uint32_t c = 0; c < config.pe_array_size; c++) {
            pes[r][c] = ProcessingElement(r, c);
        }
    }
}

void PEArray::reset() {
    current_cycle = 0;
    for (uint32_t r = 0; r < config.pe_array_size; r++) {
        for (uint32_t c = 0; c < config.pe_array_size; c++) {
            pes[r][c].reset();
        }
    }
}

ProcessingElement& PEArray::getPE(uint32_t row, uint32_t col) {
    return pes[row][col];
}

const ProcessingElement& PEArray::getPE(uint32_t row, uint32_t col) const {
    return pes[row][col];
}

void PEArray::tick() {
    current_cycle++;
    for (uint32_t r = 0; r < config.pe_array_size; r++) {
        for (uint32_t c = 0; c < config.pe_array_size; c++) {
            pes[r][c].tick();
        }
    }
}

bool PEArray::allIdle() const {
    for (uint32_t r = 0; r < config.pe_array_size; r++) {
        for (uint32_t c = 0; c < config.pe_array_size; c++) {
            if (pes[r][c].isBusy()) {
                return false;
            }
        }
    }
    return true;
}

bool PEArray::allComplete() const {
    for (uint32_t r = 0; r < config.pe_array_size; r++) {
        for (uint32_t c = 0; c < config.pe_array_size; c++) {
            if (!pes[r][c].isComplete()) {
                return false;
            }
        }
    }
    return true;
}

uint64_t PEArray::waitUntilIdle() {
    uint64_t cycles_waited = 0;
    while (!allIdle()) {
        tick();
        cycles_waited++;
    }
    return cycles_waited;
}

uint32_t PEArray::getMemoryLoadDelay(uint32_t pe_row, uint32_t pe_col) const {
    // Dummy function: returns constant delay regardless of PE position
    // In a more sophisticated model, this could depend on:
    // - PE position (edge PEs might have faster access)
    // - Memory bank conflicts
    // - Bus contention
    (void)pe_row;  // Unused for now
    (void)pe_col;  // Unused for now
    return config.mem_load_delay;
}

uint32_t PEArray::getMemoryWriteDelay(uint32_t pe_row, uint32_t pe_col) const {
    // Dummy function: returns constant delay regardless of PE position
    (void)pe_row;  // Unused for now
    (void)pe_col;  // Unused for now
    return config.mem_write_delay;
}

void PEArray::collectStats(SimStats& stats) const {
    for (uint32_t r = 0; r < config.pe_array_size; r++) {
        for (uint32_t c = 0; c < config.pe_array_size; c++) {
            stats.total_pe_active_cycles += pes[r][c].active_cycles;
            stats.total_pe_possible_cycles += pes[r][c].active_cycles + pes[r][c].idle_cycles;
        }
    }
}

void PEArray::printState() const {
    std::cout << "Cycle " << current_cycle << " PE Array State:\n";
    std::cout << "    ";
    for (uint32_t c = 0; c < config.pe_array_size; c++) {
        std::cout << std::setw(6) << "C" << c;
    }
    std::cout << "\n";
    
    for (uint32_t r = 0; r < config.pe_array_size; r++) {
        std::cout << "R" << r << ": ";
        for (uint32_t c = 0; c < config.pe_array_size; c++) {
            std::cout << std::setw(6) << peStateToString(pes[r][c].state);
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}

uint32_t PEArray::countActivePEs() const {
    uint32_t count = 0;
    for (uint32_t r = 0; r < config.pe_array_size; r++) {
        for (uint32_t c = 0; c < config.pe_array_size; c++) {
            if (pes[r][c].state == PEState::COMPUTING || 
                pes[r][c].state == PEState::LOADING ||
                pes[r][c].state == PEState::WRITING) {
                count++;
            }
        }
    }
    return count;
}
