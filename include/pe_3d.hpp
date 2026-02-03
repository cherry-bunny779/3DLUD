#ifndef PE_3D_HPP
#define PE_3D_HPP

#include "pe.hpp"
#include <cstdint>

/**
 * 3D Processing Element class
 * Extends base PE with Z-dimension (layer) information and TSV communication
 */
class ProcessingElement3D : public ProcessingElement {
public:
    // Layer position in 3D stack
    uint32_t layer;
    
    // TSV communication state
    bool tsv_sending;           // Currently sending data via TSV
    bool tsv_receiving;         // Currently receiving data via TSV
    uint32_t tsv_cycles_remaining;
    
    // TSV data registers
    float tsv_send_data;
    float tsv_recv_data;
    
    // Statistics
    uint64_t tsv_active_cycles;
    
    // Constructor
    ProcessingElement3D(uint32_t r = 0, uint32_t c = 0, uint32_t z = 0);
    
    // Reset PE state (overrides base)
    void reset();
    
    // Start TSV send operation to another layer
    void startTSVSend(float data, uint32_t latency);
    
    // Start TSV receive operation from another layer
    void startTSVReceive(uint32_t latency);
    
    // Set received TSV data (called when transfer completes)
    void setTSVReceivedData(float data);
    
    // Get TSV send data
    float getTSVSendData() const { return tsv_send_data; }
    
    // Get TSV received data
    float getTSVReceivedData() const { return tsv_recv_data; }
    
    // Check if TSV operation is in progress
    bool isTSVBusy() const { return tsv_sending || tsv_receiving; }
    
    // Advance one cycle (overrides base to handle TSV)
    bool tick();
    
    // Check if PE is busy (including TSV operations)
    bool isBusy() const;
};

#endif // PE_3D_HPP
