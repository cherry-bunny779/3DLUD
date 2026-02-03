#ifndef PE_HPP
#define PE_HPP

#include "config.hpp"
#include <cstdint>

/**
 * Processing Element (PE) class
 * Each PE can hold one matrix element and perform MAC/DIV operations
 */
class ProcessingElement {
public:
    // Position in the PE array
    uint32_t row;
    uint32_t col;
    
    // Local registers
    float reg_a;        // Input/accumulator register
    float reg_b;        // Second operand register
    float reg_result;   // Result register
    
    // State tracking
    PEState state;
    uint32_t cycles_remaining;  // Cycles until current operation completes
    
    // Statistics
    uint64_t active_cycles;
    uint64_t idle_cycles;
    
    // Constructor
    ProcessingElement(uint32_t r = 0, uint32_t c = 0);
    
    // Reset PE state
    void reset();
    
    // Start a MAC operation: result = a - b * c (for LU: A[i,j] -= L[i,k] * U[k,j])
    void startMAC(float a, float b, float c, uint32_t latency);
    
    // Start a division operation: result = a / b (for LU: L[j,k] = A[j,k] / U[k,k])
    void startDIV(float a, float b, uint32_t latency);
    
    // Start memory load
    void startLoad(uint32_t delay);
    
    // Start memory write
    void startWrite(uint32_t delay);
    
    // Set to waiting state
    void setWaiting();
    
    // Set to idle state
    void setIdle();
    
    // Advance one cycle, returns true if operation completed this cycle
    bool tick();
    
    // Check if PE is busy
    bool isBusy() const;
    
    // Check if PE just completed an operation
    bool isComplete() const;
    
    // Get the result value
    float getResult() const { return reg_result; }
};

#endif // PE_HPP
