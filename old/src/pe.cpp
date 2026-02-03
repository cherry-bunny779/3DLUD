#include "pe.hpp"

ProcessingElement::ProcessingElement(uint32_t r, uint32_t c)
    : row(r)
    , col(c)
    , reg_a(0.0f)
    , reg_b(0.0f)
    , reg_result(0.0f)
    , state(PEState::IDLE)
    , cycles_remaining(0)
    , active_cycles(0)
    , idle_cycles(0)
{}

void ProcessingElement::reset() {
    reg_a = 0.0f;
    reg_b = 0.0f;
    reg_result = 0.0f;
    state = PEState::IDLE;
    cycles_remaining = 0;
    active_cycles = 0;
    idle_cycles = 0;
}

void ProcessingElement::startMAC(float a, float b, float c, uint32_t latency) {
    // Compute: result = a - b * c
    reg_a = a;
    reg_b = b;
    reg_result = a - b * c;  // Actual computation (result available after latency)
    state = PEState::COMPUTING;
    cycles_remaining = latency;
}

void ProcessingElement::startDIV(float a, float b, uint32_t latency) {
    // Compute: result = a / b
    reg_a = a;
    reg_b = b;
    if (b != 0.0f) {
        reg_result = a / b;
    } else {
        reg_result = 0.0f;  // Handle division by zero gracefully
    }
    state = PEState::COMPUTING;
    cycles_remaining = latency;
}

void ProcessingElement::startLoad(uint32_t delay) {
    state = PEState::LOADING;
    cycles_remaining = delay;
}

void ProcessingElement::startWrite(uint32_t delay) {
    state = PEState::WRITING;
    cycles_remaining = delay;
}

void ProcessingElement::setWaiting() {
    state = PEState::WAITING;
    cycles_remaining = 0;
}

void ProcessingElement::setIdle() {
    state = PEState::IDLE;
    cycles_remaining = 0;
}

bool ProcessingElement::tick() {
    bool completed = false;
    
    switch (state) {
        case PEState::IDLE:
            idle_cycles++;
            break;
            
        case PEState::LOADING:
        case PEState::COMPUTING:
        case PEState::WRITING:
            active_cycles++;
            if (cycles_remaining > 0) {
                cycles_remaining--;
            }
            if (cycles_remaining == 0) {
                completed = true;
                state = PEState::IDLE;
            }
            break;
            
        case PEState::WAITING:
            idle_cycles++;  // Waiting counts as idle for utilization
            break;
    }
    
    return completed;
}

bool ProcessingElement::isBusy() const {
    return state != PEState::IDLE;
}

bool ProcessingElement::isComplete() const {
    return state == PEState::IDLE && cycles_remaining == 0;
}
