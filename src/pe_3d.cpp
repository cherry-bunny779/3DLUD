#include "pe_3d.hpp"

ProcessingElement3D::ProcessingElement3D(uint32_t r, uint32_t c, uint32_t z)
    : ProcessingElement(r, c)
    , layer(z)
    , tsv_sending(false)
    , tsv_receiving(false)
    , tsv_cycles_remaining(0)
    , tsv_send_data(0.0f)
    , tsv_recv_data(0.0f)
    , tsv_active_cycles(0)
{}

void ProcessingElement3D::reset() {
    ProcessingElement::reset();
    layer = 0;
    tsv_sending = false;
    tsv_receiving = false;
    tsv_cycles_remaining = 0;
    tsv_send_data = 0.0f;
    tsv_recv_data = 0.0f;
    tsv_active_cycles = 0;
}

void ProcessingElement3D::startTSVSend(float data, uint32_t latency) {
    tsv_send_data = data;
    tsv_sending = true;
    tsv_receiving = false;
    tsv_cycles_remaining = latency;
}

void ProcessingElement3D::startTSVReceive(uint32_t latency) {
    tsv_sending = false;
    tsv_receiving = true;
    tsv_cycles_remaining = latency;
}

void ProcessingElement3D::setTSVReceivedData(float data) {
    tsv_recv_data = data;
}

bool ProcessingElement3D::tick() {
    // First handle base PE tick
    bool base_completed = ProcessingElement::tick();
    
    // Handle TSV operations
    if (tsv_sending || tsv_receiving) {
        tsv_active_cycles++;
        if (tsv_cycles_remaining > 0) {
            tsv_cycles_remaining--;
        }
        if (tsv_cycles_remaining == 0) {
            tsv_sending = false;
            tsv_receiving = false;
        }
    }
    
    return base_completed;
}

bool ProcessingElement3D::isBusy() const {
    return ProcessingElement::isBusy() || tsv_sending || tsv_receiving;
}
