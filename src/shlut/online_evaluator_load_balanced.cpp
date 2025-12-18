#include "online_evaluator.h"

#include "../utils/helpers.h"

namespace shlut
{
// Helper function to reconstruct shares via king party or direct all-to-all
/*
void OnlineEvaluator::reconstruct(int nP, int pid, std::shared_ptr<io::NetIOMP> network,
                                    const std::vector<Ring>& shares_list,
                                    std::vector<Ring>& reconstructed_list,
                                    bool via_pking, int latency) {
    int pKing = 1;
    size_t num_shares = shares_list.size();
    reconstructed_list.resize(num_shares, 0);
    
    if (via_pking) {
        // Reconstruction via king party
        if (pid != pKing) {
            network->send(pKing, shares_list.data(), shares_list.size() * sizeof(Ring));
            usleep(latency);
            network->recv(pKing, reconstructed_list.data(), reconstructed_list.size() * sizeof(Ring));
        } else {
            std::vector<std::vector<Ring>> share_recv(nP);
            share_recv[pKing - 1] = shares_list;
            usleep(latency);
            
            // Receive from all parties (not parallelized as recv is blocking)
            for (int p = 1; p <= nP; ++p) {
                if (p != pKing) {
                    share_recv[p - 1].resize(num_shares);
                    network->recv(p, share_recv[p - 1].data(), share_recv[p - 1].size() * sizeof(Ring));
                }
            }
            
            // Aggregate shares
            for (int p = 0; p < nP; ++p) {
                for (size_t i = 0; i < num_shares; ++i) {
                    reconstructed_list[i] += share_recv[p][i];
                }
            }
            
            // Send result to all parties (sequential for now to avoid race conditions)
            for (int p = 1; p <= nP; ++p) {
                if (p != pKing) {
                    network->send(p, reconstructed_list.data(), reconstructed_list.size() * sizeof(Ring));
                    network->flush(p);
                }
            }
        }
    } else {
        // Direct reconstruction (all parties exchange shares)
        std::vector<std::vector<Ring>> share_recv(nP);
        share_recv[pid - 1] = shares_list;
        
        // Send to all parties (sequential for now to avoid race conditions)
        for (int p = 1; p <= nP; ++p) {
            if (p != pid) {
                network->send(p, shares_list.data(), shares_list.size() * sizeof(Ring));
            }
        }
        
        usleep(latency);
        
        // Receive from all parties
        for (int p = 1; p <= nP; ++p) {
            if (p != pid) {
                share_recv[p - 1].resize(num_shares);
                network->recv(p, share_recv[p - 1].data(), share_recv[p - 1].size() * sizeof(Ring));
            }
        }
        
        // Aggregate shares
        for (int p = 0; p < nP; ++p) {
            for (size_t i = 0; i < num_shares; ++i) {
                reconstructed_list[i] += share_recv[p][i];
            }
        }
    }
}

OnlineEvaluator::OnlineEvaluator(int nP, int id, std::shared_ptr<io::NetIOMP> network,
                                 PreprocCircuit<Ring> preproc,
                                 common::utils::LevelOrderedCircuit circ,
                                 int threads, int seed, int latency, bool use_pking)
    : nP_(nP),
        id_(id),
        latency_(latency),
        use_pking_(use_pking),
        rgen_(id, seed),
        network_(std::move(network)),
        preproc_(std::move(preproc)),
        circ_(std::move(circ)),
        wires_(circ.num_wires)
    {
        tpool_ = std::make_shared<common::utils::ThreadPool>(threads);
    }

OnlineEvaluator::OnlineEvaluator(int nP, int id, std::shared_ptr<io::NetIOMP> network,
                                 PreprocCircuit<Ring> preproc,
                                 common::utils::LevelOrderedCircuit circ,
                                 std::shared_ptr<common::utils::ThreadPool> tpool, int seed, int latency, bool use_pking)
    : nP_(nP),
        id_(id),
        latency_(latency),
        use_pking_(use_pking),
        rgen_(id, seed),
        network_(std::move(network)),
        preproc_(std::move(preproc)),
        circ_(std::move(circ)),
        tpool_(std::move(tpool)),
        wires_(circ.num_wires) {}

    
        void OnlineEvaluator::setInputs(const std::unordered_map<common::utils::wire_t, Ring> &inputs) {
            // Input gates have depth 0
            for (auto &g : circ_.gates_by_level[0]) {
                if (g->type == common::utils::GateType::kInp) {
                    auto *pre_input = static_cast<PreprocInput<Ring> *>(preproc_.gates[g->out].get());
                    auto pid = pre_input->pid;
                    if (id_ != 0) {
                        if (pid == id_) {
                            Ring accumulated_val = Ring(0);
                            for (size_t i = 1; i <= nP_; i++) {
                                if (i != pid) {
                                    Ring rand_sh;
                                    rgen_.pi(i).random_data(&rand_sh, sizeof(Ring));
                                    accumulated_val += rand_sh;
                                }
                            }
                            wires_[g->out] = inputs.at(g->out) - accumulated_val;
                        } else {
                            rgen_.pi(id_).random_data(&wires_[g->out], sizeof(Ring));
                        }
                    }
                }
            }
        }

        void OnlineEvaluator::setRandomInputs() {
            // Input gates have depth 0.
            for (auto &g : circ_.gates_by_level[0]) {
                if (g->type == common::utils::GateType::kInp) {
                    rgen_.pi(id_).random_data(&wires_[g->out], sizeof(Ring));
                }
            }
        }

        void OnlineEvaluator::multEvaluate(const std::vector<common::utils::FIn2Gate> &mult_gates) {
            if (id_ == 0) { return; }
            size_t num_mult_gates = mult_gates.size();
            std::vector<Ring> shares_to_send(2 * num_mult_gates);

            // Compute masked values: e = x - a, d = y - b
            #pragma omp parallel for
            for (size_t i = 0; i < num_mult_gates; ++i) {
                auto &mult_gate = mult_gates[i];
                auto *pre_mult = static_cast<PreprocMultGate<Ring> *>(preproc_.gates[mult_gate.out].get());
                shares_to_send[2*i] = wires_[mult_gate.in1] - pre_mult->triple_a.valueAt();
                shares_to_send[2*i + 1] = wires_[mult_gate.in2] - pre_mult->triple_b.valueAt();
            }

            // Reconstruct the masked values
            std::vector<Ring> reconstructed(2 * num_mult_gates, 0);
            reconstruct(nP_, id_, network_, shares_to_send, reconstructed, use_pking_, latency_);

            // Compute output shares: z = c + e*b + d*a + e*d (only party 1 adds e*d)
            #pragma omp parallel for
            for (size_t i = 0; i < num_mult_gates; ++i) {
                auto &mult_gate = mult_gates[i];
                auto *pre_mult = static_cast<PreprocMultGate<Ring> *>(preproc_.gates[mult_gate.out].get());
                Ring e = reconstructed[2*i];
                Ring d = reconstructed[2*i + 1];
                Ring share_z = pre_mult->triple_c.valueAt() + e * pre_mult->triple_b.valueAt() + d * pre_mult->triple_a.valueAt();
                if (id_ == 1) {
                    share_z += e * d;
                }
                wires_[mult_gate.out] = share_z;
            }
        }

        void OnlineEvaluator::eqzEvaluate(const std::vector<common::utils::FIn1Gate> &eqz_gates) {
            if (id_ == 0) { return; }
            int pKing = 1; // Designated king party
            size_t num_eqz_gates = eqz_gates.size();
            std::vector<Ring> r1_send(num_eqz_gates);
            std::vector<Ring> r2_send(num_eqz_gates);

            // Compute share of m1 = input + random_value r1
            #pragma omp parallel for
            for (size_t i = 0; i < num_eqz_gates; ++i) {
                auto &eqz_gate = eqz_gates[i];
                auto *pre_eqz = static_cast<PreprocEqzGate<Ring> *>(preproc_.gates[eqz_gate.out].get());
                Ring share_m1 = wires_[eqz_gate.in] + pre_eqz->share_r1.valueAt();
                r1_send[i] = share_m1;
            }

            // Reconstruct the masked input m1
            std::vector<Ring> recon_m1(num_eqz_gates, 0);
            reconstruct(nP_, id_, network_, r1_send, recon_m1, use_pking_, latency_);

            // Compute hamming distance between bits of m1 and bits of r1
            std::vector<Ring> share_m2(num_eqz_gates, 0);
            #pragma omp parallel for
            for (int i = 0; i < num_eqz_gates; ++i) {
                auto *pre_eqz = static_cast<PreprocEqzGate<Ring> *>(preproc_.gates[eqz_gates[i].out].get());
                std::vector<Ring> m1_bits(RINGSIZEBITS);
                m1_bits = bitDecomposeToInt(recon_m1[i]);
                std::vector<Ring> r1_bits(RINGSIZEBITS);
                for (int j = 0; j < RINGSIZEBITS; ++j) {
                    r1_bits[j] = pre_eqz->share_r1_bits[j].valueAt();
                }
                if (id_ == 1) {
                    for (int j = 0; j < RINGSIZEBITS; ++j) {
                        share_m2[i] += m1_bits[j] + r1_bits[j] - 2 * m1_bits[j] * r1_bits[j];
                    }
                    share_m2[i] += pre_eqz->share_r2.valueAt();
                }
                else{
                    for (int j = 0; j < RINGSIZEBITS; ++j) {
                        share_m2[i] += r1_bits[j] - 2 * m1_bits[j] * r1_bits[j];
                    }
                    share_m2[i] += pre_eqz->share_r2.valueAt();
                }
            }


            // Reconstruct the masked input m2
            std::vector<Ring> recon_m2(num_eqz_gates, 0);
            reconstruct(nP_, id_, network_, share_m2, recon_m2, use_pking_, latency_);
    

            // Compute final output share
            std::vector<Ring> recon_out(num_eqz_gates, 0);
            #pragma omp parallel for
            for (int i = 0; i < num_eqz_gates; ++i) {
                auto *pre_eqz = static_cast<PreprocEqzGate<Ring> *>(preproc_.gates[eqz_gates[i].out].get());
                std::vector<Ring> r2_bits(RINGSIZEBITS);
                for (int j = 0; j < RINGSIZEBITS; ++j) {
                    r2_bits[j] = pre_eqz->share_r2_bits[j].valueAt();
                }
                recon_out[i] = r2_bits[recon_m2[i]%RINGSIZEBITS];
                wires_[eqz_gates[i].out] = recon_out[i]; // Reconstructed output
            }
        }

        void OnlineEvaluator::recEvaluate(const std::vector<common::utils::FIn1Gate> &rec_gates) {
            if (id_ == 0) { return; }
            size_t num_rec_gates = rec_gates.size();
            std::vector<Ring> shares_to_send(num_rec_gates);

            // Gather shares to reconstruct
            #pragma omp parallel for
            for (size_t i = 0; i < num_rec_gates; ++i) {
                auto &rec_gate = rec_gates[i];
                shares_to_send[i] = wires_[rec_gate.in];
            }

            // Reconstruct the values
            std::vector<Ring> reconstructed(num_rec_gates, 0);
            
            reconstruct(nP_, id_, network_, shares_to_send, reconstructed, use_pking_, latency_);

            // Store reconstructed values in output wires
            #pragma omp parallel for
            for (size_t i = 0; i < num_rec_gates; ++i) {
                wires_[rec_gates[i].out] = reconstructed[i];
            }
        }

            std::vector<Ring> OnlineEvaluator::getOutputs() {
        std::vector<Ring> outvals(circ_.outputs.size());
        if (circ_.outputs.empty()) {
            return outvals;
        }
        if (id_ != 0) {
            std::vector<std::vector<Ring>> output_shares(nP_, std::vector<Ring>(circ_.outputs.size()));
            for (size_t i = 0; i < circ_.outputs.size(); ++i) {
                auto wout = circ_.outputs[i];
                output_shares[id_ - 1][i] = wires_[wout];
            }
            for (int pid = 1; pid <= nP_; ++pid) {
                if (pid != id_) {
                    network_->send(pid, output_shares[id_ - 1].data(), output_shares[id_ - 1].size() * sizeof(Ring));
                }
            }
            usleep(latency_);
            // #pragma omp parallel for
            for (int pid = 1; pid <= nP_; ++pid) {
                if (pid != id_) {
                    network_->recv(pid, output_shares[pid - 1].data(), output_shares[pid - 1].size() * sizeof(Ring));
                }
            }
            for (size_t i = 0; i < circ_.outputs.size(); ++i) {
                Ring outmask = Ring(0);
                for (int pid = 1; pid <= nP_; ++pid) {
                    outmask += output_shares[pid - 1][i];
                }
                outvals[i] = outmask;
            }
        }
        return outvals;
    }

    std::vector<Ring> OnlineEvaluator::evaluateCircuit(const std::unordered_map<common::utils::wire_t, Ring> &inputs) {
        setInputs(inputs);
        for (size_t i = 0; i < circ_.gates_by_level.size(); ++i) {
            evaluateGatesAtDepth(i);
        }
        return getOutputs();
    }


void OnlineEvaluator::evaluateGatesAtDepth(size_t depth) {
    
std::vector<common::utils::FIn2Gate> mult_gates;
std::vector<common::utils::FIn2Gate> mult_gates_RSS;
std::vector<common::utils::FIn1Gate> eqz_gates;
std::vector<common::utils::FIn1Gate> rec_gates;
std::vector<common::utils::FIn1Gate> rec_gates_RSS;

for (auto &gate : circ_.gates_by_level[depth]) {
    switch (gate->type) {
        
        case common::utils::GateType::kMul: {
            auto *g = static_cast<common::utils::FIn2Gate *>(gate.get());
            mult_gates.push_back(*g);
            break;
        }
        
        case common::utils::GateType::kEqz: {
            auto *g = static_cast<common::utils::FIn1Gate *>(gate.get());
            eqz_gates.push_back(*g);
            break;
        }
        case common::utils::GateType::kRec: {
            auto *g = static_cast<common::utils::FIn1Gate *>(gate.get());
            rec_gates.push_back(*g);
            break;
        }
        
       
       case common::utils::GateType::kRecRSS: {
           auto *g = static_cast<common::utils::FIn1Gate *>(gate.get());
           rec_gates_RSS.push_back(*g);
           break;
        }
        case common::utils::GateType::kMulRSS: {
            auto *g = static_cast<common::utils::FIn2Gate *>(gate.get());
            mult_gates_RSS.push_back(*g);
            break;
        }
        default:
        break;
    }
}
if (!mult_gates.empty()) { multEvaluate(mult_gates); }
if (!eqz_gates.empty()) { eqzEvaluate(eqz_gates); }
if (!rec_gates.empty()) { recEvaluate(rec_gates); }

if (!mult_gates_RSS.empty()) { multEvaluateRSS(mult_gates_RSS); }
// if (!rec_gates_RSS.empty()) { recEvaluateRSS(rec_gates_RSS); }

// Second pass: handle locally evaluable gates.
for (auto &gate : circ_.gates_by_level[depth]) {
    switch (gate->type) {
        case common::utils::GateType::kAdd:
        case common::utils::GateType::kAddRSS: {
            auto *g = static_cast<common::utils::FIn2Gate *>(gate.get());
            wires_[g->out] = wires_[g->in1] + wires_[g->in2];
            break;
        }
        case common::utils::GateType::kSub:
        case common::utils::GateType::kSubRSS: {
            auto *g = static_cast<common::utils::FIn2Gate *>(gate.get());
            wires_[g->out] = wires_[g->in1] - wires_[g->in2];
            break;
        }
        case common::utils::GateType::kConstAdd:
        case common::utils::GateType::kConstAddRSS: {
            auto *g = static_cast<common::utils::ConstOpGate *>(gate.get());
            if (id_ == 1) {
                wires_[g->out] = wires_[g->in] + g->cval;
            } else {
                wires_[g->out] = wires_[g->in];
            }
            break;
        }
        case common::utils::GateType::kConstMul:
        case common::utils::GateType::kConstMulRSS: {
            auto *g = static_cast<common::utils::ConstOpGate *>(gate.get());
            wires_[g->out] = wires_[g->in] * g->cval;
            break;
        }
        default:
        break;
    }
}
}

*/

// Make function which will reconstruct in parallel.

void OnlineEvaluator::multipleAbortReconstructRSS(int pid, std::shared_ptr<io::NetIOMP> network,
                                     std::vector<PartyShare<Ring>>& shares_list, 
                                     int latency, std::vector<Ring>& reconstructed_list){
    size_t num_shares = shares_list.size();
    reconstructed_list.resize(num_shares, 0);

    std::vector<Ring> commonvalue1_list(num_shares, 0);
    std::vector<Ring> commonvalue2_list(num_shares, 0);
    std::vector<Ring> other_commonvalue1_list(num_shares, 0);
    std::vector<Ring> other_commonvalue2_list(num_shares, 0);
    Ring send_value1{}, send_value2{};

    // Send both shares to all parties
    for (int p = 0; p < 3; ++p) {
        if (p != pid) {
            for (size_t i = 0; i < num_shares; ++i) {
                send_value1 = shares_list[i].value1();
                send_value2 = shares_list[i].value2();
                network->send(p, &send_value1, sizeof(Ring));
                network->send(p, &send_value2, sizeof(Ring));
            }
        }
    }
    usleep(latency);
    // Receive from all parties
    for (size_t i = 0; i < num_shares; ++i) {
        network->recv((pid + 1) % 3, &commonvalue1_list[i], sizeof(Ring));
        network->recv((pid + 1) % 3, &other_commonvalue1_list[i], sizeof(Ring));
        network->recv((pid + 2) % 3, &other_commonvalue2_list[i], sizeof(Ring));
        network->recv((pid + 2) % 3, &commonvalue2_list[i], sizeof(Ring));
    }

    // Checks and reconstruction
    for (size_t i = 0; i < num_shares; ++i) {
        if (commonvalue1_list[i] != shares_list[i].value2() || 
            commonvalue2_list[i] != shares_list[i].value1() || 
            other_commonvalue1_list[i] != other_commonvalue2_list[i]){
            throw std::runtime_error("Abort reconstruction failed: inconsistent shares received.");
        }
        reconstructed_list[i] = shares_list[i].value1() + shares_list[i].value2() + commonvalue1_list[i];
    }
} //Check this function once                

OnlineEvaluator::OnlineEvaluator(int id, std::shared_ptr<io::NetIOMP> network,
                                     PreprocCircuit<Ring> preproc,
                                     common::utils::LevelOrderedCircuit circ,
                                     int threads, int seed, int latency, bool use_pking)
    : id_(id),
        latency_(latency),
        rgen_(id, seed),
        use_pking_(use_pking),
        network_(std::move(network)),
        preproc_(std::move(preproc)),
        circ_(std::move(circ)),
        wires_(circ.num_wires)
    {
        tpool_ = std::make_shared<common::utils::ThreadPool>(threads);
    }

OnlineEvaluator::OnlineEvaluator(int id, std::shared_ptr<io::NetIOMP> network,
                                     PreprocCircuit<Ring> preproc,
                                     common::utils::LevelOrderedCircuit circ,
                                     std::shared_ptr<common::utils::ThreadPool> tpool, int seed, int latency, bool use_pking)
    : id_(id),
        latency_(latency),
        rgen_(id, seed),
        use_pking_(use_pking),
        network_(std::move(network)),
        preproc_(std::move(preproc)),
        circ_(std::move(circ)),
        tpool_(std::move(tpool)),
        wires_(circ.num_wires) {}


        void OnlineEvaluator::evaluateGatesAtDepth(size_t depth) {
            std::vector<common::utils::FIn2Gate> mult_gates_additive;
            std::vector<common::utils::FIn1Gate> add2rss_gates;
            // Probably should have something for Augmented Sharing as well.
            
            for (auto &gate : circ_.gates_by_level[depth]){
                switch (gate->type) {
                    case common::utils::GateType::kMul: {
                        auto *g = static_cast<common::utils::FIn2Gate *>(gate.get());
                        mult_gates_additive.push_back(*g);
                        break;
                    }

                    case common::utils::GateType::kAdd2RSS: {
                        auto *g = static_cast<common::utils::FIn1Gate *>(gate.get());
                        add2rss_gates.push_back(*g);
                        break;
                    }
                    default:
                        break;
            }
        }

        if (!mult_gates_additive.empty()) { multEvaluateAdditive(mult_gates_additive); }
        if (!add2rss_gates.empty()) { additive2RSS(add2rss_gates); }

        // Second pass: handle locally evaluable gates.
        for (auto &gate : circ_.gates_by_level[depth]) {
            switch (gate->type) {
                case common::utils::GateType::kAdd:{
                    auto *g = static_cast<common::utils::FIn2Gate *>(gate.get());
                    wires_[g->out] = {wires_[g->in1][0] + wires_[g->in2][0]};
                    break;
                }

                case common::utils::GateType::kAddRSS:{
                    auto *g = static_cast<common::utils::FIn2Gate *>(gate.get());
                    wires_[g->out] = {wires_[g->in1][0] + wires_[g->in2][0], 
                                     wires_[g->in1][1] + wires_[g->in2][1]};
                    break;
                }

                case common::utils::GateType::kSub:{
                    auto *g = static_cast<common::utils::FIn2Gate *>(gate.get());
                    wires_[g->out] = {wires_[g->in1][0] - wires_[g->in2][0]};
                    break;
                }

                case common::utils::GateType::kSubRSS:{
                    auto *g = static_cast<common::utils::FIn2Gate *>(gate.get());
                    wires_[g->out] = {wires_[g->in1][0] - wires_[g->in2][0], 
                                     wires_[g->in1][1] - wires_[g->in2][1]};
                    break;
                }

                case common::utils::GateType::kConstAdd:{
                    auto *g = static_cast<common::utils::ConstOpGate *>(gate.get());
                    if (id_ == 1) {
                        wires_[g->out] = {wires_[g->in][0] + g->cval};
                    } else {
                        wires_[g->out] = {wires_[g->in][0]};
                    }
                    break;
                }

                case common::utils::GateType::kConstAddRSS:{
                    auto *g = static_cast<common::utils::ConstOpGate *>(gate.get());
                    if (id_ == 1) {
                        wires_[g->out] = {wires_[g->in][0] + g->cval, wires_[g->in][1]};
                    } else if (id_ == 2) {
                        wires_[g->out] = {wires_[g->in][0], wires_[g->in][1] + g->cval};
                    } else {
                        wires_[g->out] = {wires_[g->in][0], wires_[g->in][1]};
                    }
                    break;
                }

                case common::utils::GateType::kConstMul:{
                    auto *g = static_cast<common::utils::ConstOpGate *>(gate.get());
                    wires_[g->out] = {wires_[g->in][0] * g->cval};
                    break;
                }

                case common::utils::GateType::kConstMulRSS:{
                    auto *g = static_cast<common::utils::ConstOpGate *>(gate.get());
                    wires_[g->out] = {wires_[g->in][0] * g->cval, wires_[g->in][1] * g->cval};
                    break;
                }

                case common::utils::GateType::kMulRSS:{
                    auto *g = static_cast<common::utils::FIn2Gate *>(gate.get());
                    wires_[g->out] = {wires_[g->in1][0] * wires_[g->in2][0] + wires_[g->in1][0] * wires_[g->in2][1] + wires_[g->in1][1] * wires_[g->in2][0]}
                }
                default:
                    break;
            }
        }
    }



        void OnlineEvaluator::multiReconstructAdditive(int pid, std::shared_ptr<io::NetIOMP> network,
                                     std::vector<Ring>& shares_list, 
                                     int latency, std::vector<Ring>& reconstructed_list){
            size_t num_shares = shares_list.size();
            reconstructed_list.resize(num_shares, 0);
            Ring send_value{};
            // Send shares to all parties
            for (int p = 1; p <= 3; ++p) {
                if (p != pid) {
                    for (size_t i = 0; i < num_shares; ++i) {
                        send_value = shares_list[i];
                        network->send(p, &send_value, sizeof(Ring));
                    }
                }
            }
            usleep(latency);
            // Receive from all parties
            for (size_t i = 0; i < num_shares; ++i) {
                Ring recv_value1{}, recv_value2{};
                network->recv((pid + 1) % 3, &recv_value1, sizeof(Ring));
                network->recv((pid + 2) % 3, &recv_value2, sizeof(Ring));
                reconstructed_list[i] = shares_list[i] + recv_value1 + recv_value2;
            }
        }≈

        std::vector<Ring> OnlineEvaluator::getOutputsRSS(){
            std::vector<Ring> outvals(circ_.outputs.size());
            if (circ_.outputs.empty()) {
                return outvals;
            }
            std::vector<PartyShare<Ring>> output_shares(circ_.outputs.size());
            for (size_t i = 0; i < circ_.outputs.size(); ++i) {
                auto wout = circ_.outputs[i];
                output_shares[i] = wires_[wout];
            }
            multipleAbortReconstructRSS(id_, network_, output_shares, latency_, outvals);

            return outvals;
        }

        void OnlineEvaluator::setInputsAdditive(const std::unordered_map<common::utils::wire_t, Ring> &inputs)
        {
            for(auto &g : circ_.gates_by_level[0]){
                if(g->type == common::utils::GateType::kInpAdd){
                    auto *pre_input = static_cast<PreprocInputAdd<Ring> *>(preproc_.gates[g->out].get());
                    auto pid = pre_input->pid;
                    if (pid == id_){
                        Ring value_ = Ring(0);
                        std::vector<Ring> val = {inputs.at(g->out)};
                        wires_[g->out] = inputs.at(g->out);
                        rgen_.pi((pid + 1) % 3).random_data(&value_, sizeof(Ring));
                        val[0] = val[0] - value_;
                        rgen_.pi((pid + 2) % 3).random_data(&value_, sizeof(Ring));
                        val[0] = val[0] - value_;
                        // Make a vector of size 1 and keep it as wires_[g->out]
                        wires_[g->out] = val;
                        
                    } else {
                        Ring value_ = Ring(0);
                        rgen_.pi(pid).random_data(&value_, sizeof(Ring));
                        wires_[g->out] = {value_};
                    }
                }
            }
        }


        void OnlineEvaluator::setInputsRSS(const std::unordered_map<common::utils::wire_t, Ring> &inputs) {
            // Input gates have depth 0
            for (auto &g : circ_.gates_by_level[0]) {
                if (g->type == common::utils::GateType::kInpRSS) {
                    auto *pre_input = static_cast<PreprocInput<Ring> *>(preproc_.gates[g->out].get());
                    auto pid = pre_input->pid;
                    if (pid == id_){
                        Ring value1_ = Ring(0);
                        Ring value2_ = Ring(0);
                        rgen_.pi((pid + 2) % 3).random_data(&value1_, sizeof(Ring));
                        rgen_.all().random_data(&value2_, sizeof(Ring));
                        value2_ = inputs.at(g->out) - value1_ - value2_;
                        // Send value2_ to party (pid + 1) % 3
                        network_->send((pid + 1) % 3, &value2_, sizeof(Ring));
                        wires_[g->out] = {value1_, value2_};
                    } else {
                        Ring value1_ = Ring(0);
                        Ring value2_ = Ring(0);
                        if (id_ == (pid + 2) % 3){
                            rgen_.pi((pid + 2) % 3).random_data(&value2_, sizeof(Ring));
                            rgen_.all().random_data(&value1_, sizeof(Ring));
                        } else {
                            network_->recv(pid, &value1_, sizeof(Ring));
                            rgen_.all().random_data(&value2_, sizeof(Ring));
                        }
                        wires_[g->out] = {value1_, value2_};
                    }
                }
            }        
        }

        void OnlineEvaluator::setInputsAug(const std::unordered_map<common::utils::wire_t, Ring> &inputs) {

        }
        
        void OnlineEvaluator::setRandomInputsAdditive() {
            for (auto &g : circ_.gates_by_level[0]) {
                if (g->type == common::utils::GateType::kInpAdd) {
                    Ring value_ = Ring(0);
                    rgen_.pi(id_).random_data(&value_, sizeof(Ring));
                    wires_[g->out] = {value_};
                }
            }
        }

        void OnlineEvaluator::setRandomInputsRSS() {
            // Input gates have depth 0.
            for (auto &g : circ_.gates_by_level[0]) {
                if (g->type == common::utils::GateType::kInp) {
                    Ring value1_ = Ring(0);
                    Ring value2_ = Ring(0);
                    rgen_.pi((id_ + 2) % 3).random_data(&value1_, sizeof(Ring));
                    rgen_.pi((id_ + 1) % 3).random_data(&value2_, sizeof(Ring));
                    wires_[g->out] = {value1_, value2_};
                }
            }
        }

        void OnlineEvaluator::setRandomInputsAug() {
            for (auto &g : circ_.gates_by_level[0]) {
                if (g->type == common::utils::GateType::kInpAug) {
                    Ring value1_ = Ring(0);
                    Ring value2_ = Ring(0);
                    Ring value3_ = Ring(0);
                    rgen_.all().random_data(&value1_, sizeof(Ring));
                    rgen_.pi((id + 1) % 3).random_data(&value1_, sizeof(Ring));
                    rgen_.pi((id + 2) % 3).random_data(&value1_, sizeof(Ring));
                    wires_[g->out] = {value1_, value2_, value3_};
                }
            }
        }

        void OnlineEvaluator::multEvaluateAdditive(const std::vector<common::utils::FIn2Gate> &mult_gates_additive) {
            size_t num_mult_gates = mult_gates_additive.size();
            std::vector<Ring> shares_to_send(2 * num_mult_gates);

            // #pragma omp parallel for
            for(size_t i = 0; i < num_mult_gates; ++i){
                auto &mult_gate = mult_gates_additive[i];
                auto *pre_mult = static_cast<PreprocMultGateAdd<Ring> *>(preproc_.gates[mult_gate.out].get()); // Need to change the name here and in preproc.h for all sharings
                
                Ring x_share = wires_[mult_gate.in1][0];
                Ring y_share = wires_[mult_gate.in2][0];
                
                shares_to_send[2*i] = x_share - pre_mult->triple_a.valueAt();
                shares_to_send[2*i + 1] = y_share - pre_mult->triple_b.valueAt();
            }

            std::vector<Ring> reconstructed(2 * num_mult_gates, 0);
            multiReconstructAdditive(id_, network_, shares_to_send, latency_, reconstructed);

            // #pragma omp parallel for
            for(size_t i = 0; i < num_mult_gates; ++i){
                auto &mult_gate = mult_gates_additive[i];
                auto *pre_mult = static_cast<PreprocMultGateAdd<Ring> *>(preproc_.gates[mult_gate.out].get());
                
                Ring e = reconstructed[2*i];
                Ring d = reconstructed[2*i + 1];
                Ring share_z = pre_mult->triple_c.valueAt() + e * pre_mult->triple_b.valueAt() + d * pre_mult->triple_a.valueAt();
                if (id_ == 1) {
                    share_z += e * d;
                }
                wires_[mult_gate.out] = {share_z};
            }
        }

        void OnlineEvaluator::additive2RSS(const std::vector<common::utils::FIn1Gate> &additive2RSS_gates) {
            size_t num_gates = additive2RSS_gates.size();

            // #pragma omp parallel for
            for (size_t i = 0; i < num_gates; ++i){
                auto &gate = additive2RSS_gates[i];
                Ring share = wires_[gate.in][0];
                // Share the value into RSS
                Ring value1_ = Ring(0);
                Ring value2_ = Ring(0);
                rgen_.pi((id_ + 2) % 3).random_data(&value1_, sizeof(Ring));
                rgen_.all().random_data(&value2_, sizeof(Ring));
                value2_ = share - value1_ - value2_;
                // Send value2_ to party (id_ + 1) % 3
                network_->send((id_ + 1) % 3, &value2_, sizeof(Ring));
                wires_[gate.out] = {value1_, value2_};

                //Receive from other parties
                usleep(latency_);
                Ring other_value1_ = Ring(0);
                Ring other_value2_ = Ring(0);
                rgen_.pi((id_ + 1) % 3).random_data(&other_value2_, sizeof(Ring));
                rgen_.all().random_data(&other_value1_, sizeof(Ring));
                wires_[gate.out][0] += other_value1_;
                wires_[gate.out][1] += other_value2_;
                
                usleep(latency_);
                other_value1_ = Ring(0);
                other_value2_ = Ring(0);
                network_->recv((id_ + 2) % 3, &other_value1_, sizeof(Ring));
                rgen_.all().random_data(&other_value2_, sizeof(Ring));
                wires_[gate.out][0] += other_value1_;
                wires_[gate.out][1] += other_value2_;
            }
        }



        /*
        void OnlineEvaluator::publicLUTSemiRSS(int dealer_pid, const std::vector<Ring>& lut, 
                                            const std::vector<common::utils::FIn1Gate>& lut_gates){
            size_t num_lut_gates = lut_gates.size();

            #pragma omp parallel for
            for (size_t i = 0; i < num_lut_gates; ++i){
                auto &lut_gate = lut_gates[i];
                PartyShare<Ring> input_share = wires_[lut_gate.in];
                int pid = input_share.pid();
                wires_[lut_gate.out] = PartyShare<Ring>(pid, Ring(0), Ring(0)); // Placeholder
                Ring value1_ = Ring(0);
                Ring value2_ = Ring(0);

                if(pid == dealer_pid){
                    Ring r = (input_share.value1() + input_share.value2());
                    K_a, K_b = DPFGen(r, lut.size()) // Need to implement DPFGen
                    network_->send((dealer_pid + 1) % 3, K_b);
                    network_->send((dealer_pid + 2) % 3, K_c);

                    value1_ = Ring(0);
                    value2_ = Ring(0);
                    rgen_.pi((pid + 1) % 3).random_data(&value1_, sizeof(Ring));
                    rgen_.all().random_data(&value2_, sizeof(Ring));
                    wires_[lut_gate.out] += PartyShare<Ring>(pid, value1_, value2_);

                    value1_ = Ring(0);
                    value2_ = Ring(0);
                    network_->recv((pid + 2) % 3, &value1_, sizeof(Ring));
                    rgen_.all().random_data(&value2_, sizeof(Ring));
                    wires_[lut_gate.out] += PartyShare<Ring>(pid, value1_, value2_);
                }
                else{
                    DPFKey K;
                    network_->recv(dealer_pid, &K);
                    std::vector<Ring> lut_shares = EvalAll(K, lut.size()); // Need to implement EvalAll
                    int val = (dealer_pid - pid + 3) % 3;
                    Ring num_rotate = val == 1 ? input_share.value1() : input_share.value2();
                    // Rotate lut_shares by num_rotate
                    std::vector<Ring> rotated_lut_shares(lut.size());
                    for (size_t j = 0; j < lut.size(); ++j){
                        rotated_lut_shares[j] = lut_shares[(j + num_rotate) % lut.size()];
                    }
                    // Take the inner product of rotated_lut_shares and lut
                    Ring share_lut = Ring(0);
                    for (size_t j = 0; j < lut.size(); ++j){
                        share_lut += rotated_lut_shares[j] * lut[j];
                    }
                    
                    Ring value1_ = Ring(0);
                    Ring value2_ = Ring(0);
                    rgen_.pi((pid + 2) % 3).random_data(&value1_, sizeof(Ring));
                    rgen_.all().random_data(&value2_, sizeof(Ring));
                    value2_ = share_lut - value1_ - value2_;
                    // Send value2_ to party (pid + 1) % 3
                    network_->send((pid + 1) % 3, &value2_, sizeof(Ring));
                    wires_[lut_gate.out] += PartyShare<Ring>(pid, value1_, value2_);

                    if ((pid + 1) % 3 != dealer_pid){
                        value1_ = Ring(0);
                        value2_ = Ring(0);
                        rgen_.pi((pid + 1) % 3).random_data(&value2_, sizeof(Ring));
                        rgen_.all().random_data(&value1_, sizeof(Ring));
                        wires_[lut_gate.out] += PartyShare<Ring>(pid, value1_, value2_);
                    }

                    if ((pid + 2) % 3 != dealer_pid){
                        value1_ = Ring(0);
                        value2_ = Ring(0);
                        network_->recv((pid + 2) % 3, &value1_, sizeof(Ring));
                        rgen_.all().random_data(&value2_, sizeof(Ring));
                        wires_[lut_gate.out] += PartyShare<Ring>(pid, value1_, value2_);
                    }
                }
            }
        }
        
        void OnlineEvaluator::privateLUTSemiRSS(const std::vector<PartyShare<Ring>>& lut_shares, 
                                             const std::vector<common::utils::FIn1Gate>& lut_gates){
            // To be done later after some clarifications.
        }
        */
}; // namespace shlut

// NOTE TO SELF: Check the multevaluate and publicLUTSemi in the sharing part.

