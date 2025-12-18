#pragma once

#include <memory>
#include <unordered_map>
#include <vector>

#include "../io/netmp.h"
#include "../utils/circuit.h"
#include "../utils/thread_pool.h"
#include "preproc.h"
#include "rand_gen_pool.h"
#include "sharing.h"
#include "../utils/types.h"

using namespace common::utils;

namespace shlut {
  class OnlineEvaluator {
    int id_;
    int latency_;  // Network latency in microseconds
    bool use_pking_;  // Use king party for reconstruction
    RandGenPool rgen_;
    std::shared_ptr<io::NetIOMP> network_;
    PreprocCircuit<Ring> preproc_;
    common::utils::LevelOrderedCircuit circ_;
    std::vector<PartyShare<Ring>> wires_;
    std::shared_ptr<common::utils::ThreadPool> tpool_;

    // Helper function to reconstruct shares via king party or direct all-to-all
    /* static void reconstruct(int nP, int pid, std::shared_ptr<io::NetIOMP> network,
                           const std::vector<Ring>& shares_list,
                           std::vector<Ring>& reconstructed_list,
                           bool via_pking, int latency); */
    


  public:
    /* OnlineEvaluator(int nP, int id, std::shared_ptr<io::NetIOMP> network,
                    PreprocCircuit<Ring> preproc,
                    common::utils::LevelOrderedCircuit circ,
                    int threads, int seed = 200, int latency = 100, bool use_pking = true);

    OnlineEvaluator(int nP, int id, std::shared_ptr<io::NetIOMP> network,
                    PreprocCircuit<Ring> preproc,
                    common::utils::LevelOrderedCircuit circ,
                    std::shared_ptr<common::utils::ThreadPool> tpool, int seed = 200, int latency = 100, bool use_pking = true); */
    
    OnlineEvaluator(int pid, std::shared_ptr<io::NetIOMP> network,
                    PreprocCircuit<Ring> preproc,
                    common::utils::LevelOrderedCircuit circ,
                    int threads, int seed = 200, int latency = 100, bool use_pking = true);
      
    OnlineEvaluator(int pid, std::shared_ptr<io::NetIOMP> network, 
                    PreprocCircuit<Ring> preproc,
                    common::utils::LevelOrderedCircuit circ,
                    std::shared_ptr<common::utils::ThreadPool> tpool, int seed = 200, int latency = 100, bool use_pking = true);

    /* void setInputs(const std::unordered_map<common::utils::wire_t, Ring> &inputs);

    void setRandomInputs();

    void multEvaluate(const std::vector<common::utils::FIn2Gate> &mult_gates);

    void eqzEvaluate(const std::vector<common::utils::FIn1Gate> &eqz_gates);

    void recEvaluate(const std::vector<common::utils::FIn1Gate> &rec_gates);
    
    Ring reconstruct(AddShare<Ring> &shares);
    
    std::vector<Ring> evaluateCircuit(const std::unordered_map<common::utils::wire_t, Ring> &inputs);

    std::vector<Ring> getOutputs();
 */
    void evaluateGatesAtDepth(size_t depth);

    void multiReconstructAdditive(int pid, std::shared_ptr<io::NetIOMP> network,
      std::vector<Ring>& shares_list, 
      int latency, std::vector<Ring>& reconstructed_list);

    std::vector<Ring> getOutputsRSS();

    void setInputsAdditive(const std::unordered_map<common::utils::wire_t, Ring> &inputs);

    void setInputsRSS(const std::unordered_map<common::utils::wire_t, Ring> &inputs);

    void setInputsAug(const std::unordered_map<common::utils::wire_t, Ring> &inputs);

    void setRandomInputsAdditive();

    void setRandomInputsRSS();

    void setRandomInputsAug();

    void multEvaluateAdditive(const std::vector<common::utils::FIn2Gate> &mult_gates_additive);

    void additive2RSS(const std::vector<common::utils::FIn1Gate> &additive2RSS_gates);

    void multipleAbortReconstructRSS(int pid, std::shared_ptr<io::NetIOMP> network,
      std::vector<PartyShare<Ring>>& shares_list, 
      int latency, std::vector<Ring>& reconstructed_list);
      
      /*     void publicLUTSemiRSS(int dealer_pid, const std::vector<Ring> &lut, const std::vector<common::utils::FIn1Gate> &lut_gates);
      
          void privateLUTSemiRSS(const std::vector<PartyShare<Ring>> &lut_shares, const std::vector<common::utils::FIn1Gate> &lut_gates); */


  };

}; // namespace shlut
