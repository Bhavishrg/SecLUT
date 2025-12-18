#pragma once

#include <emp-tool/emp-tool.h>

#include <algorithm>
#include <memory>
#include <string>
#include <unordered_map>

#include "../io/netmp.h"
#include "../utils/circuit.h"
#include "../utils/thread_pool.h"


#include "preproc.h"
#include "shlut/rand_gen_pool.h"
#include "sharing.h"
#include "../utils/types.h"

using namespace common::utils;

namespace shlut {
class OfflineEvaluator {
  int id_;
  int latency_;  // Network latency in microseconds
  bool use_pking_;  // Use king party for reconstruction
  RandGenPool rgen_;
  std::shared_ptr<io::NetIOMP> network_;
  common::utils::LevelOrderedCircuit circ_;
  std::shared_ptr<common::utils::ThreadPool> tpool_;
  PreprocCircuit<Ring> preproc_;

  // Used for running common coin protocol. Returns common random PRG key which
  // is then used to generate randomness for common coin output.
  //emp::block commonCoinKey();

  // Helper function to reconstruct shares via king party or direct all-to-all
  static void reconstruct(int nP, int pid, std::shared_ptr<io::NetIOMP> network,
                         const std::vector<Ring>& shares_list, 
                         std::vector<Ring>& reconstructed_list,
                         bool via_pking, int latency);

  public:
  
  OfflineEvaluator(int my_id, std::shared_ptr<io::NetIOMP> network, PreprocCircuit<Ring> preproc,
                   common::utils::LevelOrderedCircuit circ, int threads, int seed = 200, int latency = 100, bool use_pking = true);

  // Generate sharing of a random unknown value.
  static void randomAdditiveShare(RandGenPool& rgen, Share<Ring, 1>& share);

  static void randomReplicatedShare(RandGenPool& rgen, Share<Ring, 2>& share);

  static void randomAugmentedShare(RandGenPool& rgen, Share<Ring, 3>& share);
  // Generate sharing of a random value known to party. Should be called by
  // dealer when other parties call other variant.
  static void AdditiveShareSecret(RandGenPool& rgen, Share<Ring, 1>& share, Ring secret, int pid);
  
  static void ReplicatedShareSecret(RandGenPool& rgen, Share<Ring, 2>& share, Ring secret, int pid);

  static void AugmentedShareSecret(RandGenPool& rgen, Share<Ring, 3>& share, Ring secret, int pid);
  // Following methods implement various preprocessing subprotocols.

  // Set masks for each wire. Should be called before running any of the other
  // subprotocols.
  void setWireMasksParty(const std::unordered_map<common::utils::wire_t, int>& input_pid_map,
                         std::vector<Ring>& rand_sh_sec);

  // void setWireMasks(const std::unordered_map<common::utils::wire_t, int>& input_pid_map);

  PreprocCircuit<Ring> getPreproc();

  // Efficiently runs above subprotocols.
  PreprocCircuit<Ring> run(const std::unordered_map<common::utils::wire_t, int>& input_pid_map);
};

};  // namespace shlut
