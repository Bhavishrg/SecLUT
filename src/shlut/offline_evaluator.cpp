#include "offline_evaluator.h"

#include <NTL/BasicThreadPool.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <thread>

// #include "../utils/helpers.h"

namespace shlut {
OfflineEvaluator::OfflineEvaluator(int nP, int my_id,
                                   std::shared_ptr<io::NetIOMP> network,
                                   common::utils::LevelOrderedCircuit circ,
                                   int threads, int seed, int latency, bool use_pking)
    : nP_(nP),
      id_(my_id),
      latency_(latency),
      use_pking_(use_pking),
      rgen_(my_id, seed), 
      network_(std::move(network)),
      circ_(std::move(circ))
      // preproc_(circ.num_gates)

      { } // tpool_ = std::make_shared<ThreadPool>(threads); }




    void OfflineEvaluator::randomShare(int nP, int pid, RandGenPool& rgen, AddShare<Ring>& share, TPShare<Ring>& tpShare) {
      Ring val = Ring(0);
      if (pid == 0) {
        share.pushValue(Ring(0));
        tpShare.pushValues(Ring(0));
        for (int i = 1; i <= nP; i++) {
          rgen.pi(i).random_data(&val, sizeof(Ring));
          tpShare.pushValues(val);
        }
      } else {
        rgen.p0().random_data(&val, sizeof(Ring));
        share.pushValue(val);
      }
    }

    void OfflineEvaluator::randomShareSecret(int nP, int pid, RandGenPool& rgen,
                                            AddShare<Ring>& share, TPShare<Ring>& tpShare, Ring secret,
                                            std::vector<Ring>& rand_sh_sec, size_t& idx_rand_sh_sec) {
      if (pid == 0) {
        Ring val = Ring(0);
        Ring valn = Ring(0);
        share.pushValue(Ring(0));
        tpShare.pushValues(Ring(0));
        for (int i = 1; i < nP; i++) {
          rgen.pi(i).random_data(&val, sizeof(Ring));
          tpShare.pushValues(val);
          valn += val;
        }
        valn = secret - valn;
        tpShare.pushValues(valn);
        rand_sh_sec.push_back(valn);
      } else {
        if (pid != nP) {
          Ring val;
          rgen.p0().random_data(&val, sizeof(Ring));
          share.pushValue(val);
        } else {
          share.pushValue(rand_sh_sec[idx_rand_sh_sec]);
          idx_rand_sh_sec++;
        }
      }
    }




    void OfflineEvaluator::setWireMasksParty(const std::unordered_map<common::utils::wire_t, int>& input_pid_map, 
                                            std::vector<Ring>& rand_sh_sec) {
      size_t idx_rand_sh_sec = 0;
      size_t idx_delta_sh = 0;
      size_t b_idx_rand_sh_sec = 0;

      for (const auto& level : circ_.gates_by_level) {
        for (const auto& gate : level) {
          switch (gate->type) {
            case common::utils::GateType::kInp: {
              auto pregate = std::make_unique<PreprocInput<Ring>>();
              auto pid = input_pid_map.at(gate->out);
              pregate->pid = pid;
              preproc_.gates[gate->out] = std::move(pregate);
              break;
            }

            case common::utils::GateType::kRec: {
              auto pregate = std::make_unique<PreprocRecGate<Ring>>();
              // King party (party 1) receives the reconstructed value
              bool is_king = (id_ == 1);
              pregate->Pking = is_king;
              preproc_.gates[gate->out] = std::move(pregate);
              break;
            }

            case common::utils::GateType::kMul: {
              AddShare<Ring> triple_a; // Holds one beaver triple share of a random value a
              TPShare<Ring> tp_triple_a; // Holds all the beaver triple shares of a random value a
              AddShare<Ring> triple_b; // Holds one beaver triple share of a random value b
              TPShare<Ring> tp_triple_b; // Holds all the beaver triple shares of a random value b
              AddShare<Ring> triple_c; // Holds one beaver triple share of c=a*b
              TPShare<Ring> tp_triple_c; // Holds all the beaver triple shares of c=a*b
              randomShare(nP_, id_, rgen_, triple_a, tp_triple_a);
              randomShare(nP_, id_, rgen_, triple_b, tp_triple_b);
              Ring tp_prod;
              if (id_ == 0) { tp_prod = tp_triple_a.secret() * tp_triple_b.secret(); }
              randomShareSecret(nP_, id_, rgen_, triple_c, tp_triple_c, tp_prod, rand_sh_sec, idx_rand_sh_sec);
              preproc_.gates[gate->out] =
                  std::move(std::make_unique<PreprocMultGate<Ring>>(triple_a, tp_triple_a, triple_b, tp_triple_b, triple_c, tp_triple_c));
              break;
            }

            case common::utils::GateType::kEqz: {
              AddShare<Ring> share_r1;
              TPShare<Ring> tp_share_r1;
              AddShare<Ring> share_r2;
              TPShare<Ring> tp_share_r2;
              std::vector<AddShare<Ring>> share_r1_bits(RINGSIZEBITS);
              std::vector<TPShare<Ring>> tp_share_r1_bits(RINGSIZEBITS);
              std::vector<AddShare<Ring>> share_r2_bits(RINGSIZEBITS);
              std::vector<TPShare<Ring>> tp_share_r2_bits(RINGSIZEBITS);
              Ring tp_r1 = Ring(0);
              Ring tp_r2 = Ring(0);
              std::vector<Ring> tp_r1_bits(RINGSIZEBITS);
              std::vector<Ring> tp_r2_bits(RINGSIZEBITS);

              // sharing r1 and r1_bits
              randomShare(nP_, id_, rgen_, share_r1, tp_share_r1);
              
              if (id_ == 0) {
                tp_r1 = tp_share_r1.secret();
                tp_r1_bits = bitDecomposeToInt(tp_r1);
              }
              for (int i = 0; i < RINGSIZEBITS; ++i) {
                  randomShareSecret(nP_, id_, rgen_, share_r1_bits[i], tp_share_r1_bits[i], tp_r1_bits[i],
                                                        rand_sh_sec, idx_rand_sh_sec);                                      
              }

              // sharing r2 and r2_bits
              if (id_ == 0) {
                rgen_.p0().random_data(&tp_r2, sizeof(Ring));
                tp_r2 = tp_r2 % RINGSIZEBITS; // make sure r2 is in [0, RINGSIZEBITS-1]
              }
              randomShareSecret(nP_, id_, rgen_, share_r2, tp_share_r2, tp_r2, rand_sh_sec, idx_rand_sh_sec);

              if (id_ == 0) {
                tp_r2 = tp_share_r2.secret();
                for (int i = 0; i < RINGSIZEBITS; ++i) {
                  if (i == tp_r2 % RINGSIZEBITS) {
                    tp_r2_bits[i] = 1;
                  } else {
                    tp_r2_bits[i] = 0;
                  }
                }
              }

              for (int i = 0; i < RINGSIZEBITS; ++i) {
                randomShareSecret(nP_, id_, rgen_, share_r2_bits[i], tp_share_r2_bits[i], tp_r2_bits[i],
                                                        rand_sh_sec, idx_rand_sh_sec);
              }
              preproc_.gates[gate->out] =
                  std::make_unique<PreprocEqzGate<Ring>>(share_r1, tp_share_r1, share_r2, tp_share_r2, share_r1_bits, tp_share_r1_bits, share_r2_bits, tp_share_r2_bits);
              break;
            }

            default: {
              break;
            }
          }
        }
      }
    }


    void OfflineEvaluator::setWireMasks(const std::unordered_map<common::utils::wire_t, int>& input_pid_map) {
      std::vector<Ring> rand_sh_sec;

      if (id_ == 0) {
        setWireMasksParty(input_pid_map, rand_sh_sec);

        size_t rand_sh_sec_num = rand_sh_sec.size();
        size_t arith_comm = rand_sh_sec_num;
        std::vector<size_t> lengths(2);
        lengths[0] = arith_comm;
        lengths[1] = rand_sh_sec_num;

        network_->send(nP_, lengths.data(), sizeof(size_t) * lengths.size());

        std::vector<Ring> offline_arith_comm(arith_comm);

        for (size_t i = 0; i < rand_sh_sec_num; i++) {
          offline_arith_comm[i] = rand_sh_sec[i];
        }
        network_->send(nP_, offline_arith_comm.data(), sizeof(Ring) * arith_comm);

      } else if (id_ != nP_) {

        setWireMasksParty(input_pid_map, rand_sh_sec);

      } else {

        std::vector<size_t> lengths(2);
        usleep(latency_);
        network_->recv(0, lengths.data(), sizeof(size_t) * lengths.size());
        size_t arith_comm = lengths[0];
        size_t rand_sh_sec_num = lengths[1];

        std::vector<Ring> offline_arith_comm(arith_comm);
        network_->recv(0, offline_arith_comm.data(), sizeof(Ring) * arith_comm);

        rand_sh_sec.resize(rand_sh_sec_num);
        for (int i = 0; i < rand_sh_sec_num; i++) {
          rand_sh_sec[i] = offline_arith_comm[i];
        }

        setWireMasksParty(input_pid_map, rand_sh_sec);
      }
    }

    PreprocCircuit<Ring> OfflineEvaluator::getPreproc() {
      return std::move(preproc_);
    }

    PreprocCircuit<Ring> OfflineEvaluator::run(const std::unordered_map<common::utils::wire_t, int>& input_pid_map) {
      setWireMasks(input_pid_map);
      return std::move(preproc_);
    }

};  // namespace shlut
