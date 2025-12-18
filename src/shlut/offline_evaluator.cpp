#include "offline_evaluator.h"

#include <NTL/BasicThreadPool.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <thread>

// #include "../utils/helpers.h"

namespace shlut {
OfflineEvaluator::OfflineEvaluator(int my_id,
                                   std::shared_ptr<io::NetIOMP> network, PreprocCircuit<Ring> preproc,
                                   common::utils::LevelOrderedCircuit circ,
                                   int threads, int seed, int latency, bool use_pking)
    : id_(my_id),
      latency_(latency),
      use_pking_(use_pking),
      rgen_(my_id, seed), 
      network_(std::move(network)),
      circ_(std::move(circ)),
      preproc_(std::move(preproc))

      {
        tpool_ = std::make_shared<common::utils::ThreadPool>(threads);
      }

OfflineEvaluator::OfflineEvaluator(int my_id,
                                   std::shared_ptr<io::NetIOMP> network, PreprocCircuit<Ring> preproc,
                                   common::utils::LevelOrderedCircuit circ, std::shared_ptr<common::utils::ThreadPool> tpool,
                                   int threads, int seed, int latency, bool use_pking)
    : id_(my_id),
      latency_(latency),
      use_pking_(use_pking),
      rgen_(my_id, seed), 
      network_(std::move(network)),
      circ_(std::move(circ)),
      preproc_(std::move(preproc)),
      tpool_(std::move(tpool)) {}

    void OfflineEvaluator::randomAdditiveShare(RandGenPool& rgen, Share<Ring, 1>& share) {
      Ring val = Ring(0);
      rgen.pi(id_).random_data(&val, sizeof(Ring));
      share[0] = val;
    }

    void OfflineEvaluator::randomReplicatedShare(RandGenPool& rgen, Share<Ring, 2>& share){
      Ring val1 = Ring(0);
      Ring val2 = Ring(0);
      rgen.pi((id_ + 1) % 3).random_data(&val1, sizeof(Ring));
      rgen.pi((id_ + 2) % 3).random_data(&val2, sizeof(Ring));
      share[0] = val1;
      share[1] = val2;
    }

    void OfflineEvaluator::randomAugmentedShare(RandGenPool& rgen, Share<Ring, 3>& share){
      Ring val1 = Ring(0);
      Ring val2 = Ring(0);
      Ring val3 = Ring(0);
      rgen.all().random_data(&val1, sizeof(Ring));
      rgen.pi((id_ + 1) % 3).random_data(&val2, sizeof(Ring));
      rgen.pi((id_ + 2) % 3).random_data(&val3, sizeof(Ring));
      share[0] = val1;
      share[1] = val2;
      share[2] = val3;
    }

    void OfflineEvaluator::AdditiveShareSecret(RandGenPool& rgen, Share<Ring, 1>& share, Ring secret, int pid) {
      if (id_ == pid) {
        Ring val = Ring(0);
        share[0] = Ring(0);
        for (int i = 0; i < nP_; i++) {
          if (i != pid) {
            rgen.pi(i).random_data(&val, sizeof(Ring));
            share[0] += val;
          }
        }
        share[0] += secret - share[0];
      } else {
        Ring val;
        rgen.pi(pid).random_data(&val, sizeof(Ring));
        share[0] = val;
      }
    }

    void OfflineEvaluator::ReplicatedShareSecret(RandGenPool& rgen, Share<Ring, 2>& share, Ring secret, int pid) {
      if (id_ == pid) {
        Ring val1 = Ring(0);
        Ring val2 = Ring(0);
        // Generate a random mask by secret + val
        rgen.pi((id_ + 1) % 3).random_data(&val1, sizeof(Ring));
        rgen.pi((id_ + 2) % 3).random_data(&val2, sizeof(Ring));
        Ring val = secret - val1 - val2;
        // Distribute val to other parties
        network_->send((id_ + 1) % 3, &val, sizeof(Ring));
        network_->send((id_ + 2) % 3, &val, sizeof(Ring));

        share[1] = val1;
        share[0] = val2;
      } else {
        Ring val_recv;
        network_->recv(pid, &val_recv, sizeof(Ring));
        Ring val = Ring(0);
        rgen.pi(pid).random_data(&val, sizeof(Ring));

        if ((id_ + 1) % 3 == pid) {
          share[1] = val;
          share[0] = val_recv;
        } else {
          share[0] = val;
          share[1] = val_recv;
        }
      }
    }

      void OfflineEvaluator::AugmentedShareSecret(RandGenPool& rgen, Share<Ring, 3>& share, Ring secret, int pid) {
        if(id_ == pid){
          Ring val = Ring(0);
          Ring val1 = Ring(0);
          Ring val2 = Ring(0);
          Ring val_send = Ring(0);
          
          rgen.self().random_data(&val, sizeof(Ring)); // Will be a secret share
          val += secret;
          rgen.pi((id_ + 1) % 3).random_data(&val1, sizeof(Ring));
          rgen.pi((id_ + 2) % 3).random_data(&val2, sizeof(Ring));

          val_send = secret - val1 - val2;
          network_->send((id_ + 1) % 3, &val_send, sizeof(Ring));
          network_->send((id_ + 2) % 3, &val_send, sizeof(Ring));
          network_->send((id_ + 1) % 3, &val, sizeof(Ring));
          network_->send((id_ + 2) % 3, &val, sizeof(Ring));

          share[2] = val1;
          share[1] = val2;
          share[0] = val;
        } else {
          Ring val_recv;
          network_->recv(pid, &val_recv, sizeof(Ring));
          Ring val = Ring(0);
          rgen.pi(pid).random_data(&val, sizeof(Ring));

          if ((id_ + 1) % 3 == pid) {
            share[2] = val;
            share[1] = val_recv;
          } else if ((id_ + 2) % 3 == pid) {
            share[1] = val;
            share[2] = val_recv;
          }

          Ring val;
          network_->recv(pid, &val, sizeof(Ring));
          share[0] = val;
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
            case common::utils::GateType::kInp:
            case common::utils::GateType::kInpRSS:
            case common::utils::GateType::kInpAug: {
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
            } // Check if I need this for RSS and Aug Rec as well 

            case common::utils::GateType::kMul: {
              Share<Ring, 1> triple_a; // Holds one beaver triple share of a random value a
              Share<Ring, 1> triple_b; // Holds one beaver triple share of a random value b
              Share<Ring, 1> triple_c; // Holds one beaver triple share of c=a*b
              randomAdditiveShare(rgen_, triple_a);
              randomAdditiveShare(rgen_, triple_b);
              Ring tp_prod;
              if (id_ == 0) { tp_prod = triple_a.secret() * triple_b.secret(); } // This needs to be changed as there is no dealer here
              AdditiveShareSecret(rgen_, triple_c, tp_prod, id_);
              preproc_.gates[gate->out] =
                  std::move(std::make_unique<PreprocMultGate<Ring>>(triple_a, triple_b, triple_c));
              break;
            } //Ask Bhavish how to make Beaver triples here



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

/*
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
    } // See what is the use of this function, mostly I can remove it
*/

    PreprocCircuit<Ring> OfflineEvaluator::getPreproc() {
      return std::move(preproc_);
    }

    PreprocCircuit<Ring> OfflineEvaluator::run(const std::unordered_map<common::utils::wire_t, int>& input_pid_map) {
      std::vector<Ring> rand_sh_sec;
      setWireMasksParty(input_pid_map, rand_sh_sec);
      return std::move(preproc_);
    }

};  // namespace shlut


/*
  I need to add preproc to my code for preprocessing. I also need to think about the gates which I want to use.
*/