#pragma once

#include "../utils/circuit.h"
#include "sharing.h"
#include "../utils/types.h"
#include <unordered_map>

using namespace common::utils;

namespace shlut {
// Preprocessed data for a gate.
template <class R>
struct PreprocGate {
  PreprocGate() = default;
  virtual ~PreprocGate() = default;
};

template <class R>
using preprocg_ptr_t = std::unique_ptr<PreprocGate<R>>;

template <class R>
struct PreprocInput : public PreprocGate<R> {
  // ID of party providing input on wire.
  int pid{};
  PreprocInput() = default;
  PreprocInput(int pid) 
      : PreprocGate<R>(), pid(pid) {}
  PreprocInput(const PreprocInput<R>& pregate) 
      : PreprocGate<R>(), pid(pregate.pid) {}
};

template <class R>
struct PreprocRecGate : public PreprocGate<R> {
  bool Pking = false;
  PreprocRecGate() = default;
  PreprocRecGate(bool Pking)
    : PreprocGate<R>(), Pking(Pking) {}
};

template <class R>
struct PreprocMultGateAdd : public PreprocGate<R> {
  // Secret shared product of inputs masks.
  Share<R, 1> triple_a; // Holds one beaver triple share of a random value a
  Share<R, 1> triple_b; // Holds one beaver triple share of a random value b
  Share<R, 1> triple_c; // Holds one beaver triple share of c=a*b
  PreprocMultGateAdd() = default;
  PreprocMultGateAdd(const Share<R, 1>& triple_a,
                  const Share<R, 1>& triple_b,
                  const Share<R, 1>& triple_c)
      : PreprocGate<R>(), triple_a(triple_a),
        triple_b(triple_b),
        triple_c(triple_c) {}
};

template <class R>
struct PreprocEqzGate : public PreprocGate<R> {
  AddShare<R> share_r1;
  TPShare<R> tp_share_r1;
  AddShare<R> share_r2;
  TPShare<R> tp_share_r2;
  std::vector<AddShare<R>> share_r1_bits;
  std::vector<TPShare<R>> tp_share_r1_bits;
  std::vector<AddShare<R>> share_r2_bits;
  std::vector<TPShare<R>> tp_share_r2_bits;
  PreprocEqzGate() = default;
  PreprocEqzGate(const AddShare<R> &share_r1, const TPShare<R> &tp_share_r1,
                 const AddShare<R> &share_r2, const TPShare<R> &tp_share_r2,
                 const std::vector<AddShare<R>> &share_r1_bits, const std::vector<TPShare<R>> &tp_share_r1_bits,
                 const std::vector<AddShare<R>> &share_r2_bits, const std::vector<TPShare<R>> &tp_share_r2_bits)
    : PreprocGate<R>(), share_r1(share_r1), tp_share_r1(tp_share_r1), share_r2(share_r2), tp_share_r2(tp_share_r2), share_r1_bits(share_r1_bits), tp_share_r1_bits(tp_share_r1_bits), share_r2_bits(share_r2_bits), tp_share_r2_bits(tp_share_r2_bits) {}
};

// Preprocessed data for the circuit.
template <class R>
struct PreprocCircuit {
  std::unordered_map<wire_t, preprocg_ptr_t<R>> gates;
  PreprocCircuit() = default;
};
};  // namespace shlut
