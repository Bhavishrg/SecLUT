#pragma once

#include <algorithm>
#include <array>
#include <boost/format.hpp>
#include <cmath>
#include <iostream>
#include <memory>
#include <optional>
#include <stdexcept>
#include <type_traits>
#include <unordered_map>
#include <vector>

#include "helpers.h"
#include "types.h"

namespace common::utils {

using wire_t = size_t;

enum GateType {
  kInp,
  kRec,
  kAdd,
  kMul,
  kSub,
  kEqz,
  kConstAdd,
  kConstMul,
  kInvalid,
  NumGates
};

std::ostream& operator<<(std::ostream& os, GateType type);

// Gates represent primitive operations.
// All gates have one output.
struct Gate {
  GateType type{GateType::kInvalid};
  int owner;
  wire_t out;
  std::vector<wire_t> outs;
  std::vector<std::vector<wire_t>> multi_outs;

  Gate() = default;
  Gate(GateType type, wire_t out);
  Gate(GateType type, int owner, wire_t out, std::vector<wire_t> outs);
  Gate(GateType type, int owner, wire_t out, std::vector<std::vector<wire_t>> multi_outs);

  virtual ~Gate() = default;
};

// Represents a gate with fan-in 2.
struct FIn2Gate : public Gate {
  wire_t in1{0};
  wire_t in2{0};

  FIn2Gate() = default;
  FIn2Gate(GateType type, wire_t in1, wire_t in2, wire_t out);
};

struct FIn3Gate : public Gate {
  wire_t in1{0};
  wire_t in2{0};
  wire_t in3{0};

  FIn3Gate() = default;
  FIn3Gate(GateType type, wire_t in1, wire_t in2, wire_t in3, wire_t out);
};

struct FIn4Gate : public Gate {
  wire_t in1{0};
  wire_t in2{0};
  wire_t in3{0};
  wire_t in4{0};

  FIn4Gate() = default;
  FIn4Gate(GateType type, wire_t in1, wire_t in2, wire_t in3, wire_t in4, wire_t out);
};

// Represents a gate with fan-in 1.
struct FIn1Gate : public Gate {
  wire_t in{0};

  FIn1Gate() = default;
  FIn1Gate(GateType type, wire_t in, wire_t out);
};

// Represents a gate used to denote SIMD operations.
// These type is used to represent operations that take vectors of inputs but
// might not necessarily be SIMD e.g., dot product.
struct SIMDGate : public Gate {
  std::vector<wire_t> in1{0};
  std::vector<wire_t> in2{0};

  SIMDGate() = default;
  SIMDGate(GateType type, std::vector<wire_t> in1, std::vector<wire_t> in2, wire_t out);
};

// Represents gates where one input is a constant.
template <class R>
struct ConstOpGate : public Gate {
  wire_t in{0};
  R cval;

  ConstOpGate() = default;
  ConstOpGate(GateType type, wire_t in, R cval, wire_t out)
      : Gate(type, out), in(in), cval(std::move(cval)) {}
};

using gate_ptr_t = std::shared_ptr<Gate>;

// Gates ordered by multiplicative depth.
//
// Addition gates are not considered to increase the depth.
// Moreover, if gates_by_level[l][i]'s output is input to gates_by_level[l][j]
// then i < j.
struct LevelOrderedCircuit {
  size_t num_gates;
  size_t num_wires;
  std::array<uint64_t, GateType::NumGates> count;
  std::vector<wire_t> outputs;
  std::vector<std::vector<gate_ptr_t>> gates_by_level;

  friend std::ostream& operator<<(std::ostream& os, const LevelOrderedCircuit& circ);
};

// Represents an arithmetic circuit.
template <class R>
class Circuit {
  std::vector<wire_t> outputs_;
  std::vector<gate_ptr_t> gates_;
  size_t num_wires;

  bool isWireValid(wire_t wid) { return wid < num_wires; }

 public:
  Circuit() : num_wires(0) {}

  // Methods to manually build a circuit.
  wire_t newInputWire() {
    wire_t wid = num_wires;
    gates_.push_back(std::make_shared<Gate>(GateType::kInp, wid));
    num_wires += 1;
    return wid;
  }

  void setAsOutput(wire_t wid) {
    if (!isWireValid(wid)) {
      throw std::invalid_argument("Invalid wire ID.");
    }

    outputs_.push_back(wid);
  }

  // Function to add a gate with fan-in 2.
  wire_t addGate(GateType type, wire_t input1, wire_t input2) {
    if (type != GateType::kAdd && type != GateType::kMul &&
        type != GateType::kSub) {
      throw std::invalid_argument("Invalid gate type.");
    }

    if (!isWireValid(input1) || !isWireValid(input2)) {
      throw std::invalid_argument("Invalid wire ID.");
    }

    wire_t output = num_wires;
    gates_.push_back(std::make_shared<FIn2Gate>(type, input1, input2, output));
    num_wires += 1;

    return output;
  }

  // Function to add a gate with one input from a wire and a second constant
  // input.
  wire_t addConstOpGate(GateType type, wire_t wid, R cval) {
    if (type != kConstAdd && type != kConstMul) {
      throw std::invalid_argument("Invalid gate type.");
    }

    if (!isWireValid(wid)) {
      throw std::invalid_argument("Invalid wire ID.");
    }

    wire_t output = num_wires;
    gates_.push_back(std::make_shared<ConstOpGate<R>>(type, wid, cval, output));
    num_wires += 1;

    return output;
  }

  // Function to add a single input gate.
  wire_t addGate(GateType type, wire_t input) {
    if (type != GateType::kEqz && type != GateType::kRec) {
      throw std::invalid_argument("Invalid gate type.");
    }

    if (!isWireValid(input)) {
      throw std::invalid_argument("Invalid wire ID.");
    }

    wire_t output = num_wires;
    gates_.push_back(std::make_shared<FIn1Gate>(type, input, output));
    num_wires += 1;

    return output;
  }





  
  // Level ordered gates are helpful for evaluation.
  [[nodiscard]] LevelOrderedCircuit orderGatesByLevel() const {
    LevelOrderedCircuit res;
    res.outputs = outputs_;
    res.num_gates = gates_.size();
    res.num_wires = num_wires;

    // Map from output wire id to multiplicative depth/level.
    // Input gates have a depth of 0.
    std::vector<size_t> gate_level(num_wires, 0);
    size_t depth = 0;

    // This assumes that if gates_[i]'s output is input to gates_[j] then
    // i < j.
    for (const auto& gate : gates_) {
      switch (gate->type) {
        case GateType::kRec: {
          const auto* g = static_cast<FIn1Gate*>(gate.get());
          gate_level[g->out] = gate_level[g->in] + 1;
          depth = std::max(depth, gate_level[gate->out]);
          break;
        }
        case GateType::kAdd:
        case GateType::kSub: {
          const auto* g = static_cast<FIn2Gate*>(gate.get());
          gate_level[g->out] = std::max(gate_level[g->in1], gate_level[g->in2]);
          depth = std::max(depth, gate_level[gate->out]);
          break;
        }

        case GateType::kMul: {
          const auto* g = static_cast<FIn2Gate*>(gate.get());
          gate_level[g->out] = std::max(gate_level[g->in1], gate_level[g->in2]) + 1;
          depth = std::max(depth, gate_level[gate->out]);
          break;
        }

        case GateType::kConstAdd:
        case GateType::kConstMul: {
          const auto* g = static_cast<ConstOpGate<R>*>(gate.get());
          gate_level[g->out] = gate_level[g->in];
          depth = std::max(depth, gate_level[gate->out]);
          break;
        }

        case GateType::kEqz: {
          const auto* g = static_cast<FIn1Gate*>(gate.get());
          gate_level[g->out] = gate_level[g->in] + 1;
          depth = std::max(depth, gate_level[gate->out]);
          break;
        }

        case GateType::kInp:
          // Input gates have depth 0.
          break;

        default:
          break;
      }
    }

    std::fill(res.count.begin(), res.count.end(), 0);

    std::vector<std::vector<gate_ptr_t>> gates_by_level(depth + 1);
    for (const auto& gate : gates_) {
      res.count[gate->type]++;
      // if (gate->type == GateType::kShuffle || gate->type == GateType::kPermAndSh || gate->type == GateType::kPublicPerm || gate->type == GateType::kCompact || gate->type == GateType::kGroupwiseIndex || gate->type == GateType::kGroupwisePropagate || gate->type == GateType::kSort || gate->type == GateType::kRewire) {
      //   gates_by_level[gate_level[gate->outs[0]]].push_back(gate);
      // } else {
        gates_by_level[gate_level[gate->out]].push_back(gate);
      // } 
    }

    res.gates_by_level = std::move(gates_by_level);

    return res;
  }

  // Evaluate circuit on plaintext inputs.
  [[nodiscard]] std::vector<R> evaluate(const std::unordered_map<wire_t, R>& inputs) const {
    auto level_circ = orderGatesByLevel();
    std::vector<R> wires(level_circ.num_gates);

    auto num_inp_gates = level_circ.count[GateType::kInp];
    if (inputs.size() != num_inp_gates) {
      throw std::invalid_argument(boost::str(
          boost::format("Expected %1% inputs but received %2% inputs.") %
          num_inp_gates % inputs.size()));
    }

    for (const auto& level : level_circ.gates_by_level) {
      for (const auto& gate : level) {
        switch (gate->type) {
          case GateType::kInp: {
            wires[gate->out] = inputs.at(gate->out);
            break;
          }

          case GateType::kMul: {
            auto* g = static_cast<FIn2Gate*>(gate.get());
            wires[g->out] = wires[g->in1] * wires[g->in2];
            break;
          }

          case GateType::kAdd: {
            auto* g = static_cast<FIn2Gate*>(gate.get());
            wires[g->out] = wires[g->in1] + wires[g->in2];
            break;
          }

          case GateType::kSub: {
            auto* g = static_cast<FIn2Gate*>(gate.get());
            wires[g->out] = wires[g->in1] - wires[g->in2];
            break;
          }

          case GateType::kConstAdd: {
            auto* g = static_cast<ConstOpGate<R>*>(gate.get());
            wires[g->out] = wires[g->in] + g->cval;
            break;
          }

          case GateType::kConstMul: {
            auto* g = static_cast<ConstOpGate<R>*>(gate.get());
            wires[g->out] = wires[g->in] * g->cval;
            break;
          }

          case GateType::kEqz: {
            auto* g = static_cast<FIn1Gate*>(gate.get());
            if (wires[g->in] == 0) {
              wires[g->out] = 1;
            }
            else {
              wires[g->out] = 0;
            }
            break;
          }

          default: {
            throw std::runtime_error("Invalid gate type.");
          }
        }
      }
    }

    std::vector<R> outputs;
    for (auto i : level_circ.outputs) {
      outputs.push_back(wires[i]);
    }

    return outputs;
  }

};
};  // namespace common::utils
