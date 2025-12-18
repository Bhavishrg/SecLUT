#pragma once
#include <emp-tool/emp-tool.h>

#include <vector>
#include <algorithm>
#include "../utils/helpers.h"

using namespace common::utils;

namespace shlut {

// Collection of PRGs.
class RandGenPool {
  int id_;

  emp::PRG k_all;
  std::vector<emp::PRG> k_pi;
  

 public:
  explicit RandGenPool(int my_id, uint64_t seed = 200);
  
  emp::PRG& all(); // { return k_all; }
  emp::PRG& pi(int i);
};
};  // namespace shlut
