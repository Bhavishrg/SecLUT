#include "rand_gen_pool.h"

#include <algorithm>

#include "../utils/helpers.h"

namespace shlut {

  RandGenPool::RandGenPool(int my_id, uint64_t seed) 
    : id_{my_id}, k_pi(3) { 
    
    
    auto seed_block = emp::makeBlock(seed, 0); 
    k_all.reseed(&seed_block, 0);

    for (int i = 0; i < 3; i++) {
      k_pi[i].reseed(&seed_block, 0); 
    }
    
  }


emp::PRG& RandGenPool::all() { return k_all; }

emp::PRG& RandGenPool::pi(int i) { return k_pi[i]; }

}  // namespace shlut
