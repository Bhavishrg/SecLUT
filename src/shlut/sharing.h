#pragma once

#include <emp-tool/emp-tool.h>

#include <array>
#include <vector>

#include "../utils/helpers.h"
#include "../utils/types.h"

using namespace common::utils;

namespace shlut {

// Maybe I want to remove size_t N here
template <class R, size_t N>
class Share {
  std::array<R, N> value;

  public:
    Share() = default;
    explicit Share(std::vector value)
        : value_(value) {}
    
    void value(std::array<R, N>& array) {array = value;}

    //Finish the operators which can be done.
    Share<R, N>& operator+=(const Share<R, N>& rhs) {
      for (size_t i = 0; i < value.size(); i++) {
        value[i] += rhs.value[i];
      }
      return *this;
    }

    friend Share<R, N> operator+(Share<R, N> lhs, const Share<R, N>& rhs) {
      lhs += rhs;
      return lhs;
    }

    Share<R, N>& operator*=(const R& rhs) {
      for (size_t i = 0; i < value.size(); i++) {
        value[i] *= rhs;
      }
      return *this;
    }

    friend Share<R, N> operator*(Share<R, N> lhs, const R& rhs) {
      lhs *= rhs;
      return lhs;
    }

    Share<R, N>& operator-=(const Share<R, N>& rhs) {
      (*this) += (rhs * R(-1));
      return *this;
    }

    friend Share<R, N> operator-(Share<R, N> lhs, const Share<R, N>& rhs) {
      lhs -= rhs;
      return lhs;
    }
};

template <class R, size_t N>
class MultiShare {
  std::vector<std::array<R, N>> values_;

  public:
  MultiShare() = default;
  explicit MultiShare(std::vector<std::array<R, N>> value)
      : values_{std::move(value)} {}

  // Access share elements.
  std::array<R, N>& operator[](size_t idx) { return values_.at(idx); }
  std::array<R, N> operator[](size_t idx) const { return values_.at(idx); }

  void pushValues(std::array<R, N> val) { values_.push_back(val); }

  // Arithmetic operators.
  MultiShare<R, N>& operator+=(const MultiShare<R, N>& rhs) {
    for (size_t i = 0; i < values_.size(); i++) {
      for (size_t j = 0; j < N; j++) {
        values_[i][j] += rhs.values_[i][j];
      }
    }
    return *this;
  }

  friend MultiShare<R, N> operator+(MultiShare<R, N> lhs,
                                      const MultiShare<R, N>& rhs) {
    lhs += rhs;
    return lhs;
  }

  MultiShare<R, N>& operator*=(const R& rhs) {
    for (size_t i = 0; i < values_.size(); i++) {
      for (size_t j = 0; j < N; j++) {
        values_[i][j] *= rhs;
      }
    }
    return *this;
  }

  friend MultiShare<R, N> operator*(MultiShare<R, N> lhs, const R& rhs) {
    lhs *= rhs;
    return lhs;
  }

  MultiShare<R, N>& operator-=(const MultiShare<R, N>& rhs) {
    (*this) += (rhs * R(-1));
    return *this;
  }

  friend MultiShare<R, N> operator-(MultiShare<R, N> lhs,
                                      const MultiShare<R, N>& rhs) {
    lhs -= rhs;
    return lhs;
  }
};

// Make a class for party-local additive shares.
template <class R>
class PartyShare {
  int pid_;
  R value1_;
  R value2_;

public:
  PartyShare() = default;
  explicit PartyShare(int pid, R value1, R value2)
      : pid_(pid), value1_(value1), value2_(value2) {}

  int pid() const { return pid_; }
  R value1() const { return value1_; }
  R value2() const { return value2_; }

  //Arithmetic operators.
  PartyShare<R>& operator+=(const PartyShare<R>& rhs) {
    value1_ += rhs.value1_;
    value2_ += rhs.value2_;
    return *this;
  }

  friend PartyShare<R> operator+(PartyShare<R> lhs,
                                      const PartyShare<R>& rhs) {
    lhs += rhs;
    return lhs;
  }

  PartyShare<R>& operator+=(const R& rhs) {
    if (pid_ == 1) {
      value2_ += rhs;
    } else if (pid_ == 2) {
      value1_ += rhs;
    }

    return *this;
  }

  friend PartyShare<R> operator+(PartyShare<R> lhs, const R& rhs) {
    lhs += rhs;
    return lhs;
  }
  
  PartyShare<R>& operator*=(const R& rhs) {
    value1_ *= rhs;
    value2_ *= rhs;
    return *this;
  }
  
  friend PartyShare<R> operator*(PartyShare<R> lhs, const R& rhs) {
    lhs *= rhs;
    return lhs;
  }

  friend PartyShare<R> operator*(const R& lhs, PartyShare<R> rhs) {
    rhs *= lhs;
    return rhs;
  }
  
  PartyShare<R>& operator-=(const PartyShare<R>& rhs) {
    (*this) += (rhs * R(-1));
    return *this;
  }

  friend PartyShare<R> operator-(PartyShare<R> lhs,
                                      const PartyShare<R>& rhs) {
    lhs -= rhs;
    return lhs;
  }

  // Need to finish few more functions later.
};

template <class R>
class AddShare {
  R value_;
  
 public:
  AddShare() = default;
  explicit AddShare(R value)
      : value_{value} {}

  void randomize(emp::PRG& prg) {
    randomizeZZp(prg, value_.data(), sizeof(R));
  }

  R& valueAt() { return value_; }

  void pushValue(R val) { value_ = val; }
  
  R valueAt() const { return value_; }

  // Arithmetic operators.
  AddShare<R>& operator+=(const AddShare<R>& rhs) {
    value_ += rhs.value_;
    return *this;
  }

  // what is "friend"?
  friend AddShare<R> operator+(AddShare<R> lhs,
                                      const AddShare<R>& rhs) {
    lhs += rhs;
    return lhs;
  }

  AddShare<R>& operator-=(const AddShare<R>& rhs) {
    (*this) += (rhs * R(-1));
    return *this;
  }

  friend AddShare<R> operator-(AddShare<R> lhs,
                                      const AddShare<R>& rhs) {
    lhs -= rhs;
    return lhs;
  }

  AddShare<R>& operator*=(const R& rhs) {
    value_ *= rhs;
    return *this;
  }

  friend AddShare<R> operator*(AddShare<R> lhs, const R& rhs) {
    lhs *= rhs;
    return lhs;
  }

  AddShare<R>& operator<<=(const int& rhs) {
    uint64_t value = conv<uint64_t>(value_);
    value <<= rhs;
    value_ = value;
    return *this;
  }

  friend AddShare<R> operator<<(AddShare<R> lhs, const int& rhs) {
    lhs <<= rhs;
    return lhs;
  }

  AddShare<R>& operator>>=(const int& rhs) {
    uint64_t value = conv<uint64_t>(value_);
    value >>= rhs;
    value_ = value;
    return *this;
  }

  friend AddShare<R> operator>>(AddShare<R> lhs, const int& rhs) {
    lhs >>= rhs;
    return lhs;
  }

  AddShare<R>& add(R val, int pid) {
    if (pid == 1) {
      value_ += val;
    }
    return *this;
  }

  AddShare<R>& addWithAdder(R val, int pid, int adder) {
    if (pid == adder) {
      value_ += val;
    }
    return *this;
  }

  AddShare<R>& shift() {
    auto bits = bitDecomposeTwo(value_);
    if (bits[63] == 1)
      value_ = 1;
    else
      value_ = 0;
    return *this;
  }
  
};

template <class R>
class TPShare {
  std::vector<R> values_;

  public:
  TPShare() = default;
  explicit TPShare(std::vector<R> value)
      : values_{std::move(value)} {}

  // Access share elements.
  // idx = i retreives value common with party having i.
  R& operator[](size_t idx) { return values_.at(idx); }
  // idx = i retreives tag common with party having i.
  //R& operator()(size_t idx) { return tags_.at(idx); }
  
  R operator[](size_t idx) const { return values_.at(idx); }
  //R operator()(size_t idx) { return tags_.at(idx); }

  R& commonValueWithParty(int pid) {
    return values_.at(pid);
  }

  [[nodiscard]] R commonValueWithParty(int pid) const {
    return values_.at(pid);
  }

  void pushValues(R val) { values_.push_back(val); }

  [[nodiscard]] R secret() const { 
    R res = values_[0];
    for (int i = 1; i < values_.size(); i++)
     res += values_[i];
    return res;
  }
  // Arithmetic operators.
  TPShare<R>& operator+=(const TPShare<R>& rhs) {
    for (size_t i = 1; i < values_.size(); i++) {
      values_[i] += rhs.values_[i];
    }
    return *this;
  }

  friend TPShare<R> operator+(TPShare<R> lhs, const TPShare<R>& rhs) {
    lhs += rhs;
    return lhs;
  }

  TPShare<R>& operator-=(const TPShare<R>& rhs) {
    (*this) += (rhs * R(-1));
    return *this;
  }

  friend TPShare<R> operator-(TPShare<R> lhs, const TPShare<R>& rhs) {
    lhs -= rhs;
    return lhs;
  }

  TPShare<R>& operator*=(const R& rhs) {
    for (size_t i = 1; i < values_.size(); i++) {
      values_[i] *= rhs;
    }
    return *this;
  }

  friend TPShare<R> operator*(TPShare<R> lhs, const R& rhs) {
    lhs *= rhs;
    return lhs;
  }

  TPShare<R>& operator<<=(const int& rhs) {
    for (size_t i = 1; i < values_.size(); i++) {
        uint64_t value = conv<uint64_t>(values_[i]);
        value <<= rhs;
        values_[i] = value;
    }
    return *this;
  }

  friend TPShare<R> operator<<(TPShare<R> lhs, const int& rhs) {
    lhs <<= rhs;
    return lhs;
  }

  TPShare<R>& operator>>=(const int& rhs) {
    for (size_t i = 1; i < values_.size(); i++) {
        uint64_t value = conv<uint64_t>(values_[i]);
        value >>= rhs;
        values_[i] = value;
    }
    return *this;
  }

  friend TPShare<R> operator>>(TPShare<R> lhs, const int& rhs) {
    lhs >>= rhs;
    return lhs;
  }

  AddShare<R> getAS(size_t pid) {
    return AddShare<R>({values_.at(pid)});
  }

  TPShare<R>& shift() {
    for (size_t i = 1; i < values_.size(); i++) {
      auto bits = bitDecomposeTwo(values_[i]);
      if (bits[63] == 1)
        values_[i] = 1;
      else 
        values_[i] = 0;
    }
    return *this;
  }

  //Add above
  
};

template <>
void AddShare<BoolRing>::randomize(emp::PRG& prg);
//add the constructor above



// Contains all elements of a secret sharing. Used only for generating dummy
// preprocessing data.
/*
template <class R>
struct DummyShare { 
  // number of components will depent upon number of parties
  std::array<R, 6> share_elements;

  DummyShare() = default;

  explicit DummyShare(std::array<R, 6> share_elements)
      : share_elements(std::move(share_elements)) {}

  DummyShare(R secret, emp::PRG& prg) {
    prg.random_data(share_elements.data(), sizeof(R) * 5);

    R sum = share_elements[0];
    for (int i = 1; i < 5; ++i) {
      sum += share_elements[i];
    }
    share_elements[5] = secret - sum;
  }

  void randomize(emp::PRG& prg) {
    prg.random_data(share_elements.data(), sizeof(R) * 6);
  }

  [[nodiscard]] R secret() const {
    R sum = share_elements[0];
    for (size_t i = 1; i < 6; ++i) {
      sum += share_elements[i];
    }

    return sum;
  }

  DummyShare<R>& operator+=(const DummyShare<R>& rhs) {
    for (size_t i = 0; i < 6; ++i) {
      share_elements[i] += rhs.share_elements[i];
    }

    return *this;
  }

  friend DummyShare<R> operator+(DummyShare<R> lhs, const DummyShare<R>& rhs) {
    lhs += rhs;
    return lhs;
  }

  DummyShare<R>& operator-=(const DummyShare<R>& rhs) {
    for (size_t i = 0; i < 6; ++i) {
      share_elements[i] -= rhs.share_elements[i];
    }

    return *this;
  }

  friend DummyShare<R> operator-(DummyShare<R> lhs, const DummyShare<R>& rhs) {
    lhs -= rhs;
    return lhs;
  }

  DummyShare<R>& operator*=(const R& rhs) {
    for (size_t i = 0; i < 6; ++i) {
      share_elements[i] *= rhs;
    }

    return *this;
  }

  friend DummyShare<R> operator*(DummyShare<R> lhs, const R& rhs) {
    lhs *= rhs;
    return lhs;
  }

  friend DummyShare<R> operator*(const R& lhs, DummyShare<R> rhs) {
    // Assumes abelian ring.
    rhs *= lhs;
    return rhs;
  }

  //Dealer<R> getRSS(size_t pid) {
  //  return Dealer<R>({getShareElement(pid, pidFromOffset(pid, 1)),
  //                             getShareElement(pid, pidFromOffset(pid, 2)),
  //                             getShareElement(pid, pidFromOffset(pid, 3))});
  //}

  R getShareElement(size_t i, size_t j) {
    return share_elements.at(upperTriangularToArray(i, j));
  }
};*/

//template <>
//void AddShare<BoolRing>::randomize(emp::PRG& prg);

//template <>
//TPShare<BoolRing>::TPShare(BoolRing secret, emp::PRG& prg);

//template <>
//void TPShare<BoolRing>::randomize(emp::PRG& prg);
};  // namespace shlut
