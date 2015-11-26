// Fast Hadamard Transform
#ifndef __FHT_HPP__
#define __FHT_HPP__
#include <cstddef>
#include <cmath>
#include <vector>

namespace fast_hadamard_transform {

namespace detail {
size_t power_of_two(size_t sz) {
  size_t n = 0;
  size_t x = 1;
  while (x < sz) {
    ++n;
    x *= 2;
    if (x > x*2) break; // overflow
  }
  return n;
}
}

template <typename T>
bool fht(T* __restrict dst, T* __restrict src, size_t sz, bool unitary=false) {
  if (!dst) return false;
  if (!src) return false;
  const size_t n = detail::power_of_two(sz);
  if (std::pow(2, n) != sz) return false;
  std::vector<T> tmp0(sz), tmp1(sz);
  std::copy(src, src + sz, tmp0.begin());
  for (size_t step=sz/2; step>0; step/=2) {
    const size_t skip = step*2;
    for (size_t i=0; i<sz; i+=skip) {
      for (size_t j=0; j<step; ++j) {
        tmp1[i        + j] = tmp0[i + j] + tmp0[i + step + j];
        tmp1[i + step + j] = tmp0[i + j] - tmp0[i + step + j];
      }
    }
    std::swap(tmp0, tmp1);
  }
  if (unitary) {
    const T scale = 1.0/std::pow(std::sqrt(2.0), n);
    for (size_t i=0; i<sz; ++i) {
      dst[i] = tmp0[i]*scale;
    }
  } else {
    std::copy(tmp0.begin(), tmp0.end(), dst);
  }
  return true;
}

template <typename T>
bool ifht(T* __restrict dst, T* __restrict src, size_t sz, bool unitary=false) {
  if (unitary) {
    return fht(dst, src, sz, true);
  } else {
    if (!fht(dst, src, sz, false)) return false;
    for (size_t i=0; i<sz; ++i) {
      dst[i] = dst[i]/T(sz);
    }
    return true;
  }
}

}

#endif // __FHT_HPP__

