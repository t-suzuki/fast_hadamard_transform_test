#include "fast_hadamard_transform.hpp"
#include <vector>
#include <iostream>

int main(int argc, char* argv[]) {
  const size_t N = 8;

  auto print_buf = [N](std::string name, float* p) {
    std::cout << name << ": ";
    for (size_t i=0; i<N; ++i) {
      std::cout << p[i] << " ";
    }
      std::cout << std::endl;
  };
  float buf[N] = {1, 0, 1, 0, 0, 1, 1, 0};

  using namespace fast_hadamard_transform;

  { // non-unitary
    std::vector<float> buf_ht(N);
    fht(buf_ht.data(), buf, N);

    std::vector<float> buf_ht_ht(N);
    fht(buf_ht_ht.data(), buf_ht.data(), N);

    std::vector<float> buf_ht_iht(N);
    ifht(buf_ht_iht.data(), buf_ht.data(), N);

    std::cout << "Non-unitary.." << std::endl;
    print_buf("buf", buf);
    print_buf("HT(buf)", buf_ht.data());
    print_buf("HT(HT(buf))", buf_ht_ht.data());
    print_buf("IHT(HT(buf))", buf_ht_iht.data());
  }

  { // unitary
    std::vector<float> bufu_ht(N);
    fht(bufu_ht.data(), buf, N, true);

    std::vector<float> bufu_ht_ht(N);
    fht(bufu_ht_ht.data(), bufu_ht.data(), N, true);

    std::vector<float> bufu_ht_iht(N);
    ifht(bufu_ht_iht.data(), bufu_ht.data(), N, true);

    std::cout << "Unitary.." << std::endl;
    print_buf("buf", buf);
    print_buf("HT(buf)", bufu_ht.data());
    print_buf("HT(HT(buf))", bufu_ht_ht.data());
    print_buf("IHT(HT(buf))", bufu_ht_iht.data());
  }


  return 0;
}
