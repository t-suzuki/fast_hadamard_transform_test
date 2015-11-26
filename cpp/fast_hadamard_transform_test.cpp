#include "fast_hadamard_transform.hpp"
#include <vector>
#include <iostream>
#include <chrono>

template <typename Func>
void timer(std::string name, size_t n_iter, size_t n_burnin, Func&& func) {
  for (size_t i=0; i<n_burnin; ++i) func();
  auto t = std::chrono::system_clock::now();
  for (size_t i=0; i<n_iter; ++i) {
    func();
  }
  auto diff = std::chrono::system_clock::now() - t;
  auto elapsed_us = std::chrono::duration_cast<std::chrono::microseconds>(diff).count();
  std::cout << "[" << name << "] " << elapsed_us << " us. " << (elapsed_us/double(n_iter)) << " us/iter." << std::endl;
}

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

  { // benchmark
    const size_t N_benchmark = 1024;
    const size_t n_iter = 100;
    const size_t n_burnin = 10;

    std::vector<float> buf_src(N_benchmark);
    std::vector<float> buf_dst(N_benchmark);
    for (size_t i=0; i<N_benchmark; ++i) buf_src[i] = i*10.5f/(i + 100.0f);

    timer("ad-hoc HT", n_iter, n_burnin, [&]{ 
        fht(buf_dst.data(), buf_src.data(), N_benchmark);
        ifht(buf_src.data(), buf_dst.data(), N_benchmark);
    });

    fht_plan<float> plan(N_benchmark);
    timer("planned iHT", n_iter, n_burnin, [&]{ 
        plan.fht(buf_dst.data(), buf_src.data());
        plan.ifht(buf_src.data(), buf_dst.data());
    });
  }

  return 0;
}
