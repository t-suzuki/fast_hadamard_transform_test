CXXFLAGS=-std=c++11 -Wall -Ofast -march=native

fast_hadamard_transform: fast_hadamard_transform_test.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

.PHONY: clean

clean:
	rm -f fast_hadamard_transform
