# Build: make            (Linux, g++ with OpenMP)
#        make OMP=0      (no OpenMP, e.g. Apple clang)
CXX ?= g++
OMP ?= 1
CXXFLAGS = -O3 -march=native -std=c++17 -Wall -Wextra -Wno-sign-compare
ifeq ($(OMP),1)
CXXFLAGS += -fopenmp
endif

BIN = bin/htq bin/gen
HDR = src/common.hpp src/td.hpp src/sig.hpp src/baselines.hpp

all: $(BIN)

bin/htq: src/htq.cpp $(HDR)
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) -o $@ src/htq.cpp

bin/gen: src/gen.cpp src/common.hpp
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) -o $@ src/gen.cpp

test: all
	bash tests/run_tests.sh

clean:
	rm -rf bin

.PHONY: all test clean model-tools

model-tools: bin/nmf bin/reference_common_fill bin/semantics_fixed bin/model_case

bin/nmf: model/nmf.cpp model/common.hpp
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) -o $@ $<

bin/reference_common_fill bin/semantics_fixed bin/model_case: bin/%: model/%.cpp model/common.hpp model/td.hpp model/sig.hpp model/baselines.hpp
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) -o $@ $<
