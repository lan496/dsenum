CXX	= g++
INCDIR = -I ../../extern/TdZdd/include -I ..
WARNINGFLAGS = -Wall -Wcast-qual -Winit-self -Wlogical-op -Wmissing-include-dirs -Wredundant-decls
CXXFLAGS	= -fopenmp $(WARNINGFLAGS) -std=c++14 $(INCDIR)

.DEFAULT_GOAL := test

.PONY: all test clean

all: test clean

test:
	$(CXX) $(CXXFLAGS) -O3 test_choice.cpp -o choice.out
	./choice.out
	$(CXX) $(CXXFLAGS) -O3 test_combination.cpp -o combination.out
	./combination.out
	$(CXX) $(CXXFLAGS) -O3 test_simpath.cpp -o simpath.out
	./simpath.out
	$(CXX) $(CXXFLAGS) -O3 test_stpath.cpp -o stpath.out
	./stpath.out
	$(CXX) $(CXXFLAGS) -O3 test_spanning_forest.cpp -o forest.out
	./forest.out
	$(CXX) $(CXXFLAGS) -O3 test_permutation.cpp -o permutation.out
	./permutation.out
	$(CXX) $(CXXFLAGS) -O3 test_superperiodic.cpp -o superperiodic.out
	./superperiodic.out
	$(CXX) $(CXXFLAGS) -O3 test_structure_enumeration.cpp -o structure.out
	./structure.out
	echo "[success]"

clean:  ## clean
	$(RM) *.out

help: ## Self-documented Makefile
		@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) \
			| sort \
			| awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}'
