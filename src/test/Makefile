CXX	= g++
INCDIR = -I ../../extern/TdZdd/include -I ..
WARNINGFLAGS = -Wall -Wcast-qual -Winit-self -Wlogical-op -Wmissing-include-dirs -Wredundant-decls
CXXFLAGS	= -fopenmp $(WARNINGFLAGS) -std=c++14 $(INCDIR)

.DEFAULT_GOAL := test

.PONY: all test clean

all: test clean

test:
	$(CXX) $(CXXFLAGS) -O3 test_simpath.cpp -o simpath.out
	./simpath.out
	$(CXX) $(CXXFLAGS) -O3 test_stpath.cpp -o stpath.out
	./stpath.out

clean:  ## clean
	$(RM) *.out


help: ## Self-documented Makefile
		@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) \
			| sort \
			| awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}'
