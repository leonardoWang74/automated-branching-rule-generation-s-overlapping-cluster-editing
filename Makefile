.PHONY: all help readme verify clusterEditing run sOverlappingClusterEditing clean

CXX = g++
CXXFLAGS = -O3 -std=c++17 -Wall -Wextra

all: run

help:
	cat README.md
readme: help

out:
	-mkdir out

Graph.o: Graph.cpp Graph.h
	$(CXX) $(CXXFLAGS) -c Graph.cpp

# nauty-geng plugin: enumerate non-isomorphic graphs of a given order
NAUTY_VERSION = nauty2_9_3
NAUTY_GENG = ./$(NAUTY_VERSION)/geng

$(NAUTY_GENG):
	wget --no-check-certificate https://pallini.di.uniroma1.it/$(NAUTY_VERSION).tar.gz
	tar xvzf $(NAUTY_VERSION).tar.gz && cd $(NAUTY_VERSION) && ./configure && make

# branchingAutomated script
branchingAutomated.o: branchingAutomated.cpp Graph.h
	$(CXX) $(CXXFLAGS) -c branchingAutomated.cpp

branchingAutomated-compile: branchingAutomated.o Graph.o out
	$(CXX) $(CXXFLAGS) branchingAutomated.o Graph.o -o out/branchingAutomated

TIME ?= /usr/bin/time -v

verify: branchingAutomated-compile $(NAUTY_GENG)
	# fast algorithm tests -> use only one worker
	$(NAUTY_GENG) -cq 3 | $(TIME) ./out/branchingAutomated -p 1 -s 1 -u > result-s1-n3.log 2>&1
	$(NAUTY_GENG) -cq 4 | $(TIME) ./out/branchingAutomated -p 1 -s 1 -u > result-s1-n4.log 2>&1
	$(NAUTY_GENG) -cq 5 | $(TIME) ./out/branchingAutomated -p 1 -s 1 -u > result-s1-n5.log 2>&1

	# on n=6 every worker uses around 5GB of memory
	$(NAUTY_GENG) -cq 6 | $(TIME) ./out/branchingAutomated -p 4 -s 1 -u > result-s1-n6.log 2>&1

	# on n=7 every worker uses around 30GB of memory -> use only one single worker
	# $(NAUTY_GENG) -cq 7 | ./out/branchingAutomated -p 1 -s 1 -u > result-s1-n7.log 2>&1

# alias for verify
clusterEditing: verify

run: branchingAutomated-compile $(NAUTY_GENG)
	$(NAUTY_GENG) -cq 5 | $(TIME) ./out/branchingAutomated -p 1 -s 2 -u > result-s2-n5.log 2>&1
	$(NAUTY_GENG) -cq 6 | $(TIME) ./out/branchingAutomated -p 4 -s 2 -u > result-s2-n6.log 2>&1

	# $(NAUTY_GENG) -cq 7 | ./out/branchingAutomated -p 1 -s 2 -u  > result-s2-n7.log 2>&1

# run s=2 and n=6 while measuring average CPU clock speed and memory usage (requires sudo and turbostat from linux-tools)
run-turbostat: branchingAutomated-compile $(NAUTY_GENG)
	sudo turbostat --quiet sh -c '$(NAUTY_GENG) -cq 6 | $(TIME) ./out/branchingAutomated -p 4 -s 2 -u' > result-s2-n6.log

# alias for run
sOverlappingClusterEditing: run

clean:
	rm -rf out
	rm -f *.o
