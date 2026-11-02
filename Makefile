.PHONY: all help readme verify clusterEditing run sOverlappingClusterEditing clean

CXX = g++
CXXFLAGS = -O3 -std=c++17 -Wall -Wextra

all: run

help:
	cat README.md
readme: help

out-directory:
	-mkdir out

Graph.o: Graph.cpp Graph.h
	$(CXX) $(CXXFLAGS) -c Graph.cpp

# nauty-geng plugin: enumerate non-isomorphic graphs of a given order
NAUTY_VERSION = nauty2_9_3
$(NAUTY_VERSION):
	wget --no-check-certificate https://pallini.di.uniroma1.it/$(NAUTY_VERSION).tar.gz
	tar xvzf $(NAUTY_VERSION).tar.gz && cd $(NAUTY_VERSION) && ./configure && make

NAUTY_GENG = ./$(NAUTY_VERSION)/geng

# branchingAutomated script
branchingAutomated.o: branchingAutomated.cpp Graph.h
	$(CXX) $(CXXFLAGS) -c branchingAutomated.cpp

branchingAutomated-compile: branchingAutomated.o Graph.o out-directory
	$(CXX) $(CXXFLAGS) branchingAutomated.o Graph.o -o out/branchingAutomated

verify: branchingAutomated-compile $(NAUTY_VERSION)
	# ./out/branchingAutomated -t

	# fast algorithm tests -> use only one worker
	$(NAUTY_GENG) -cq 3 | /usr/bin/time -v ./out/branchingAutomated -p 1 -s 1 -u > result-s1-n3.log
	$(NAUTY_GENG) -cq 4 | /usr/bin/time -v ./out/branchingAutomated -p 1 -s 1 -u > result-s1-n4.log
	$(NAUTY_GENG) -cq 5 | /usr/bin/time -v ./out/branchingAutomated -p 1 -s 1 -u > result-s1-n5.log

	# on n=6 every worker uses around 5GB of memory
	$(NAUTY_GENG) -cq 6 | /usr/bin/time -v ./out/branchingAutomated -p 4 -s 1 -u > result-s1-n6.log

	# on n=7 every worker uses around 30GB of memory -> use only one single worker
	# $(NAUTY_GENG) -cq 7 | ./out/branchingAutomated -p 1 -s 1 -u > result-s1-n7.log

# alias for verify
clusterEditing: verify

run: branchingAutomated-compile $(NAUTY_VERSION)
	$(NAUTY_GENG) -cq 5 | /usr/bin/time -v ./out/branchingAutomated -p 1 -s 2 -u > result-s2-n5.log

	# measure average CPU clock speed and memory usage
	sudo turbostat --quiet sh -c '$(NAUTY_GENG) -cq 6 | /usr/bin/time -v ./out/branchingAutomated -p 4 -s 2 -u' > result-s2-n6.log

	# $(NAUTY_GENG) -cq 7 | ./out/branchingAutomated -p 1 -s 2 -u  > result-s2-n7.log

# alias for run
sOverlappingClusterEditing: run

clean:
	rm -rf out
	rm -f *.o
