# Automated Branching Rule Generation for s-Overlapping Cluster Editing
is an implementation based on the framework by [Gramm Guo Hueffner Niedermeier 2004 - Automated generation of search tree algorithms for hard graph modiﬁcation problems], that can calculate optimal branching rules for s-Overlapping Cluster Editing and Cluster Editing.

### How to run
Requirements: installed `make, python3`, `g++`, `gcc, wget` commands, bash commands like `cat`, `cd`, and `mkdir`, 8 GB free disk space, and 6 GB memory (RAM) for $n \leq 6$ per worker.

Notable Makefile targets:
- Run `make clusterEditing` to run the branching rule generation algorithm for $s=1$ (cluster editing) on $n \in [3,6]$. The results are saved in files `result-s1-n*.log`.
- Run `make run` to run the branching rule generation algorithm for $s=2$ on $n \in [5,6]$. The results are saved in files `result-s2-n*.log`.

The line "Maximum resident set size (kbytes)" (if you have the `time` package installed) is the maximum RAM usage (per worker).

The program (after compiling with `make branchingAutomated-compile`) can also be run directly with `./nauty2_9_3/geng -cq 5 | /usr/bin/time -v ./out/branchingAutomated -p 1 -s 1 -u > result-s1-n5.log` where `p` is the number of processes, `s` is $s$, and `-u` turns on the greedy scheme.
