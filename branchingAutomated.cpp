/*

Automated generation of branching rules, based on
[Gramm Guo Hueffner Niedermeier 2004 - Automated generation of search tree algorithms for hard graph modiﬁcation problems]

Run using nauty geng
nauty-geng -c 6 | ./branchingAutomated
*/
#include <unistd.h>
#include <sys/wait.h>

#include <cerrno>
#include <cstring>
#include <ctime>
#include <iostream>
#include <string>
#include <iomanip>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <functional>
#include <cmath>
#include <unordered_map>
#include <unordered_set>
#include "Graph.h"

#define ANNOTATION_UNDEF 0
#define ANNOTATION_PERMANENT 1
#define ANNOTATION_FORBIDDEN 2

// optional flags
#define DEDUPLICATION_ON 1 // if defined, use deduplication (memoization of \pi -> B)

// global variables
bool debug = false;
double branchingNumberWorstCase = 0;
// if > 0 and found a branching number > this number, try the optimal algorithm
double branchingNumberTreshold = -1;
int s = 1;
int depth_cutoff = -1;
int workerId = -1;
enum AlgorithmTypes {
    GRAMM_ET_AL,
    FORBIDDEN_SUBGRAPH_SINGLE_GREEDY_ALL_PAIRS,
    FORBIDDEN_SUBGRAPH_SINGLE_GREEDY_SINGLE_PAIR,
};
std::string algorithmTypeToString(AlgorithmTypes algorithmType) {
    switch(algorithmType) {
        case GRAMM_ET_AL: return "GRAMM_ET_AL";
        case FORBIDDEN_SUBGRAPH_SINGLE_GREEDY_ALL_PAIRS: return "FORBIDDEN_SUBGRAPH_SINGLE_GREEDY_ALL_PAIRS";
        case FORBIDDEN_SUBGRAPH_SINGLE_GREEDY_SINGLE_PAIR: return "FORBIDDEN_SUBGRAPH_SINGLE_GREEDY_SINGLE_PAIR";
    }
    return "UNKNOWN";
}

// greedy algorithms only give an upper bound for the branching number of a graph.
// GRAMM_ET_AL is the only optimal one
bool isGreedy(AlgorithmTypes algorithmType) {
    return algorithmType != GRAMM_ET_AL;
}

// return `str * factor` = `factor` repetitions of `str`
std::string stringMultiply(std::string str, int factor){
    std::string ret = "";
    ret.reserve(str.size() * factor);
    for(int i=0; i<factor; ++i) ret.append(str);
    return ret;
}

// lossless text representation of a double (std::to_string only keeps 6 decimals)
std::string doubleToString(double x) {
    std::ostringstream stream;
    stream << std::setprecision(17) << x;
    return stream.str();
}

struct Pipe {
    int fd[2];
};

// type of pi: store annotations of pairs of vertices + the number of edits required
struct piType {
    std::vector<std::vector<int>> annotations;
    // number of edits = cost of this annotation
    int number_edits = 0;
    // number of fixed pairs of vertices
    int number_fixed = 0;
    // TRUE <=> there is no induced forbidden subgraph
    bool no_forbidden = false;
    // TRUE <=> this branch is infeasible (an induced forbidden subgraph with only fixed pairs of vertices exists)
    bool infeasible = false;
};
bool piNumberEditsSort(piType& a, piType& b) {
    return a.number_edits < b.number_edits;
}
typedef std::vector<int> RuleType;

#ifdef DEDUPLICATION_ON
unsigned long piHash(const piType& pi) {
    unsigned long hash = 0;
    for(unsigned int v=0; v<pi.annotations.size(); ++v) {
        auto& row = pi.annotations[v];
        for(unsigned int w=v+1; w<row.size(); ++w) {
            auto value = row[w];
            if(value == 0) continue;
            hash += value * w;
        }
        hash *= 1033;
    }
    return hash;
}
bool piEqualsAnnotations(const piType& a, const piType& b) {
    // if(a.no_forbidden != b.no_forbidden) return false;
    // if(a.number_edits != b.number_edits) return false;
    if(a.annotations.size() != b.annotations.size()) return false;

    for(unsigned int v=0; v<a.annotations.size(); ++v) {
        auto& rowA = a.annotations[v];
        auto& rowB = b.annotations[v];
        if(rowA.size() != rowB.size()) return false;

        for(unsigned int w=v+1; w<rowA.size(); ++w) {
            if(rowA[w] != rowB[w]) return false;
        }
    }

    return true;
}

// saves hash -> [list of pairs (pi, calculated rules B)]
auto deduplication_map = std::unordered_map<unsigned long, std::vector<
        std::pair<piType, std::vector<RuleType>>
    >>();
#endif

int forbiddenSubgraphSizeUpperBound() {
    if(s <= 4) return 2*s+1;

    // proven ell^{3/2} bound
    return std::min((int)std::floor(s+2 + std::sqrt(2)*std::pow(s+1, 1.5)), 1+s*(s+1));
}

// apply reduction rules
std::vector<std::pair<int,int>> br_reduce(Graph* G, MaximalCliquesInfo& cliqueInfo, piType& pi, const std::vector<int>& last_edit_vertices) {
    std::vector<std::pair<int,int>> edits_made = {};

    // check if there is a vertex in >s maximal cliques
    if(cliqueInfo.vertexInMoreThanSCliques == -1) {
        pi.no_forbidden = true;
        return edits_made;
    }

    ////////////////////////////////////////////////
    // check if there is a FIXED forbidden subgraph,
    // that is, a forbidden subgraph where all pairs of vertices are fixed (PERMANENT/FORBIDDEN)
    // then, the instance is infeasible
    if(!last_edit_vertices.empty()) {
        bool found = false;

        // lambda to receive the indices and call the solver function
        auto lambda = [G, &cliqueInfo, &found, &pi, &last_edit_vertices](size_t, std::vector<size_t> indices) {
            // check no vertex = last_edit_vertices
            for(auto v : indices) {
                for(auto w : last_edit_vertices) {
                    if(v == (size_t)w) return true;
                }
            }
            // add all last_edit_vertices to our chosen set
            for(auto w : last_edit_vertices) {
                indices.push_back(w);
            }

            // check all edits inside are fixed
            for(size_t i=0; i<indices.size(); ++i) {
                for(size_t j=i+1; j<indices.size(); ++j) {
                    int v,w;
                    if(indices[i] < indices[j]) {
                        v = indices[i];
                        w = indices[j];
                    }
                    else {
                        w = indices[i];
                        v = indices[j];
                    }

                    if(pi.annotations[v][w] == ANNOTATION_UNDEF) return true; // try a different subset
                }
            }

            // make list of the chosen vertices (faster in a new loop: most indices will contain an UNDEF annotation)
            std::vector<int> chosenVertices = {};
            chosenVertices.reserve(indices.size());
            for(auto v : indices) {
                chosenVertices.push_back((int)v);
            }

            // create subgraph
            Graph Gprime = G->getSubgraph(chosenVertices);

            // check subgraph contains a vertex
            auto cliqueInfoPrime = Gprime.getMaximalCliques(s, true);

            // subgraph contains a forbidden subgraph => infeasible
            if(cliqueInfoPrime.vertexInMoreThanSCliques >= 0) {
                found = true;
                return false; // stop looking for subsets
            }
            
            // try a different subset
            return true;
        };

        // size without last_edit_vertices: they are always added
        size_t sizeBound = std::min((size_t)forbiddenSubgraphSizeUpperBound() - last_edit_vertices.size(), (size_t)G->n());

        // check if it is possible (with pi.number_fixed) to have a fixed subgraph of this size
        // a graph of size n has $\binom{n}{2} = n(n-1)/2 = p$ pairs of vertices.
        // we have $p$ = pi.number_fixed pairs of fixed pairs of vertices.
        // Solving for $n$: $n^2-n-2p = 0$
        // $n = 1/2 +/- \sqrt{1/4 + 2p}$
        // $n \leq 1/2 + \sqrt{1/4 + 2p}$
        // subtract 2 since we fix two vertices of the enumeration
        double sizeLargestByEdits = 0.5 + std::sqrt(0.25 + 2 * pi.number_fixed) - 2;
        if(sizeBound > sizeLargestByEdits) {
            sizeBound = (size_t)( std::floor(sizeLargestByEdits) );
        }

        // look for fixed induced forbidden subgraphs
        for(size_t size=1; size <= sizeBound; ++size) {
            SubsetsOfSizeLoop(G->n(), size, lambda);

            if(found) {
                pi.infeasible = true;
                return edits_made;
            }
        }
    }

    return edits_made;
}

// prune the branching rules in $B$ by subsumption
void br_subsumption_prune(std::vector<RuleType>& B) {
    // sort the branching vector ascending
    for(auto& A : B) {
        std::sort(A.begin(), A.end());
    }

    bool restart = false;
    for(size_t i1=0; i1<B.size(); ++i1) {
        if(restart) {
            restart = false;
            i1 = 0;
        }

        auto& A1 = B[i1];
        if(A1.empty()) continue; // do not subsume empty rules

        for(size_t j2=i1+1; j2<B.size(); ++j2) {
            auto& A2 = B[j2];
            if(A2.empty()) continue; // do not subsume empty rules

            // test subsumption of A2
            if(A1.size() <= A2.size()) {
                bool subsumption = true;
                for(size_t i=0; i<A1.size(); ++i) {
                    if(A1[i] >= A2[i]) continue;
                    subsumption = false;
                    break;
                }
                if(subsumption) {
                    B.erase(B.begin()+j2);
                    --j2;
                    continue; // do not check for subsumption of A1
                }
            }

            // test subsumption of A1
            if(A1.size() >= A2.size()) {
                bool subsumption = true;
                for(size_t i=0; i<A2.size(); ++i) {
                    if(A2[i] >= A1[i]) continue;
                    subsumption = false;
                    break;
                }
                if(subsumption) {
                    B.erase(B.begin()+i1);
                    // re-do the loop with a new A1
                    if(i1>0) {
                        --i1;
                    }
                    else {
                        // start the whole loop anew
                        restart = true;
                    }
                    break; // break the j2 loop (A1 does not exist anymore)
                }
            }
        }
    }
}

// concatenate and prune branching rules (only by subsumption)
std::vector<RuleType> br_concatenate(std::vector<RuleType>& B1, std::vector<RuleType>& B2) {
    std::vector<RuleType> B = {};
    B.reserve(8192);
    
    br_subsumption_prune(B1);
    br_subsumption_prune(B2);

    // create the cartesian product B1 \cross B2
    for(auto& A1 : B1) {
        for(auto& A2 : B2) {
            // A1 \cup A2
            RuleType A(A1);
            for(auto& pi : A2) A.push_back(pi);

            B.push_back(A);
        }

        // prune early for less memory usage
        if(B.size() > 4096) br_subsumption_prune(B);
    }

    br_subsumption_prune(B);
    return B;
}

std::string ruleTypeToBranchingVector(RuleType& A) {
    std::string ret = "(";
    bool first = true;
    for (const auto& pi : A) {
        if(!first) ret += ",";
        else first = false;

        ret += std::to_string(pi);
    }
    return ret + ")";
}

// calculate the branching vector corresponding to the given branching rule
std::vector<int> branching_rule_to_vector(const RuleType& A) {
    std::vector<int> ret = {};
    ret.reserve(A.size());
    for(const auto& pi : A) {
        ret.push_back(pi);
    }
    return ret;
}

// calculate the branching number corresponding to the given branching rule
double branching_rule_to_number(const RuleType& A){
    if (A.empty()) throw std::runtime_error("Empty branching rule A.");
    for (const auto& pi : A) {
        if(pi <= 0) throw std::runtime_error("pi contains no edits.");
    }

    // binary search on the characteristic polynomial f(z) = \sum_i z^{-d_i} = 1
    double low = 1.0, high = static_cast<double>(A.size());
    for (int it = 0; it < 150; ++it) {
        const double mid = 0.5 * (low + high);

        // f(z) = \sum_i z^{-d_i}
        double f = 0;
        for (const auto& di : A) f += std::pow(mid, -di);

        // target: f(z) = 1
        if(f > 1.0) {
            low = mid;
        }
        else {
            high = mid;
        }
    }

    // return found z value
    return 0.5 * (low + high);
}

// find the smallest branching number corresponding to the branching rules in the given set of branching rules
RuleType& branching_rule_best(std::vector<RuleType>& B) {
    double br = 99999;
    if(B.empty()) {
        throw std::runtime_error("B is empty.");
    }

    size_t best = 0;
    for(size_t i=0; i<B.size(); ++i) {
        auto& A = B[i];
        auto brHere = branching_rule_to_number(A);
        if(debug) {
            std::cout << "\t" << __FILE__<<":"<<__LINE__
                <<" branching vector="<<Graph::vector_tostring(branching_rule_to_vector(A))
                <<" branching factor br="<<Graph::vector_tostring(brHere)
                <<"\n";
        }
        if(brHere >= br) continue;
        br = brHere;
        best = i;
    }
    if(debug) {
        std::cout << __FILE__<<":"<<__LINE__
            <<" best branching vector="<<Graph::vector_tostring(branching_rule_to_vector(B[best]))
            <<" branching factor br="<<Graph::vector_tostring(br)
            <<"\n";
    }
    return B[best];
}

// toggle a list of edits on the graph G
void G_edit_toggle(Graph* G, std::vector<std::pair<int,int>>& edits) {
    for(const auto& edit : edits) {
        auto v = edit.first;
        auto w = edit.second;

        if(G->edge_has(v,w)) {
            G->edge_remove(v,w);
        } else {
            G->edge_add(v,w);
        }
    }
}

// algorithm for finding a k-clique-minimal induced subgraph G[V_F]
// find V_F in the given startingCliques (=indices of cliques of cliqueInfo.cliqueSets)
std::vector<int> findForbiddenSubgraphkCliqueMinimal(Graph* G, MaximalCliquesInfo& cliqueInfo, int k, int vInAllCliques, const std::vector<size_t>& startingCliques) {
    /////////////////////////////////////////
    // trivial distinguisher set $S$, of the given startingCliques
    std::unordered_set<int> S(k*(k-1));

    size_t iterMax = std::min( (size_t)k, startingCliques.size()); // only need to distinguish $k$ cliques
    for(size_t i=0; i<iterMax; ++i) {
        for(size_t j=i+1; j<iterMax; ++j) {
            auto& cliqueA = cliqueInfo.cliqueSets[ startingCliques[i] ];
            auto& cliqueB = cliqueInfo.cliqueSets[ startingCliques[j] ];

            auto AwithoutB = Graph::set_difference(cliqueA, cliqueB);
            auto BwithoutA = Graph::set_difference(cliqueB, cliqueA);
            
            int a = *AwithoutB.begin();
            auto& neighbors = G->neighbors(a);
            auto potentialB = Graph::set_difference(cliqueB, std::unordered_set<int>(neighbors.begin(), neighbors.end()));
            int b = *potentialB.begin();

            S.insert(a);
            S.insert(b);
        }
    }

    // subgraph G' = G[S]
    Graph GprimeValue = G->getSubgraph(std::vector<int>(S.begin(), S.end()));
    Graph* Gprime = &GprimeValue;

    /////////////////////////////////////////
    // calculate maximal cliques in $S$
    auto cliqueInfoPrime = Gprime->getMaximalCliques(s, false);

    // create clique sets
    for(auto& clique : cliqueInfoPrime.cliqueList) {
        cliqueInfoPrime.cliqueSets.push_back(std::unordered_set<int>(clique.begin(), clique.end()));
    }

    /////////////////////////////////////////
    // reduce $G' = G[S]$ to a $k$-clique-minimal induced subgraph
    auto& C_F = cliqueInfoPrime.cliqueSets; // do not need to deepcopy cliqueInfoPrime.cliqueSets
    auto V_F = std::unordered_set<int>(S.size()); // unordered_set for O(1) removal
    for(int v=0; v<Gprime->n_signed(); ++v) V_F.insert(v); // IDs in the subgraph are 0 .. |V(G')|

    if((unsigned int)k > C_F.size()) throw std::runtime_error("Algorithm called findForbiddenSubgraphkCliqueMinimal with k > startingCliques "+G->to_graph6()+".");

    // remove loop
    while(true) {
        size_t cliqueCount = C_F.size();
        int dBound = 1 + (int)C_F.size() - k;

        // D(v) for every vertex. D[v][i] = v is single distinguisher for clique i.
        // vertexDpairs[v][cliqueA] = true means that v has already incremented D(v) because of cliqueA
        std::vector<std::vector<bool>> vertexDpairs(Gprime->n(), std::vector<bool>(cliqueCount, false));

        // |D(v)| for every vertex = count of cliques for which $v$ is the single distinguisher
        // = counting the entries of D(v) e.g. D(v) = {A,B,C} to |D(v)| = 3
        std::vector<int> vertexDcount(Gprime->n(), 0);

        // pairwise differences -> D(v)
        // note that the clique sets contain at most the vertices V_F
        for(size_t i=0; i<cliqueCount; ++i) {
            for(size_t j=i+1; j<cliqueCount; ++j) {
                auto& cliqueA = C_F[i];
                auto& cliqueB = C_F[j];

                auto AwithoutB = Graph::set_difference(cliqueA, cliqueB);
                auto BwithoutA = Graph::set_difference(cliqueB, cliqueA);
                
                // single distinguisher exists
                if(AwithoutB.size() == 1) {
                    int v = *AwithoutB.begin();

                    // increment D(v) if D(v) does not yet contain A (count every vertex+clique pair only once)
                    if(!vertexDpairs[v][i]){
                        vertexDpairs[v][i] = true;
                        ++vertexDcount[v];
                    }
                }

                // single distinguisher exists
                if(BwithoutA.size() == 1) {
                    int v = *BwithoutA.begin();

                    // increment D(v) if D(v) does not yet contain B (count every vertex+clique pair only once)
                    if(!vertexDpairs[v][j]){
                        vertexDpairs[v][j] = true;
                        ++vertexDcount[v];
                    }
                }
            }
        }

        // look for vertex v \in V_F with |D(v)| < 1+|C_F|-k
        int x = -1;
        for(auto w : V_F) {
            if(vertexDcount[w] >= dBound) continue;
            x = w;
            break;
        }

        // no vertex found: G[V_F] is k-clique-minimal
        if(x < 0) break;

        // remove x: update C_F and V_F
        V_F.erase(x);

        // reverse loop so removing entries preserves (lower) indices
        size_t i=cliqueCount;
        for(size_t j=0; j<cliqueCount; ++j) {
            --i; // size_t is unsigned so loop variable counts up

            // clique A_i \in D(x) => remove from C_F
            if(vertexDpairs[x][i]) {
                C_F.erase(C_F.begin() + i);
            }
            // clique A_i \not \in D(x) => keep and remove x from the clique
            else {
                C_F[i].erase(x);
            }
        }

        // we recalculate D(v) above
    }

    // after the loop: G[V_F] is k-clique-minimal.
    // return the original (mapped) vertices V_F \cup \{v\}
    std::vector<int> subgraph = {};
    subgraph.reserve(V_F.size()+1);
    subgraph.push_back(vInAllCliques);
    for(auto w : V_F) {
        subgraph.push_back(Gprime->id_get(w));
    }
    return subgraph;
}

// find a set of vertices inducing a minimal forbidden subgraph
std::vector<int> findForbiddenSubgraph(Graph* G, MaximalCliquesInfo& cliqueInfo) {
    // no vertex in s+1 maximal cliques
    if(cliqueInfo.vertexInMoreThanSCliques < 0) return {};

    int v = cliqueInfo.vertexInMoreThanSCliques;

    auto& vertexCliques = cliqueInfo.vertexCliques[v];
    
    // s=1: trivial distinguishers = P_3
    if(s == 1) {
        auto& cliqueA = cliqueInfo.cliqueSets[vertexCliques[0]];
        auto cliqueB = cliqueInfo.cliqueSets[vertexCliques[1]];

        auto AwithoutB = Graph::set_difference(cliqueA, cliqueB);
        auto BwithoutA = Graph::set_difference(cliqueB, cliqueA);
        
        auto a = *AwithoutB.begin();
        auto& neighbors = G->neighbors(a);
        auto potentialB = Graph::set_difference(cliqueB, std::unordered_set<int>(neighbors.begin(), neighbors.end()));
        int b = *potentialB.begin();

        return {v,a,b};
    }

    // otherwise, find a (s+1)-clique-minimal induced subgraph in the cliques of v
    return findForbiddenSubgraphkCliqueMinimal(G, cliqueInfo, s+1, v, vertexCliques);
}

// compute branching rules
// based on [Gramm et al 2004 - Automated generation of search tree algorithms for hard graph modiﬁcation problems]
std::vector<RuleType> br_compute(Graph* G, piType& piOriginal, AlgorithmTypes algorithmType, const std::vector<int>& last_edit_vertices={}) {
    std::vector<RuleType> B = {};

    // copy pi
    piType pi = piOriginal;

#ifdef DEDUPLICATION_ON
    // check deduplication
    auto hash = piHash(pi);
    auto hashIt = deduplication_map.find(hash);
    auto hashCollision = hashIt != deduplication_map.end();

    // hash collision: check if real collision or only hash collision
    if(hashCollision) {
        for(const auto& dedup_pair : hashIt->second) {
            // real collision: already had this branch
            if(piEqualsAnnotations(pi, dedup_pair.first)) {
                // return the already calculated rule set
                return dedup_pair.second;
            }
        }
    }
#endif

    std::vector<int> verticesToIteratePairsOn = {};
    std::vector<std::pair<int,int>> edge_edits;

    // new block to free cliqueInfo after
    {
        // get all cliques + create clique sets
        auto cliqueInfo = G->getMaximalCliques(s, false);
        for(auto& clique : cliqueInfo.cliqueList) {
            cliqueInfo.cliqueSets.push_back(std::unordered_set<int>(clique.begin(), clique.end()));
        }
        
        // apply reduction rules
        edge_edits = br_reduce(G, cliqueInfo, pi, last_edit_vertices);

        // infeasible: found forbidden subgraph where no edit is possible = all pairs of vertices are fixed => no children branch can fix this
        if(pi.infeasible) {
            #ifdef DEDUPLICATION_ON
            // update deduplication map with (pi,output)
            deduplication_map[hash].push_back({pi, {}});
            #endif

            if(debug) {
                std::cout << stringMultiply("    ", pi.number_fixed) << __FILE__<<":"<<__LINE__
                    <<" current graph="<<G->to_graph6()
                    <<" edits="<<Graph::vector_tostring(pi.number_edits)
                    <<" INFEASIBLE "
                    <<"\n";
            }

            // undo edits done in br_reduce
            G_edit_toggle(G, edge_edits);
            return {};
        }

        // find branching pairs
        if(!pi.no_forbidden) {
            // check all pairs
            if(algorithmType == GRAMM_ET_AL) {
                verticesToIteratePairsOn.reserve(G->n());
                for(int v=0; v<G->n_signed(); ++v) verticesToIteratePairsOn.push_back(v);
            }
            // branch on all pairs of any forbidden subgraph
            else if(algorithmType == FORBIDDEN_SUBGRAPH_SINGLE_GREEDY_ALL_PAIRS || algorithmType == FORBIDDEN_SUBGRAPH_SINGLE_GREEDY_SINGLE_PAIR) {
                verticesToIteratePairsOn = findForbiddenSubgraph(G, cliqueInfo);
                std::sort(verticesToIteratePairsOn.begin(), verticesToIteratePairsOn.end());
            }

            // !pi.noforbidden -> we should find vertices to branch on
            if(verticesToIteratePairsOn.empty()) {
                throw std::runtime_error("Did not find vertices to iterate on graph="+G->to_graph6()+" but there should be a forbidden subgraph.");
            }
        }
        else if(debug) {
            std::cout << stringMultiply("    ", pi.number_fixed) << __FILE__<<":"<<__LINE__
                <<" current graph="<<G->to_graph6()
                <<" edits="<<Graph::vector_tostring(pi.number_edits)
                <<" SOLUTION"
                <<"\n";
        }
                    
    }

    // branch: find branching rules = try pairs of vertices to branch on
    if(!verticesToIteratePairsOn.empty()) {
        // set to true if there is a non-annotated pair
        bool branched = false;

        for(size_t i=0; i<verticesToIteratePairsOn.size(); ++i) {
            auto v = verticesToIteratePairsOn[i];
            for(size_t j=i+1; j<verticesToIteratePairsOn.size(); ++j) {
                auto w = verticesToIteratePairsOn[j];

                if(pi.annotations[v][w] != ANNOTATION_UNDEF) continue;
                branched = true;

                if(debug) {
                    std::cout << stringMultiply("    ", pi.number_fixed) << __FILE__<<":"<<__LINE__
                        <<" current graph="<<G->to_graph6()
                        <<" edits="<<Graph::vector_tostring(pi.number_edits)
                        <<" branching on v="<<Graph::vector_tostring(v)
                        <<" w="<<Graph::vector_tostring(w)
                        <<" verticesToIteratePairsOn="<<Graph::vector_tostring(verticesToIteratePairsOn)
                        <<"\n";
                }

                // permanent
                piType pi1 = pi; // copy
                pi1.annotations[v][w] = ANNOTATION_PERMANENT;
                ++pi1.number_fixed;

                // forbidden
                piType pi2 = pi; // copy
                pi2.annotations[v][w] = ANNOTATION_FORBIDDEN;
                ++pi2.number_fixed;

                std::vector<RuleType> B1 = {};
                std::vector<RuleType> B2 = {};
                if(G->edge_has(v,w)) {
                    // edit graph for pi2
                    ++pi2.number_edits;
                    G->edge_remove(v,w);
                    B2 = br_compute(G, pi2, algorithmType, {v,w});
                    G->edge_add(v,w); // undo edit

                    // graph is unchanged for pi1
                    B1 = br_compute(G, pi1, algorithmType, {v,w});
                }
                else {
                    // edit graph for pi1
                    ++pi1.number_edits;
                    G->edge_add(v,w);
                    B1 = br_compute(G, pi1, algorithmType, {v,w});
                    G->edge_remove(v,w); // undo edit

                    // graph is unchanged for pi2
                    B2 = br_compute(G, pi2, algorithmType, {v,w});
                }

                // both branches are infeasible => this branch is infeasible
                if(B1.empty() && B2.empty()) {
                    /*std::cout << "\t" <<__FILE__<<":"<<__LINE__
                        <<" v="<<v
                        <<" w="<<w
                        <<" both branches infeasible for "
                        <<"\n";*/

                    #ifdef DEDUPLICATION_ON
                    // update deduplication map with (pi,output)
                    deduplication_map[hash].push_back({pi, {}});
                    #endif

                    if(debug) {
                        std::cout << stringMultiply("    ", pi.number_fixed) << __FILE__<<":"<<__LINE__
                            <<" current graph="<<G->to_graph6()
                            <<" edits="<<Graph::vector_tostring(pi.number_edits)
                            <<" INFEASIBLE (both children infeasible) "
                            <<"\n";
                    }
                    
                    // undo edits done in br_reduce
                    G_edit_toggle(G, edge_edits);
                    return {}; 
                }
                // one branch is infeasible === "reduction rule" => we do not have to branch on the infeasible branch
                else if(B1.empty() || B2.empty()) {
                    auto& BnonEmpty = B1.empty() ? B2 : B1;
                    br_subsumption_prune(BnonEmpty);
                    for(auto& Bi : BnonEmpty) {
                        B.push_back(Bi);
                    }
                }
                // both branches feasible: concatenate (cross product)
                else {
                    auto Bconcat = br_concatenate(B1, B2);
                    for(auto& Bi : Bconcat) {
                        B.push_back(Bi);
                    }
                }

                if(algorithmType == FORBIDDEN_SUBGRAPH_SINGLE_GREEDY_SINGLE_PAIR) break;
            }            
            if(branched && algorithmType == FORBIDDEN_SUBGRAPH_SINGLE_GREEDY_SINGLE_PAIR) break;
        }

        // every pair of the forbidden subgraph was already fixed -> must have been infeasible
        if(!branched) {
            throw std::runtime_error("Fully-fixed forbidden subgraph on graph="+G->to_graph6()+" and vertices="+Graph::vector_tostring(verticesToIteratePairsOn)+" not detected by br_reduce");
        }
    }

    // reduced pi itself is a rule
    if(pi.number_edits > 0) {
        B.push_back({pi.number_edits});
    }
    // nothing to do: emit an empty rule, so the cross product is the identity
    else if(pi.no_forbidden) {
        B.push_back({});
    }

    // undo edits done in br_reduce
    G_edit_toggle(G, edge_edits);

    // prune B
    br_subsumption_prune(B);

    #ifdef DEDUPLICATION_ON
    // update deduplication map with (pi,output)
    deduplication_map[hash].push_back({pi, B});
    #endif

    return B;
}


// type of saved graph info
struct graphSaveType {
    int graphsCount;
    // the graph in g6 format
    std::string line;
    double branching_number;
    // which algorithm has been run on this graph
    AlgorithmTypes algorithmRun = FORBIDDEN_SUBGRAPH_SINGLE_GREEDY_SINGLE_PAIR;
};

// saves graph number -> graph info, in order to run the optimal algorithm on the highest branching number graphs
auto graph_with_branching_number = std::unordered_map<int, graphSaveType>();

// "main" function for a given graph
double branchingAutomated(int graphsCount, std::string line, AlgorithmTypes algorithmType) {
    Graph GraphValue = Graph::parse_graph6(line);
    Graph* G = &GraphValue;
    
    // pi: UNDEF for every pair of vertices
    auto pi = piType();
    pi.annotations = std::vector<std::vector<int>>(G->n());
    for(size_t i=0; i<G->n(); ++i) {
        pi.annotations[i] = std::vector<int>(G->n(), ANNOTATION_UNDEF);
    }

    #ifdef DEDUPLICATION_ON
    deduplication_map.clear();
    #endif

    auto B = br_compute(G, pi, algorithmType);

    // infeasible: error
    if(B.empty()) {
        throw std::runtime_error("Algorithm returned infeasible for the graph "+line+" - but an infeasible instance does not exist.");
    }
    // trivial instance
    if(B.size()==1 && B[0].empty()) {
        std::cout << __FILE__<<":"<<__LINE__<<" graph "<<graphsCount<<": "<<line
            <<" trvial instance (no forbidden subgraph)"
            <<"\n";
        return 1;
    }

    auto& A_best = branching_rule_best(B);
    auto br = branching_rule_to_number(A_best);

    // save object
    auto save = graphSaveType();
    save.graphsCount = graphsCount;
    save.line = line;
    save.branching_number = br;
    save.algorithmRun = algorithmType;

    if(branchingNumberTreshold>0 && br >= branchingNumberTreshold) {
        if(algorithmType == FORBIDDEN_SUBGRAPH_SINGLE_GREEDY_SINGLE_PAIR) {
            std::cout << "\t" << __FILE__<<":"<<__LINE__<<" worker="<<workerId<<" graph "<<graphsCount<<": "<<line
                <<" found br="<<br
                <<" trying less-greedy algorithm"
                <<"\n";
            br = branchingAutomated(graphsCount, line, FORBIDDEN_SUBGRAPH_SINGLE_GREEDY_ALL_PAIRS);
            save.branching_number = br;
            save.algorithmRun = FORBIDDEN_SUBGRAPH_SINGLE_GREEDY_ALL_PAIRS;
        }
        else if(algorithmType == FORBIDDEN_SUBGRAPH_SINGLE_GREEDY_ALL_PAIRS) {
            std::cout << "\t" << __FILE__<<":"<<__LINE__<<" worker="<<workerId<<" graph "<<graphsCount<<": "<<line
                <<" found br="<<br
                <<" trying optimal proposition algorithm"
                <<"\n";
            br = branchingAutomated(graphsCount, line, GRAMM_ET_AL);
            std::cout << "\t" << __FILE__<<":"<<__LINE__<<" worker="<<workerId<<" graph "<<graphsCount<<": "<<line
                <<" improved br="<<br
                <<"\n";
            save.branching_number = br;
            save.algorithmRun = GRAMM_ET_AL;

            // update try optimal number - found a worse worst-case already
            if(br > branchingNumberTreshold) branchingNumberTreshold = br;
        }
    }

    graph_with_branching_number[graphsCount] = save;

    return br;
}

int test_branching_rule_to_number() {
    int failures = 0;
    std::vector<std::vector<int>> tests = {
        {1,2},
        {1,3,3},
        {1,2,5},
        {2,2,2},
        {1,2,1,2,5},
    };
    std::vector<double> results = {
        1.618,
        1.696,
        1.71,
        1.732,
        2.75,
    };
    for(size_t t=0; t<tests.size(); ++t) {
        auto result = results[t];

        RuleType A(tests[t]);
        auto br = branching_rule_to_number(A);
        if(std::abs(br-result) > 0.01) {
            ++failures;
            std::cout << "########## Test failed - branching_rule_to_number: "<<ruleTypeToBranchingVector(A)<<" branching number="<<br<<" != "<<result<<" ##########\n";
        } else {
            std::cout << "Test success - branching_rule_to_number: "<<ruleTypeToBranchingVector(A)<<" branching number="<<br<<" ~= "<<result<<"\n";
        }
    }
    return failures;
}

int test() {
    int failures = 0;
    failures += test_branching_rule_to_number();
    return failures;
}

// function run by the worker processes
void workerMain(int w, int workers, AlgorithmTypes algorithmType, Pipe& pipeReadyWorkers, std::vector<Pipe>& pipesWriteFromMain, std::vector<Pipe>& pipesWriteFromWorker) {
    workerId = w;
 
    // close unused pipes
    for(int j = 0; j < workers; ++j) {
        if(j == w) {
            close(pipesWriteFromMain[j].fd[1]); // close write end of own read pipe
            close(pipesWriteFromWorker[j].fd[0]); // close read end of own write pipe
        } else {
            // close both ends of pipes not belonging to self
            close(pipesWriteFromMain[j].fd[0]);
            close(pipesWriteFromMain[j].fd[1]);
 
            close(pipesWriteFromWorker[j].fd[0]);
            close(pipesWriteFromWorker[j].fd[1]);
        }
    }
    close(pipeReadyWorkers.fd[0]); // worker closes read end of the ready pipe
 
    const int fdReady = pipeReadyWorkers.fd[1];
    const int fdToMain = pipesWriteFromWorker[w].fd[1];
 
    auto timeStartWorker = TimeNow();
 
    FILE* stream = fdopen(pipesWriteFromMain[w].fd[0], "r");
    char buffer[4096];
 
    // false: phase 1 (given algorithm, answer with "ready"), true: phase 2 (optimal algorithm, answer with "result")
    bool optimalPhase = false;
 
    // run until main closes the pipe
    while(fgets(buffer, sizeof(buffer), stream)) {
        std::string msg(buffer);
        if(!msg.empty() && msg.back() == '\n') msg.pop_back();
        if(msg.empty()) continue;
 
        //////////////////////////////////////////
        // synchronisation: send all saves to main, then switch to the optimal algorithm
        if(msg == "SYNC") {
            long timeSpend = TimeDifference(timeStartWorker);
            std::cout << "worker id="<<w<<" finished phase 1 in "<<(timeSpend / 1000000.0)<<"s."
                <<" branchingNumberWorstCase="<<branchingNumberWorstCase
                <<" sending "<<graph_with_branching_number.size()<<" saves to main\n";
            std::cout.flush();
 
            for(const auto& entry : graph_with_branching_number) {
                const auto& save = entry.second;
                std::string toSend = "save "+std::to_string(save.graphsCount)
                    +" "+std::to_string((int)save.algorithmRun)
                    +" "+doubleToString(save.branching_number)
                    +" "+save.line+"\n";
                if(write(fdToMain, toSend.c_str(), toSend.size()) < 0) {
                    std::cout << "Failed to write to main from worker w="<<w<<" error:"<<std::strerror(errno)<<"\n";
                }
            }
            std::string toSend = "done "+std::to_string(timeSpend)+"\n";
            if(write(fdToMain, toSend.c_str(), toSend.size()) < 0) {
                std::cout << "Failed to write to main from worker w="<<w<<" error:"<<std::strerror(errno)<<"\n";
            }
 
            // main has all saves now
            graph_with_branching_number.clear();
            optimalPhase = true;
            continue;
        }
 
        //////////////////////////////////////////
        // graph message: "graphsCount graph6"
        size_t pos = msg.find(' ');
        if(pos == std::string::npos) {
            std::cout << "worker id="<<w<<" received unknown message: "<<msg<<"\n";
            continue;
        }
 
        int graphsCount = std::stoi(msg.substr(0, pos));
        std::string line = msg.substr(pos + 1);
        AlgorithmTypes algorithmHere = optimalPhase ? GRAMM_ET_AL : algorithmType;
 
        std::cout << __FILE__<<":"<<__LINE__
            <<" worker="<<workerId
            <<" starting on graph "<<graphsCount<<": "<<line
            <<" algorithm="<<algorithmTypeToString(algorithmHere)
            <<"\n";
 
        auto timeStartThisGraph = TimeNow();
 
        auto br = branchingAutomated(graphsCount, line, algorithmHere);
        if(!optimalPhase) branchingNumberWorstCase = std::max(branchingNumberWorstCase, br);
 
        long timeSpendThisGraph = TimeDifference(timeStartThisGraph);
        std::cout << __FILE__<<":"<<__LINE__
            <<" worker="<<workerId
            <<" graph "<<graphsCount<<": "<<line
            <<" br="<<br
            <<" time="<<(timeSpendThisGraph / 1000000.0)<<"s"
            <<" algorithm="<<algorithmTypeToString(algorithmHere)
            <<"\n";
        std::cout.flush();
 
        if(!optimalPhase) {
            // phase 1: send "ready" to main process
            std::string toSend = "ready "+std::to_string(workerId)+"\n";
            if(write(fdReady, toSend.c_str(), toSend.size()) < 0) {
                std::cout << "Failed to write to main from worker w="<<w<<" error:"<<std::strerror(errno)<<"\n";
            }
        }
        else {
            // phase 2: send the result directly
            std::string toSend = "result "+std::to_string(graphsCount)
                +" "+doubleToString(br)
                +" "+std::to_string(timeSpendThisGraph)+"\n";
            if(write(fdToMain, toSend.c_str(), toSend.size()) < 0) {
                std::cout << "Failed to write to main from worker w="<<w<<" error:"<<std::strerror(errno)<<"\n";
            }
        }
    }
 
    // clean-up
    std::cout << "worker id="<<w<<" stopping after "<<(TimeDifference(timeStartWorker) / 1000000.0)<<"s.\n";
    std::cout.flush();
    fclose(stream);
    close(fdToMain);
    close(fdReady);
    exit(0);
}
 
int main(int argc, char* argv[]) {
    // print starting time
    {
        time_t timestamp = time(NULL);
        std::cout << "\n\nbranchingAutomated starting at: "<<ctime(&timestamp)<<"\n";
    }
 
    /*
    Number of non-isomorphic connected graphs (counted by nauty-geng):
    n=4:              4
    n=5:             21
    n=6:            112
    n=7:            853
    n=8:         11 117
    n=9:        261 080
    n=10:    11 716 571
     */
 
    int workers = 1;
    long skipUntilNumber = 0;
    AlgorithmTypes algorithmType = GRAMM_ET_AL;
 
    // parse options
    for(int i=1; i<argc; ++i) {
        std::string option = argv[i];
 
        // p: number of processes
        if(option == "-p" && i+1 < argc) {
            workers = std::stoi(argv[++i]);
        }
        // s: $s$-Overlapping Cluster Editing
        else if(option == "-s" && i+1 < argc) {
            s = std::stoi(argv[++i]);
        }
        // k: skip until this graph
        else if(option == "-k" && i+1 < argc) {
            skipUntilNumber = std::stol(argv[++i]);
        }
        // o: threshhold for doing the optimum algorithm
        else if(option == "-o" && i+1 < argc) {
            branchingNumberTreshold = std::stod(argv[++i]);
        }
        // d: depth cutoff
        else if(option == "-d" && i+1 < argc) {
            depth_cutoff = std::stoi(argv[++i]);
        }
        // u: use proposition algorithm
        else if(option == "-u") {
            std::cout << "Using proposition algorithm FORBIDDEN_SUBGRAPH_SINGLE_GREEDY_ALL_PAIRS.\n";
            algorithmType = FORBIDDEN_SUBGRAPH_SINGLE_GREEDY_ALL_PAIRS;
        }
        // d: print debug information
        else if(option == "--debug") {
            std::cout << "Printing debug information. Recommended to be used only with one single input graph, e.g., echo \"ECZ?\" | ./out/branchingAutomated -p 1 -s 1 -u \n\n";
            debug = true;
        }
        // t: run tests instead
        else if(option == "-t") {
            std::cout << "Starting tests:\n";
            auto failures = test();
            std::cout << "Finished tests with "<<failures<<" failures\n";
            return failures ? 1 : 0;
        }
    }
    if(workers < 1) workers = 1;
 
    // workers write here if they are ready to receive their next graph (phase 1)
    Pipe pipeReadyWorkers;
    // main writes graphs here -> workers
    std::vector<Pipe> pipesWriteFromMain(workers);
    // workers write their saves (after SYNC) and phase-2 results here -> main
    std::vector<Pipe> pipesWriteFromWorker(workers);
 
    // init pipes
    if(pipe(pipeReadyWorkers.fd) != 0) {
        std::cout << "Could not create pipe\n";
        exit(1);
    }
    for(int i=0;i<workers;i++) {
        if(pipe(pipesWriteFromMain[i].fd) != 0) {
            std::cout << "Could not create pipe\n";
            exit(1);
        }
        if(pipe(pipesWriteFromWorker[i].fd) != 0) {
            std::cout << "Could not create pipe\n";
            exit(1);
        }
    }
 
    std::cout << "\n";
    std::cout.flush(); // flush before fork, otherwise the buffer is duplicated into every worker
 
    // fork workers
    for(int w=0; w<workers; ++w) {
        pid_t pid = fork();
        if(pid < 0) {
            std::cout << "Could not fork worker "<<w<<" error:"<<std::strerror(errno)<<"\n";
            exit(1);
        }
 
        // child: never returns
        if(pid == 0) {
            workerMain(w, workers, algorithmType, pipeReadyWorkers, pipesWriteFromMain, pipesWriteFromWorker);
        }
 
        close(pipesWriteFromMain[w].fd[0]); // parent closes read end of the write pipe
        close(pipesWriteFromWorker[w].fd[1]); // parent closes write end of the read pipe
    }
    close(pipeReadyWorkers.fd[1]); // parent closes write end of the ready pipe
 
    // streams to read from the workers, used in both phases
    std::vector<FILE*> workerStreams(workers);
    for(int w=0; w<workers; ++w) {
        workerStreams[w] = fdopen(pipesWriteFromWorker[w].fd[0], "r");
    }
 
    auto timeStart = TimeNow();
    std::string line;
    long graphsCount = 0;
    long graphsSent = 0;
 
    // skip first graphs
    if(skipUntilNumber > 0) {
        std::cout << "SKIPPING GRAPHS UNTIL graphCount="<<skipUntilNumber<<"\n";
 
        while (std::getline(std::cin, line)) {
            if(line.empty()) continue;
 
            ++graphsCount;
            if(graphsCount+1 >= skipUntilNumber) break;
        }
    }
 
    //////////////////////////////////////////////////////////////
    // PHASE 1: distribute all input graphs on demand

    // function to read stdin and send a worker the next graph
    // returns TRUE if a graph was in stdin, and FALSE if no next graph exists
    auto sendNextGraph = [&](int w) -> bool {
        std::string line;
        while (std::getline(std::cin, line)) {
            if (line.empty()) continue;
            ++graphsCount;
            ++graphsSent;
            std::string toSend = std::to_string(graphsCount) + " " + line + "\n";
            if(write(pipesWriteFromMain[w].fd[1], toSend.c_str(), toSend.size()) < 0) {
                std::cout << "Failed to write to worker "<<w<<" error:"<<std::strerror(errno)<<"\n";
            }
            return true;
        }
        return false;
    };
 
    // number of workers that received "SYNC"
    int synchronisedWorkers = 0;

    // function to send the SYNC message to workers
    auto sendSync = [&pipesWriteFromMain, &synchronisedWorkers](int w) {
        std::string toSend = "SYNC\n";
        if(write(pipesWriteFromMain[w].fd[1], toSend.c_str(), toSend.size()) < 0) {
            std::cout << "Failed to write to worker "<<w<<" error:"<<std::strerror(errno)<<"\n";
        }
        ++synchronisedWorkers;
    };
 
    // initial graph for workers (round-robin the first round)
    for(int w=0; w<workers; ++w) {
        if(!sendNextGraph(w)) sendSync(w);
    }
 
    // wait for the next worker to be ready, until every worker received "SYNC"
    // (the workers stay alive, so we cannot wait for EOF on the ready pipe)
    FILE* readyStream = fdopen(pipeReadyWorkers.fd[0], "r");
    char buffer[4096];
    while(synchronisedWorkers < workers && fgets(buffer, sizeof(buffer), readyStream)) {
        std::string msg(buffer);
        if(!msg.empty() && msg.back() == '\n') msg.pop_back();
        if(msg.empty()) continue;
 
        // split at first space
        size_t pos = msg.find(' ');
        if(pos == std::string::npos) continue;
 
        std::string msgType = msg.substr(0, pos);
 
        // worker is ready to receive their next graph
        if(msgType == "ready") {
            int w = std::stoi(msg.substr(pos + 1));
            if(!sendNextGraph(w)) sendSync(w);
        }
        else {
            std::cout << "Unknown worker message type="<<msgType<<"\n";
        }
    }
    if(synchronisedWorkers < workers) {
        std::cout << "Ready pipe closed before all workers were synchronised (did a worker crash?)\n";
        exit(1);
    }
 
    //////////////////////////////////////////////////////////////
    // SYNCHRONISATION: receive the saves (graph + found branching number of greedy) of all workers
 
    long totalRunningTimeOfWorkers = 0L;
    auto& saves = graph_with_branching_number; // empty in main so far
    for(int w=0; w<workers; ++w) {
        bool done = false;
        while(fgets(buffer, sizeof(buffer), workerStreams[w])) {
            std::string msg(buffer);
            if(!msg.empty() && msg.back() == '\n') msg.pop_back();
            std::istringstream iss(msg);
            std::string msgType;
            iss >> msgType;
 
            if(msgType == "save") {
                graphSaveType save;
                int algorithmRun;
                iss >> save.graphsCount >> algorithmRun >> save.branching_number >> save.line;
                save.algorithmRun = (AlgorithmTypes)algorithmRun;
                saves[save.graphsCount] = save;
            }
            else if(msgType == "done") {
                long microseconds;
                iss >> microseconds;
                totalRunningTimeOfWorkers += microseconds;
                done = true;
                break;
            }
            else {
                std::cout << "Unknown message from worker "<<w<<": "<<msg<<"\n";
            }
        }
        if(!done) {
            std::cout << "Worker "<<w<<" closed its pipe before sending all saves (did it crash?)\n";
            exit(1);
        }
    }
 
    long timePhase1 = TimeDifference(timeStart);
    std::cout << "\nPhase 1 finished in "<<(timePhase1 / 1000000.0)<<"s."
        <<" Received "<<saves.size()<<" saves ("<<(graphsSent - (long)saves.size())<<" trivial graphs without forbidden subgraph).\n\n";
    std::cout.flush();
 
    //////////////////////////////////////////////////////////////
    // PHASE 2: run the optimal algorithm on the graphs with the highest branching numbers
    // until the worst-case branching number is certified by the optimal algorithm
 
    auto timeStartPhase2 = TimeNow();
    int rounds = 0;
    long optimalRuns = 0;
    long improvedRuns = 0;
    while(true) {
        // find the largest optimal branching number
        double maxOptimal = 0;
        for(const auto& entry : saves) {
            const auto& save = entry.second;
            if(!isGreedy(save.algorithmRun)) maxOptimal = std::max(maxOptimal, save.branching_number);
        }
 
        // find candidates: greedy results that could still be the worst case
        std::vector<graphSaveType*> candidates;
        for(auto& entry : saves) {
            auto& save = entry.second;
            if(isGreedy(save.algorithmRun) && save.branching_number > maxOptimal) {
                candidates.push_back(&save);
            }
        }
 
        // the worst case is an optimal result: done
        if(candidates.empty()) break;
 
        // the (up to) |workers| candidates with the highest branching number (ties: smallest graph number)
        size_t batch = std::min((size_t)workers, candidates.size());
        std::partial_sort(candidates.begin(), candidates.begin() + batch, candidates.end(),
            [](const graphSaveType* a, const graphSaveType* b) {
                if(a->branching_number != b->branching_number) return a->branching_number > b->branching_number;
                return a->graphsCount < b->graphsCount;
            });
 
        ++rounds;
        std::cout << "Phase 2 round "<<rounds<<": "<<candidates.size()<<" candidates above the largest optimal br="<<maxOptimal
            <<", running the optimal algorithm on "<<batch<<" graphs\n";
        std::cout.flush();
 
        // send candidates to workers
        for(size_t i=0; i<batch; ++i) {
            std::string toSend = std::to_string(candidates[i]->graphsCount) + " " + candidates[i]->line + "\n";
            if(write(pipesWriteFromMain[i].fd[1], toSend.c_str(), toSend.size()) < 0) {
                std::cout << "Failed to write to worker "<<i<<" error:"<<std::strerror(errno)<<"\n";
            }
        }
 
        // collect the results (one per worker that received a graph)
        for(size_t i=0; i<batch; ++i) {
            if(!fgets(buffer, sizeof(buffer), workerStreams[i])) {
                std::cout << __FILE__ <<":" << __LINE__ << " could not send to worker "<<i<<"\n";
                exit(1);
            }
            std::string msg(buffer);
            if(!msg.empty() && msg.back() == '\n') msg.pop_back();
 
            std::istringstream iss(msg);
            std::string msgType;
            int resultGraphsCount;
            double br;
            long microseconds;
            iss >> msgType >> resultGraphsCount >> br >> microseconds;
            if(msgType != "result") {
                std::cout << "Unknown message from worker "<<i<<": "<<msg<<"\n";
                exit(1);
            }
 
            auto& save = saves.at(resultGraphsCount);
            std::cout << "\tgraph "<<resultGraphsCount<<": "<<save.line
                <<" br="<<save.branching_number<<" ("<<algorithmTypeToString(save.algorithmRun)<<")"
                <<" -> optimal br="<<br
                <<" time="<<(microseconds / 1000000.0)<<"s\n";
 
            save.algorithmRun = GRAMM_ET_AL;
            totalRunningTimeOfWorkers += microseconds;

            if(save.branching_number-br > 0.000001) {
                save.branching_number = br;
                ++improvedRuns;
            }
            ++optimalRuns;
        }
    }
    long timePhase2 = TimeDifference(timeStartPhase2);
 
    // stop workers: they exit when their pipe is closed
    for(int w=0; w<workers; ++w) {
        close(pipesWriteFromMain[w].fd[1]);
    }
    for(int w=0; w<workers; ++w) {
        wait(nullptr);
    }
    for(int w=0; w<workers; ++w) {
        fclose(workerStreams[w]);
    }
    fclose(readyStream);
 
    //////////////////////////////////////////////////////////////
    // final report
 
    long timeSpend = TimeDifference(timeStart);
 
    const graphSaveType* worst = nullptr;
    double sumBranchingNumbers = 0;
    std::unordered_map<int, long> countPerAlgorithm;
    for(const auto& entry : saves) {
        const auto& save = entry.second;
        sumBranchingNumbers += save.branching_number;
        ++countPerAlgorithm[(int)save.algorithmRun];
        // ties: prefer the non-greedy result, since it certifies the worst case
        if(worst == nullptr || save.branching_number > worst->branching_number
            || (save.branching_number == worst->branching_number && isGreedy(worst->algorithmRun) && !isGreedy(save.algorithmRun))) {
            worst = &save;
        }
    }
 
    std::cout << "\nbranchingAutomated finished."
        <<"\n\tReal time: "<<(timeSpend / 1000000.0)<<"s"
        <<" (phase 1: "<<(timePhase1 / 1000000.0)<<"s, phase 2: "<<(timePhase2 / 1000000.0)<<"s,"
        <<" "<<rounds<<" rounds, "<<optimalRuns<<" optimal runs, "<<improvedRuns<<" optimal IMPROVED runs)"
        <<"\n\tTotal aggregated worker time: "<<(totalRunningTimeOfWorkers / 1000000.0)<<"s"
        <<"\n\tGraphs: "<<graphsSent<<" ("<<saves.size()<<" non-trivial)";
    for(const auto& entry : countPerAlgorithm) {
        std::cout << "\n\t\t"<<algorithmTypeToString((AlgorithmTypes)entry.first)<<": "<<entry.second;
    }
    if(worst != nullptr) {
        std::cout
            <<"\n\tWorst branching number="<<doubleToString(worst->branching_number)
            <<" (graph "<<worst->graphsCount<<": "<<worst->line<<", "<<algorithmTypeToString(worst->algorithmRun)<<")"
            <<"\n\tAverage branching number="<<doubleToString(sumBranchingNumbers / saves.size())
            <<" (over the "<<saves.size()<<" non-trivial graphs)";
    }
    else {
        std::cout << "\n\tNo non-trivial graphs.";
    }
    std::cout << "\n\n";
    return 0;
}







