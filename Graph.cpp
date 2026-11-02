#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <chrono>
#include <unordered_set>

#include <cerrno>
#include <cstring>
#include <stdexcept>
#include <iostream>
#include <functional>

#include "Graph.h"

#ifndef DEBUG
// use DEBUG if you want to debug
// #define DEBUG
#endif

// time tracking functions. Usage:
// auto start = TimeNow();
// timeSum += TimeDifference(start);
std::chrono::_V2::system_clock::time_point TimeNow() {
    return std::chrono::high_resolution_clock::now();
}

long TimeDifference(const std::chrono::_V2::system_clock::time_point& start) {
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start);
    return std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count();
}

/** calculate a degeneracy + degeneracy ordering of the graph O(4 * n * \Delta). If s>0, then the bound is also returned
 * [Eppstein et al. 2010 - Listing All Maximal Cliques in Sparse Graphs in Near-optimal Time, Section 2.1 before Lemma 1]
 */
DegeneracyAndOrdering Graph::getDegeneracyOrdering(int s, int k) const {
    int n = this->n();

    // Track smallest degree
    int smallestDegree = n;

    // vertexId -> current degree O(n)
    std::vector<int> verticesMap = std::vector<int>(n);

    // degree -> list of vertex pointers O(n)
    std::vector<std::unordered_set<int>> degreesMap = std::vector<std::unordered_set<int>>(n, std::unordered_set<int>(n));

    // build initial maps in O(n)
    for(int vid=0; vid<n; ++vid) {
        int degree = this->degree(vid);
        verticesMap[vid] = degree;

        smallestDegree = std::min(smallestDegree, degree);

        degreesMap[degree].insert(vid);
    }

    std::vector<int> degeneracyOrdering;
    int degeneracy = 0;
    int editBound = 0;

    // do the ordering O(n * \Delta)
    for(int i=0; i<n; ++i) {
        auto& listSmallest = degreesMap[smallestDegree];

        // get a vertex with the smallest degree
        /*int vid = listSmallest.back();
        listSmallest.pop_back();*/
        int vid = *listSmallest.begin();
        listSmallest.erase(vid);

        // removed last entry in the list - find next smallestDegree
        if(listSmallest.empty() && i<n-1) {
            // find next smallest degree O(\Delta)
            bool found = false;
            for(int d=smallestDegree; d<n; ++d) {
                if(!degreesMap.at(d).empty()) {
                    smallestDegree = d;
                    found = true;
                    break;
                }
            }
            if(!found) {
                std::cout << __FILE__ <<":"<<__LINE__
                    <<" no other smallest degree found for i=" << i << ", smallestDegree= " << smallestDegree << "\n";
                break;
            }
        }

        // "remove" from the graph
        verticesMap.at(vid) = -1;

        // decrease degree of neighbors O(\Delta)
        int degreeHere = 0;
        for(int neighborId : this->neighbors(vid)) {
            int degreePrevious = verticesMap.at(neighborId);
            if(degreePrevious < 0) continue;

            int degreeNew = degreePrevious - 1;
            ++degreeHere;

            // update neighbor degree
            verticesMap.at(neighborId) = degreeNew;

            // remove from old degree list O(1)
            degreesMap.at(degreePrevious).erase(neighborId);

            // add to new degree list O(1)
            degreesMap.at(degreeNew).insert(neighborId);

            // update smallest degree
            if(degreeNew < smallestDegree) {
                smallestDegree = degreeNew;
            }
        }

        degeneracyOrdering.push_back(vid);
        degeneracy = std::max(degeneracy, degreeHere);

        // find a minimum number of edits using the conjecture that stars are bad for all s and that you need $|leaves| - s$ edits
        // in time O(degeneracy * min(k, n/3))
        if (s > 0) {
            const auto nHere = n - i;
            if(nHere < 6) continue;
            for(auto t = std::max(3, s+1); t <= degreeHere; ++t) { // O(degeneracy of G * min(k,n/3))
                int rBound = std::min(k+1+s, 2 + (n-2) / t); // 2 + (n-2)/t since that is 1 + (n-2)/t rounded up
                for(auto r = t+1; r<rBound; ++r) {
                    int l = 1 + (nHere - 1) / (t * (r-1) + 1);
                    std::cout << __FILE__<<":"<<__LINE__<<" Trying t="<<t<<", r="<<r<<", l="<<l<<" rBound="<<rBound<<" with r-s="<<(r-s)<<", l*(t-s)="<<(l*(t-s))<<"\n";
                    editBound = std::max(editBound, std::min(r-s, l * (t-s)));
                }
            }
        }
    }

    return {degeneracy: degeneracy, ordering: degeneracyOrdering, editBound: editBound};
}

/** Bron-Kerbosch with Pivot.
 * Returns -1 if there was no vertex found in more than `s` cliques, otherwise the vertex id (if s=0 then always returns -1).
 * [Eppstein et al. 2010 - Listing All Maximal Cliques in Sparse Graphs in Near-optimal Time, Figure 2: BronKerboschPivot]
 */
int BronKerboschPivot(Graph* G, MaximalCliquesInfo& result, size_t s, bool stopAfterOneVertexInMoreThanS
#ifndef GRAPH_H_MATRIX_AND_LIST
    , std::unordered_set<int>& P, std::unordered_set<int>& R, std::unordered_set<int>& X
#else
    , std::vector<int>& P, std::vector<int>& R, std::vector<int>& X
#endif
) {
    // if $P \cup X = \emptyset$
    if(P.empty() && X.empty()) {
        // report R as a maximal clique
#ifndef GRAPH_H_MATRIX_AND_LIST
        if(result.cliqueListEnabled) result.cliqueList.push_back(std::unordered_set<int>(R));
#else
        if(result.cliqueListEnabled) result.cliqueList.push_back(std::vector<int>(R));
#endif

        // count number of cliques per vertex, if >s then return that vertex id
        if(s > 0) {
            // already found one
            if(stopAfterOneVertexInMoreThanS && result.vertexInMoreThanSCliques != -1) return result.vertexInMoreThanSCliques;

            // count
            for(auto vid : R) {
                auto& list = result.vertexCliques.at(vid);
                list.push_back(result.cliqueList.size() - 1);
                if(list.size() > s) {
                    result.vertexInMoreThanSCliques = vid;
                    if(stopAfterOneVertexInMoreThanS) return vid;
                }
            }
        }
        return -1;
    }

    // choose a pivot $u \in P \cup X$. Tomita, Tanaka, Takahashi 2006: choose $u$ to maximize $|P \cap N(u)|$ with $N(u)$ being the neighborhood of $u$
    // [Tomita, Tanaka, Takahashi 2006 - The worst-case time complexity for generating all maximal cliques and computational experiments]
    // const auto pivotCandidates = Graph::set_union(P, X);
    int pivot = !P.empty() ? *P.begin() : *X.begin();
    size_t pivotValue = 0;
    for(int vid : P) {
        const auto neighbors = G->neighbors(vid);

        // intersection cannot be larger than |neighbors|
        if(neighbors.size() < pivotValue) continue;

#ifndef GRAPH_H_MATRIX_AND_LIST
        auto value = Graph::set_intersection(neighbors, P).size();
#else
        auto value = Graph::sorted_intersection_unique(neighbors, P).size();
#endif
        if(value > pivotValue) {
            pivot = vid;
            pivotValue = value;
        }
    }
    for(int vid : X) {
        const auto neighbors = G->neighbors(vid);

        // intersection cannot be larger than |neighbors|
        if(neighbors.size() < pivotValue) continue;

#ifndef GRAPH_H_MATRIX_AND_LIST
        auto value = Graph::set_intersection(neighbors, P).size();
#else
        auto value = Graph::sorted_intersection_unique(neighbors, P).size();
#endif
        if(value > pivotValue) {
            pivot = vid;
            pivotValue = value;
        }
    }
    if(pivot<0) {
        std::cout << __FILE__<<":"<<__LINE__<< " pivot empty\n";
        return -1;
    }

    // for each vertex $v \in P \setminus N(u)$ do
#ifndef GRAPH_H_MATRIX_AND_LIST
    const auto loopSet = Graph::set_difference(P, G->neighbors(pivot));
#else
    const auto loopSet = Graph::sorted_difference(P, G->neighbors(pivot));
#endif
    for(int vid : loopSet) {
        const auto neighbors = G->neighbors(vid);

#ifndef GRAPH_H_MATRIX_AND_LIST
        auto P_new = Graph::set_intersection(P, neighbors);
        auto X_new = Graph::set_intersection(X, neighbors);
        R.insert(vid);
#else
        auto P_new = Graph::sorted_intersection_unique(P, neighbors);
        auto X_new = Graph::sorted_intersection_unique(X, neighbors);
        R.push_back(vid);
#endif

        // BronKerboschPivot(P \cap N(v), R \cup \{v\}, X \cap N(v))
        const auto vertex = BronKerboschPivot(G, result, s, stopAfterOneVertexInMoreThanS,
            P_new, R, X_new
        );
        if(vertex >= 0) {
            result.vertexInMoreThanSCliques = vertex;
            if(stopAfterOneVertexInMoreThanS) return vertex;
        }

#ifndef GRAPH_H_MATRIX_AND_LIST
        // do not have to copy R
        R.erase(vid);

        // P \leftarrow P \setminus \{v\}
        P.erase(vid);
        
        // X \leftarrow X \cup \{v\}
        X.insert(vid);
#else
        // do not have to copy R
        R.pop_back();

        // P \leftarrow P \setminus \{v\}
        Graph::sorted_remove(P, vid);
        
        // X \leftarrow X \cup \{v\}
        Graph::sorted_insert(X, vid);
#endif
    }
    return -1;
}

/** get maximal cliques using Bron-Kerbosch based on degeneracy by Eppstein, Loeffler and Strash.
 * Returns -1 if there was no vertex found in more than `s` cliques, otherwise the vertex id (if s=0 then always returns -1).
 * [Eppstein et al. 2010 - Listing All Maximal Cliques in Sparse Graphs in Near-optimal Time, Figure 4: BronKerboschDegeneracy]
 */
int BronKerboschDegeneracyByEppsteinLoefflerStrash(Graph* G, MaximalCliquesInfo& result, size_t s, bool stopAfterOneVertexInMoreThanS) {
    auto degeneracyInfo = G->getDegeneracyOrdering();

    // for each vertex vi in a degeneracy ordering $v_0, v_1, v_2, \dots$ of $(V,E)$ do
    for(size_t i = 0; i < degeneracyInfo.ordering.size(); ++i) {
        const int vid = degeneracyInfo.ordering[i];
        const auto neighbors = G->neighbors(vid);
    
        auto next = Graph::vector_slice(degeneracyInfo.ordering, i+1, degeneracyInfo.ordering.size());
#ifndef GRAPH_H_MATRIX_AND_LIST
        auto nextSet = std::unordered_set<int>(next.begin(), next.end());

        // $ P \leftarrow N(v_i) \cap \{v_{i+1}, \dots, v_{n-1}\}$
        // P = Neighborhood intersected with neighbors later in the ordering
        auto P = Graph::set_intersection(neighbors, nextSet);

        // $ X \leftarrow N(v_i) \cap \{v_0, \dots, v_{i-1}\}$
        // X = Neighborhood intersected with neighbors earlier in the ordering

        // $P \cup X = N(v_i)$ holds therefore also $X = N(v_i) \setminus P$
        // so we don't need to calculate the slice of the previous vertices in the ordering
        auto X = Graph::set_difference(neighbors, P);

        // R = \{vid\}
        std::unordered_set<int> R = std::unordered_set<int>();
        R.insert(vid);
#else
        std::sort(next.begin(), next.end());

        // $ P \leftarrow N(v_i) \cap \{v_{i+1}, \dots, v_{n-1}\}$
        // P = Neighborhood intersected with neighbors later in the ordering
        auto P = Graph::sorted_intersection_unique(neighbors, next);

        // $ X \leftarrow N(v_i) \cap \{v_0, \dots, v_{i-1}\}$
        // X = Neighborhood intersected with neighbors earlier in the ordering

        // $P \cup X = N(v_i)$ holds therefore also $X = N(v_i) \setminus P$
        // so we don't need to calculate the slice of the previous vertices in the ordering
        auto X = Graph::sorted_difference(neighbors, P);

        // R = \{vid\}
        std::vector<int> R = std::vector<int>();
        R.push_back(vid);
#endif

        const auto vertex = BronKerboschPivot(G, result, s, stopAfterOneVertexInMoreThanS, P, R, X);
        if(vertex >= 0) {
            result.vertexInMoreThanSCliques = vertex;
            if(stopAfterOneVertexInMoreThanS) return vertex;
        }
    }
    return -1;
}

/** find maximal cliques using Bron-Kerbosch based on degeneracy by Eppstein, Loeffler and Strash.
 * s>0: count number of cliques per vertex, if we found a vertex in >s maximal cliques and stopAfterOneVertexInMoreThanS=true, then stop
 * [Eppstein et al. 2010 - Listing All Maximal Cliques in Sparse Graphs in Near-optimal Time, Figure 4: BronKerboschDegeneracy]
 */
MaximalCliquesInfo Graph::getMaximalCliques(size_t s, bool stopAfterOneVertexInMoreThanS) {
    MaximalCliquesInfo info = MaximalCliquesInfo();
#ifndef GRAPH_H_MATRIX_AND_LIST
    info.cliqueList = std::vector<std::unordered_set<int>>();
#else
    info.cliqueList = std::vector<std::vector<int>>();
#endif
    info.cliqueList.reserve(this->number_vertices / 3);
    info.vertexCliques = std::vector<std::vector<size_t>>(this->n(), std::vector<size_t>());
    info.vertexInMoreThanSCliques = -1;

    BronKerboschDegeneracyByEppsteinLoefflerStrash(this, info, s, stopAfterOneVertexInMoreThanS);
    return info;
}

// get all connected components of the graph
std::vector<Graph> Graph::getComponents() const {
    // result components
    std::vector<Graph> components = std::vector<Graph>();

    // vector of vertex_id -> boolean: TRUE if we already found this vertex
    std::vector<bool> found = std::vector<bool>(this->number_vertices);
    
    for(unsigned int i=0; i<this->number_vertices; ++i) {
        if(found[i]) continue;

        // the subgraph vertex ids
        std::vector<int> vertex_ids = std::vector<int>();
        // the queued vertex ids where we still have to add the neighbors
        std::vector<int> notYetAdded = std::vector<int>();

        // add one vertex and its neighborhood
        found.push_back(i);
        notYetAdded.push_back(i);

        // add neighborhoods until nothing was added
        while(!notYetAdded.empty()) {
            const int v = notYetAdded.front();
            notYetAdded.erase(notYetAdded.begin());

            vertex_ids.push_back(v);

            // add neighbors to queue
            for(auto w : this->neighbors(v)) {
                if(found[w]) continue;
                found[w] = true;
                notYetAdded.push_back(w);
            }
        }

        components.push_back(this->getSubgraph(vertex_ids));
    }

    return components;
}

std::vector<std::vector<int>> copyVectorVectorInt(std::vector<std::vector<int>> forbidden) {
    const auto n = forbidden.size();

    std::vector<std::vector<int>> copy = std::vector<std::vector<int>>(n);

    for(size_t i=0; i<n; ++i) {
        copy[i] = std::vector<int>(forbidden[i]);
    }

    return copy;
}

// get a subgraph of the given vertex IDs
Graph Graph::getSubgraph(const std::vector<int>& vertex_ids) const {
    Graph graph(vertex_ids.size());
    
    // graph.ids = std::vector<int>(vertex_ids);
    graph.ids = std::vector<int>(vertex_ids);
    graph.ids_reverse = std::vector<int>(this->n(), -1);
    graph.ids_initialized = true;

    // set the IDs map
    for(unsigned int i=0; i<graph.number_vertices; ++i) {
        // graph.ids[i] = vertex_ids.at(i); // this would be the ids mapping
        graph.ids_reverse[vertex_ids.at(i)] = i; // ids_reverse is the reverse map of ids
    }

    // copy edges
    for(unsigned int i=0; i<graph.number_vertices; ++i) {
        for(unsigned int j=i+1; j<graph.number_vertices; ++j) {
            const int v = vertex_ids.at(i);
            const int w = vertex_ids.at(j);

            if(!this->edge_has(v, w)) continue;
            graph.edge_add(i, j);
        }
    }

    return graph;
}

// checks whether the graph has an edge: O(1)
bool Graph::edge_has(int v, int w) const {
#ifndef GRAPH_H_MATRIX_AND_LIST
    const auto& list = this->edges_list.at(v);
    return list.find(w) != list.end();
#else
    return this->edges_matrix.at(v).at(w);
#endif
}

// insert an edge into the graph: O(1)
void Graph::edge_add(int v, int w) {
#ifndef GRAPH_H_MATRIX_AND_LIST
    edges_list.at(v).insert(w);
    edges_list.at(w).insert(v);
#else
    Graph::sorted_insert(edges.at(v), w);
    Graph::sorted_insert(edges.at(w), v);
    this->edges_matrix.at(v).at(w) = 1;
    this->edges_matrix.at(w).at(v) = 1;
#endif
    ++this->number_edges;
}

// remove an edge from the graph: O(1)
void Graph::edge_remove(int v, int w) {
#ifndef GRAPH_H_MATRIX_AND_LIST
    edges_list.at(v).erase(w);
    edges_list.at(w).erase(v);
#else
    Graph::sorted_remove(edges.at(v), w);
    Graph::sorted_remove(edges.at(w), v);
    this->edges_matrix.at(v).at(w) = 0;
    this->edges_matrix.at(w).at(v) = 0;
#endif
    --this->number_edges;
}

// returns the degree of a vertex
int Graph::degree(int v) const {
#ifndef GRAPH_H_MATRIX_AND_LIST
    return this->edges_list.at(v).size();
#else
    return this->edges.at(v).size();
#endif
}

#ifndef GRAPH_H_MATRIX_AND_LIST
const std::unordered_set<int>& Graph::neighbors(int v) const {
    return edges_list.at(v);
}
#else
const std::vector<int>& Graph::neighbors(int v) const {
    return edges.at(v);
}
#endif

// returns the mapped ID of the given vertex_id. [0,this->n()] -> [0,Parent_Graph->n()]
int Graph::id_get(const int v) const {
    if(!this->ids_initialized) return -1;
    return this->ids.at(v);
}

// returns the mapped ID of the given vertex_id. [0,Parent_Graph->n()) -> [0,this->n()]
// if this is not a subgraph, returns -2
// if the id does not exist (vertex is not in the subgraph), returns -1
int Graph::id_reverse_get(const int v) const {
    if(!this->ids_initialized) return -2;
    return this->ids_reverse.at(v);
}

// returns whether there is a map of this vertices to other vertex IDs
bool Graph::id_has() const {
    return this->ids_initialized;
}

// returns the number of vertices in this graph
unsigned int Graph::n() const {
    return this->number_vertices;
}
// returns the number of vertices in this graph
int Graph::n_signed() const {
    return (int)this->number_vertices;
}

// returns the number of edges in this graph
unsigned int Graph::m() const {
    return this->number_edges;
}

// parse a graph from a graph6 string
Graph Graph::parse_graph6(const std::string& g6) {
    size_t idx = 0;

    // number of vertices = first character
    int n = g6[idx++] - 63;

    Graph G(n);

    // adjacency bits = other characters
    int bit_buffer = 0;
    int bit_count = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < i; ++j) {
            if (bit_count == 0) {
                bit_buffer = g6[idx++] - 63;
                bit_count = 6;
            }
            bit_count--;
            int bit = (bit_buffer >> bit_count) & 1;

            if(bit==1) {
                G.edge_add(i, j);
            }
        }
    }

    return G;
}

// create a graph6 from this graph
std::string Graph::to_graph6() const {
    int n = this->n();
    std::string g6;

    // graph6 only supports n <= 62 in the single-character format
    if (n > 62) {
        throw std::runtime_error("graph6 encoding for n > 62 not implemented");
    }

    // first character: number of vertices
    g6.push_back(static_cast<char>(n + 63));

    int bit_buffer = 0;
    int bit_count = 0;

    // encode upper triangle (i > j), same order as parsing
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < i; ++j) {
            int bit = this->edge_has(i, j) ? 1 : 0;

            bit_buffer = (bit_buffer << 1) | bit;
            bit_count++;

            if (bit_count == 6) {
                g6.push_back(static_cast<char>(bit_buffer + 63));
                bit_buffer = 0;
                bit_count = 0;
            }
        }
    }

    // flush remaining bits (pad with zeros on the right)
    if (bit_count > 0) {
        bit_buffer <<= (6 - bit_count);
        g6.push_back(static_cast<char>(bit_buffer + 63));
    }

    return g6;
}

// constructor initializing adjacency lists
Graph::Graph(int n) {
    this->number_vertices = n;
    this->number_edges = 0;
    this->ids_initialized = false;

    // create empty adjacency lists
#ifndef GRAPH_H_MATRIX_AND_LIST
    this->edges_list.assign(n, std::unordered_set<int>(n));
#else
    this->edges.assign(n, std::vector<int>());
    this->edges_matrix.assign(n, std::vector<bool>(n, false));
#endif
}

// copy constructor
Graph::Graph(const Graph* G) {
    this->number_vertices = G->number_vertices;
    this->number_edges = G->number_edges;

    this->ids_initialized = G->ids_initialized;
    if(G->ids_initialized) {
        this->ids = std::vector<int>(G->ids);
        this->ids_reverse = std::vector<int>(G->ids_reverse);
    }

    // copy adjacency lists and matrix
#ifndef GRAPH_H_MATRIX_AND_LIST
    this->edges_list = std::vector<std::unordered_set<int>>(this->number_vertices);
    for(size_t i=0; i<this->number_vertices; ++i) {
        this->edges_list[i] = std::unordered_set<int>(G->edges_list[i]);
    }
#else
    this->edges = std::vector<std::vector<int>>(this->number_vertices);
    this->edges_matrix = std::vector<std::vector<bool>>(this->number_vertices);
    for(size_t i=0; i<this->number_vertices; ++i) {
        this->edges[i] = std::vector<int>(G->edges[i]);
        this->edges_matrix[i] = std::vector<bool>(G->edges_matrix[i]);
    }
#endif
}

bool Graph::sorted_contains(const std::vector<int>& vec, int x) {
    return std::binary_search(vec.begin(), vec.end(), x);
}

// insert an element into a vector at the correctly sorted position
void Graph::sorted_insert(std::vector<int>& vec, int x) {
    auto it = std::lower_bound(vec.begin(), vec.end(), x);
    if (it == vec.end() || *it != x) {
        vec.insert(it, x);
    }
}

// remove an element from a sorted vector
void Graph::sorted_remove(std::vector<int>& vec, int x) {
    auto it = std::lower_bound(vec.begin(), vec.end(), x);
    if (it != vec.end() && *it == x) {
        vec.erase(it);
    }
}

// return a new set as the set union: elements in either `a` or `b` in time O(|a| + |b|)
std::unordered_set<int> Graph::set_union(const std::unordered_set<int>& a, const std::unordered_set<int>& b) {
    std::unordered_set<int> set = std::unordered_set<int>(a.size() + b.size());

    for(auto el : a) set.insert(el);
    for(auto el : b) set.insert(el);

    return set;
}
// return a new set as the set intersection: elements in both `a` and `b` in time O(min(|a|, |b|))
std::unordered_set<int> Graph::set_intersection(const std::unordered_set<int>& a, const std::unordered_set<int>& b) {
    const auto& smaller = a.size() > b.size() ? a : b;
    const auto& larger = a.size() > b.size() ? b : a;
    
    std::unordered_set<int> set = std::unordered_set<int>(smaller.size());
    for(auto el : smaller) {
        // other does not contain element - skip
        if(larger.find(el) == larger.end()) continue;
        set.insert(el);
    }

    return set;
}
// return a new set as the set difference: all the elements in `a` that are not in `b` in time O(|a|)
std::unordered_set<int> Graph::set_difference(const std::unordered_set<int>& a, const std::unordered_set<int>& b) {
    std::unordered_set<int> set = std::unordered_set<int>(a.size());
    for(auto el : a) {
        // other contains element - skip
        if(b.find(el) != b.end()) continue;
        set.insert(el);
    }
    return set;
}

// merges two sorted vectors: elements in either `a` and `b`
std::vector<int> Graph::sorted_union_unique(const std::vector<int>& a, const std::vector<int>& b) {
    return Graph::sorted_union_unique_slice(a, 0, a.size(), b, 0, b.size());
}
// intersect two sorted vectors: elements in both `a` and `b`. Assumes unique elements.
std::vector<int> Graph::sorted_intersection_unique(const std::vector<int>& a, const std::vector<int>& b) {
    return Graph::sorted_intersection_unique_slice(a, 0, a.size(), b, 0, b.size());
}
// returns a sorted vector with all the elements in `a` that are not in `b`
std::vector<int> Graph::sorted_difference(const std::vector<int>& a, const std::vector<int>& b) {
    return Graph::sorted_difference_slice(a, 0, a.size(), b, 0, b.size());
}

// merges two sorted vectors: elements in either `a` and `b`. Using the index slices [aFrom,aTo) and [bFrom,bTo).
std::vector<int> Graph::sorted_union_unique_slice(const std::vector<int>& a, size_t aFrom, size_t aTo, const std::vector<int>& b, size_t bFrom, size_t bTo) {
    std::vector<int> list = std::vector<int>();

    size_t i=aFrom;
    size_t j=bFrom;

    int elementA = i < aTo ? a[i] : INT32_MAX;
    int elementB = j < bTo ? b[j] : INT32_MAX;

    while(i < aTo || j < bTo) {
        if(elementA == elementB) {
            list.push_back(elementA);
            ++i;
            ++j;
            if(i < aTo) elementA = a[i];
            else elementA = INT32_MAX;

            if(j < bTo) elementB = b[j];
            else elementB = INT32_MAX;
        }
        else if(elementA < elementB) {
            list.push_back(elementA);
            ++i;
            if(i < aTo) elementA = a[i];
            else elementA = INT32_MAX;
        }
        // elementA > elementB
        else {
            list.push_back(elementB);
            ++j;
            if(j < bTo) elementB = b[j];
            else elementB = INT32_MAX;
        }
    }

    return list;
}

// intersect two sorted vectors: elements in both `a` and `b`. Assumes unique elements. Using the index slices [aFrom,aTo) and [bFrom,bTo).
std::vector<int> Graph::sorted_intersection_unique_slice(const std::vector<int>& a, size_t aFrom, size_t aTo, const std::vector<int>& b, size_t bFrom, size_t bTo) {
    std::vector<int> list = std::vector<int>();

    size_t i=aFrom;
    size_t j=bFrom;

    int elementA = i < aTo ? a[i] : INT32_MAX;
    int elementB = j < bTo ? b[j] : INT32_MAX;

    // && here since can only intersect as long as we have elements from both
    while(i < aTo && j < bTo) {
        if(elementA == elementB) {
            list.push_back(elementA);
            ++i;
            ++j;
            if(i < aTo) elementA = a[i];
            else elementA = INT32_MAX;

            if(j < bTo) elementB = b[j];
            else elementB = INT32_MAX;
        }
        else if(elementA < elementB) {
            ++i;
            if(i < aTo) elementA = a[i];
            else elementA = INT32_MAX;
        }
        // elementA > elementB
        else {
            ++j;
            if(j < bTo) elementB = b[j];
            else elementB = INT32_MAX;
        }
    }

    return list;
}

// returns a sorted vector with all the elements in `a` that are not in `b`. Using the index slices [aFrom,aTo) and [bFrom,bTo).
std::vector<int> Graph::sorted_difference_slice(const std::vector<int>& a, size_t aFrom, size_t aTo, const std::vector<int>& b, size_t bFrom, size_t bTo) {
    std::vector<int> list = std::vector<int>();

    size_t i=aFrom;
    size_t j=bFrom;

    int elementA = i < aTo ? a[i] : INT32_MAX;
    int elementB = j < bTo ? b[j] : INT32_MAX;

    while(i < aTo) {
        if(elementA == elementB) {
            ++i;
            ++j;
            if(i < aTo) elementA = a[i];
            else elementA = INT32_MAX;

            if(j < bTo) elementB = b[j];
            else elementB = INT32_MAX;
        }
        else if(elementA < elementB) {
            list.push_back(elementA);
            ++i;
            if(i < aTo) elementA = a[i];
            else elementA = INT32_MAX;
        }
        // elementA > elementB
        else {
            while(elementA > elementB){
                ++j;
                if(j < bTo) elementB = b[j];
                else {
                    elementB = INT32_MAX;
                    break;
                }
            }
        }
    }

    return list;
}

// loop through all subset indices `0 .. n-1` and call `function` with the indices (of size `indicesSize`) as long as `function` returns TRUE.
// Example: with n=4, indicesSize=2 the indices are [0,1], [0,2], [0,3], [1,2], [1,3], [2,3]
void SubsetsOfSizeLoop(size_t n, size_t indicesSize, std::function<bool(size_t n, std::vector<size_t>)> function) {
    std::vector<size_t> indices = std::vector<size_t>(indicesSize);
    size_t lastIndexIndex = indices.size() - 1;

    // initial indices
    for(size_t i=0; i<indices.size(); ++i) {
        indices[i] = i;
    }

    // go through all possible indices
    bool loop = true;
    while(loop) {
        // call function with indices
        auto functionCall = function(n, indices);
        if(!functionCall) break;
        
        // increment indices
        size_t incrementIndex = lastIndexIndex;
        while(true) {
            auto newValue = indices.at(incrementIndex) + 1;

            // first index can only go up so far since the other indices need space
            const auto indexBound = n - lastIndexIndex + incrementIndex;

            // wrap-around: also increment next index
            if(newValue >= indexBound) {
                // first index reached the last value: stop the outer loop
                if(incrementIndex == 0) {
                    loop = false;
                    break;
                }
                // other index reached the last value: set to value at previous index + 2 since previous will also be incremented
                else {
                    indices.at(incrementIndex) = indices.at(incrementIndex-1) + 2;
                }
                --incrementIndex;
            }
            // otherwise: only increment last index
            else {
                indices[incrementIndex] = newValue;
                break;
            }
        }

        // set new subgraph vertex for other incrementIndex and
        // adjust other indices. Example degree=8, increment [0,1,5,6,7] -> [0,2,3,7,8] -> [0,2,3,4,5]
        while(++incrementIndex < indicesSize) {
            indices[incrementIndex] = indices[incrementIndex-1] + 1;
        }
    }
}

std::string Graph::stringvector_tostring(const std::vector<std::string>& vec) {
    const auto n = vec.size();
    if(n == 0) return "[]";
    auto it = vec.begin();
    std::string s = "[" + (*it);
    for(size_t i=1; i<n; ++i) {
        ++it;
        s += "," + (*it);
    }
    return s + "]";
}
