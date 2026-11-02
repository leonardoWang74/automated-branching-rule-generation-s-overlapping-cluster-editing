#ifndef GRAPH_H
#define GRAPH_H

// if defined, then instead of std::vector<std::unordered_set<int>> use
// 1. a matrix std::vector<std::vector<bool>> and
// 2. a std::vector<std::vector<int>> as sorted adjacency list
#define GRAPH_H_MATRIX_AND_LIST

#include <functional>
#include <vector>
#include <string>
#include <optional>
#include <unordered_set>
#include <chrono>

std::chrono::_V2::system_clock::time_point TimeNow();
long TimeDifference(const std::chrono::_V2::system_clock::time_point& start);

struct MaximalCliquesInfo {
    // if FALSE do not push to cliques
    bool cliqueListEnabled = true;

    // list of cliques
#ifndef GRAPH_H_MATRIX_AND_LIST
    std::vector<std::unordered_set<int>>
#else
    std::vector<std::vector<int>> 
#endif
    cliqueList;

    std::vector<std::unordered_set<int>> cliqueSets = {};

    // map of vertex -> overlapping clique indices
    std::vector<std::vector<size_t>> vertexCliques;
    // vertex in at least s+1 cliques. If no such vertex exists, set to -1
    int vertexInMoreThanSCliques = -1;
};

struct DegeneracyAndOrdering {
    // the degeneracy of the graph
    int degeneracy;
    // the degeneracy ordering of the graph
    std::vector<int> ordering;

    // the minimum number of edits needed to solve this instance of $s$-Overlapping Cluster Editing
    int editBound;
};

class Graph {
  public:
    MaximalCliquesInfo getMaximalCliques(size_t s=0, bool stopAfterOneVertexInMoreThanS=false);
    DegeneracyAndOrdering getDegeneracyOrdering(int s=0, int k=0) const;

    std::vector<Graph> getComponents() const;
    Graph getSubgraph(const std::vector<int>& vertex_ids) const;

    bool edge_has(int v, int w) const;
    void edge_add(int v, int w);
    void edge_remove(int v, int w);

    int degree(int v) const;

    int id_get(const int v) const;
    int id_reverse_get(const int v) const;
    bool id_has() const;

    unsigned int n() const;
    int n_signed() const;
    unsigned int m() const;

    static Graph parse_graph6(const std::string& g6);
    std::string to_graph6() const;

    explicit Graph(int n);
    explicit Graph(const Graph* G);

    using iterator = std::vector<int>::iterator;
    using const_iterator = std::vector<int>::const_iterator;
    
    template <typename T> static std::string vector_tostring(const T vec);
    template <typename T> static std::string vector_tostring(const std::unordered_set<T>& vec);
    template <typename T> static std::string vector_tostring(const std::vector<T>& vec);
    static std::string stringvector_tostring(const std::vector<std::string>& vec);
    // template <typename T> static std::string vector_tostring(const std::vector<std::vector<T>>& vec);

    template <typename T> static std::vector<T> vector_slice(const std::vector<T>& vec, size_t from, size_t to);

    static bool sorted_contains(const std::vector<int>& vec, int x);
    static void sorted_insert(std::vector<int>& vec, int x);
    static void sorted_remove(std::vector<int>& vec, int x);

    static std::unordered_set<int> set_union(const std::unordered_set<int>& a, const std::unordered_set<int>& b);
    static std::unordered_set<int> set_intersection(const std::unordered_set<int>& a, const std::unordered_set<int>& b);
    static std::unordered_set<int> set_difference(const std::unordered_set<int>& a, const std::unordered_set<int>& b);

    static std::vector<int> sorted_union_unique(const std::vector<int>& a, const std::vector<int>& b);
    static std::vector<int> sorted_intersection_unique(const std::vector<int>& a, const std::vector<int>& b);
    static std::vector<int> sorted_difference(const std::vector<int>& a, const std::vector<int>& b);

    static std::vector<int> sorted_union_unique_slice(const std::vector<int>& a, size_t aFrom, size_t aTo, const std::vector<int>& b, size_t bFrom, size_t bTo);
    static std::vector<int> sorted_intersection_unique_slice(const std::vector<int>& a, size_t aFrom, size_t aTo, const std::vector<int>& b, size_t bFrom, size_t bTo);
    static std::vector<int> sorted_difference_slice(const std::vector<int>& a, size_t aFrom, size_t aTo, const std::vector<int>& b, size_t bFrom, size_t bTo);

#ifndef GRAPH_H_MATRIX_AND_LIST
    // edges of vertices as adjacency lists: O(deg(v)) enumeration of neighbors.
    // O(1) for checking if an edge exists.
    // O(1) for inserting/removing an edge
    // edges[vertex_id] is the set of adjacent vertices of the vertex. vertex_id \in [0, n)
    std::vector<std::unordered_set<int>> edges_list;
    const std::unordered_set<int>& neighbors(int v) const;
#else
    const std::vector<int>& neighbors(int v) const;

    // edges of vertices as adjacency lists: O(deg(v)) enumeration of neighbors.
    // O(log(deg(v))) for checking if an edge exists.
    // O(deg(v)) for inserting/removing an edge (shifting elements)
    // edges[vertex_id] is the set of adjacent vertices of the vertex. vertex_id \in [0, n)
    std::vector<std::vector<int>> edges;

    // edges of vertices as matrix: O(1) checking if an edge exists.
    // O(1) for inserting/removing an edge. O(n) for enumerating neighbors.
    // edges[vertex_id][vertex_id2] is TRUE if the edge exists, otherwise FALSE
    std::vector<std::vector<bool>> edges_matrix;
#endif

    // a map of vertex_id [0,H->n()] -> mapped vertex_id [0,G->n()] with H as this = subgraph
    std::vector<int> ids;
    // a map of vertex_id [0,G->n()] -> mapped vertex_id [0,H->n()] with H as this = subgraph
    std::vector<int> ids_reverse;
  private:
    unsigned int number_vertices;
    unsigned int number_edges;

    // TRUE if ids has been initialized
    bool ids_initialized;
};

template <typename T> std::vector<T> Graph::vector_slice(const std::vector<T>& vec, size_t from, size_t to) {
    const size_t n = to - from;
    std::vector<T> result = std::vector<T>(n);
    for(size_t i=0; i<n; ++i) {
        result[i] = vec[from + i];
    }
    return result;
}

template <typename T> std::string Graph::vector_tostring(const T value) {
    return std::to_string(value);
}

template <typename T> std::string Graph::vector_tostring(const std::unordered_set<T>& vec) {
    const auto n = vec.size();
    if(n == 0) return "[]";
    auto it = vec.begin();
    std::string s = "[" + Graph::vector_tostring(*it);
    for(size_t i=1; i<n; ++i) {
        ++it;
        s += "," + Graph::vector_tostring(*it);
    }
    return s + "]";
}

template <typename T> std::string Graph::vector_tostring(const std::vector<T>& vec) {
    const auto n = vec.size();
    if(n == 0) return "[]";
    auto it = vec.begin();
    std::string s = "[" + Graph::vector_tostring(*it);
    for(size_t i=1; i<n; ++i) {
        ++it;
        s += "," + Graph::vector_tostring(*it);
    }
    return s + "]";
}

void SubsetsOfSizeLoop(size_t n, size_t indicesSize, std::function<bool(size_t n, std::vector<size_t>)> function);

#endif