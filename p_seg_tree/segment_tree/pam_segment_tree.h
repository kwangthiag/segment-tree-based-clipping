#ifndef PAM_SEGTREE_H
#define PAM_SEGTREE_H

#include "edge.h"
#include "../util/util.h"
#include <pam/pam.h>
#include <parlay/sequence.h>
#include <vector>
#include <algorithm>
#include <omp.h>

// Define the Entry struct required by PAM.
struct PamEdgeSet {
    using key_t = Edge;

    // The comparison function MUST be a static member named 'comp'.
    // CORRECTED: Compare based on 'start' and 'end' members, which exist in Edge.
    static bool comp(const Edge& a, const Edge& b) {
        if (a.start != b.start) return a.start < b.start;
        if (a.end != b.end) return a.end < b.end;
        return a.id < b.id; // Use ID as a final tie-breaker
    }
};

// Create a clear type alias for our new parallel-safe set.
using pam_edge_set = pam_set<PamEdgeSet>;


// The SegTree class definition.
class SegTree {
private:
    // These are now proper C++ vector objects, not pointers.
    std::vector<pam_edge_set> coverList;
    std::vector<pam_edge_set> endList;

    // Helper function to find the leaf node for a given x-coordinate
    int find_leaf_for_point(REAL p);
    
    // CORRECTED: Pass Edge by value to avoid const-correctness issues in this recursive function.
    void insertEdgeParallel(const Edge& x, int nodeIdx);

public:
    // Public variables needed to maintain compatibility with clip.cpp
    int intervalCount;
    int treeSize;
    int numLeaves;
    int POLY_TYPE_COUNT;
    
    // These arrays are still needed for the original build logic
    REAL *elementaryEdgeArray;
    REAL *nodeIntervalArray;

    // Constructor & Destructor
    SegTree(int aNumEdges, int ic, REAL minVal, REAL maxVal);
    ~SegTree();

    // Public function declarations (the API)
    void init();
    void constructElementaryIntervalsArrayMulticore(Edge *intervals);
    void constructTreeIntervalArrayMulticore();
    
    void insertAllParallel(Edge *edges);
    void insertAllToEdgeListMulticore(Edge *edges, int bSize);

    void saveIntersectionsIdsMulticore(
        REAL *bPoly, REAL *cPoly, int bSize, int cSize,
        int bType, int cType,
        std::vector<int> &bPolLineIds, std::vector<int> &cPolLineIds, std::vector<int> &intersectTypesDup);
};

#endif // PAM_SEGTREE_H
