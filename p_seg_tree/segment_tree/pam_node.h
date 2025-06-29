#ifndef PAM_NODE_H
#define PAM_NODE_H

#include <pam/pam.h>
#include "edge.h"

// Define the comparison struct needed for the pam_set
struct edge_compare {
    bool operator()(const Edge& a, const Edge& b) const {
        return a.y < b.y;
    }
};

// Define the Entry type for our PAM set.
// It will be an ordered set of Edges.
struct PamEdgeSet {
    using key_t = Edge;
    using compare = edge_compare;
};

using pam_edge_set = pam_set<PamEdgeSet>;

// The new Node class, now using a PAM set for its lists
class Node {
public:
    Edge interval;
    
    // The coverList and endList are now both a pam_edge_set.
    // This single data structure is thread-safe and ordered.
    pam_edge_set coverList;
    pam_edge_set endList;

    Node() {}

    // No need for addToCoverList or addToEndList methods,
    // we will call pam_set::insert directly.
};

#endif // PAM_NODE_H
