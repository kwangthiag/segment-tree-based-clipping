#include "pam_segment_tree.h" // The new header file using PAM
#include <omp.h>
#include <thrust/scan.h>
#include <parlay/primitives.h>
#include <parlay/sequence.h>
#include <pam/pam.h> 
#include <algorithm>
#include <vector>

// --- Constructor, Destructor, and Initialization ---

SegTree::SegTree(int aNumEdges, int ic, REAL minVal, REAL maxVal) {
    this->intervalCount = ic;
    this->POLY_TYPE_COUNT = 0;
    this->treeSize = 0;
    this->numLeaves = 0;
    this->elementaryEdgeArray = nullptr;
    this->nodeIntervalArray = nullptr;
    // Note: The vector members coverList and endList are default-constructed.
}

SegTree::~SegTree() {
    // CORRECTED: Safely free memory allocated with malloc.
    // The vector members `coverList` and `endList` are destroyed automatically.
    if (elementaryEdgeArray) free(elementaryEdgeArray);
    if (nodeIntervalArray) free(nodeIntervalArray);
}

void SegTree::init() {
    if (treeSize > 0) {
        // CORRECTED: Correctly resize the vectors to hold the PAM sets for each node.
        coverList.resize(treeSize + 1);
        endList.resize(treeSize + 1);
    }
}


// --- Core Build Functions (Augmented with PAM) ---

void SegTree::constructElementaryIntervalsArrayMulticore(Edge* intervals) {
    auto all_endpoints = parlay::sequence<REAL>::from_function(intervalCount * 2, [&](size_t i) {
        return (i % 2 == 0) ? intervals[i / 2].start : intervals[i / 2].end;
    });

    parlay::sequence<REAL> unique_points = parlay::unique(parlay::sort(all_endpoints));
    
    int eleCount = unique_points.size();
    if (eleCount < 2) {
      numLeaves = 0;
      treeSize = 0;
      return;
    }
    // CORRECTED: Added Util namespace qualifier
    numLeaves = Util::getNPowOf2(eleCount - 1);
    treeSize = 2 * numLeaves;

    elementaryEdgeArray = (REAL*)malloc(numLeaves * 2 * sizeof(REAL));
    #pragma omp parallel for
    for (int i = 0; i < eleCount - 1; ++i) {
        elementaryEdgeArray[i * 2] = unique_points[i];
        elementaryEdgeArray[i * 2 + 1] = unique_points[i + 1];
    }
    for (int i = eleCount - 1; i < numLeaves; ++i) {
        elementaryEdgeArray[i * 2] = unique_points.back();
        elementaryEdgeArray[i * 2 + 1] = unique_points.back();
    }
}

void SegTree::constructTreeIntervalArrayMulticore() {
    if (treeSize == 0) return;
    nodeIntervalArray = (REAL*)malloc((treeSize + 1) * 2 * sizeof(REAL));

    #pragma omp parallel for
    for (int i = 0; i < numLeaves; ++i) {
        int nodeIdx = i + numLeaves;
        nodeIntervalArray[nodeIdx * 2] = elementaryEdgeArray[i * 2];
        nodeIntervalArray[nodeIdx * 2 + 1] = elementaryEdgeArray[i * 2 + 1];
    }

    for (int i = numLeaves - 1; i > 0; --i) {
        nodeIntervalArray[i * 2] = nodeIntervalArray[(i * 2) * 2];
        nodeIntervalArray[i * 2 + 1] = nodeIntervalArray[(i * 2 + 1) * 2 + 1];
    }
}

int SegTree::find_leaf_for_point(REAL p) {
    if (numLeaves == 0) return -1;
    int l = 0, r = numLeaves;
    while (r - l > 1) {
        int mid = l + (r - l) / 2;
        if (elementaryEdgeArray[mid * 2] > p) r = mid;
        else l = mid;
    }
    return numLeaves + l;
}

void SegTree::insertAllToEdgeListMulticore(Edge *edges, int bSize) {
    #pragma omp parallel for
    for (int eid = 0; eid < intervalCount; ++eid) {
        const Edge& e = edges[eid];
        // CORRECTED: Call to helper function is now a member call
        int start_node_idx = find_leaf_for_point(e.start);
        if (start_node_idx > 0 && start_node_idx <= treeSize) {
            endList[start_node_idx].insert(e);
        }
        int end_node_idx = find_leaf_for_point(e.end);
        if (end_node_idx > 0 && end_node_idx <= treeSize) {
            endList[end_node_idx].insert(e);
        }
    }

    for (int i = numLeaves - 1; i > 0; --i) {
        endList[i] = pam_edge_set::map_union(endList[i * 2], endList[i * 2 + 1]);
    }
}

void SegTree::insertEdgeParallel(const Edge& x, int nodeIdx) {
    // CORRECTED: Create a temporary Edge object for comparison
    Edge node_interval(nodeIntervalArray[nodeIdx * 2], nodeIntervalArray[nodeIdx * 2 + 1]);
    
    // CORRECTED: Pass the Edge object to the 'contains' and 'intersects' methods
    if (x.contains(node_interval)) {
        coverList[nodeIdx].insert(x);
        return;
    }
    
    int leftIdx = 2 * nodeIdx;
    if (leftIdx <= treeSize) {
        Edge left_interval(nodeIntervalArray[leftIdx * 2], nodeIntervalArray[leftIdx * 2 + 1]);
        if (x.intersects(left_interval)) {
            insertEdgeParallel(x, leftIdx);
        }
    }
    
    int rightIdx = 2 * nodeIdx + 1;
    if (rightIdx <= treeSize) {
        Edge right_interval(nodeIntervalArray[rightIdx * 2], nodeIntervalArray[rightIdx * 2 + 1]);
        if (x.intersects(right_interval)) {
            insertEdgeParallel(x, rightIdx);
        }
    }
}

void SegTree::insertAllParallel(Edge* edges) {
    #pragma omp parallel for
    for (int eid = 0; eid < intervalCount; ++eid) {
        insertEdgeParallel(edges[eid], 1);
    }
}


// --- Core Query Function ---

void SegTree::saveIntersectionsIdsMulticore(
    REAL *bPoly, REAL *cPoly, int bSize, int cSize,
    int bType, int cType,
    std::vector<int> &bPolLineIds, std::vector<int> &cPolLineIds, std::vector<int> &intersectTypesDup) {

    omp_lock_t write_lock;
    omp_init_lock(&write_lock);

    #pragma omp parallel for schedule(dynamic, 100)
    for (int nid = 1; nid <= treeSize; nid++) {
        if (coverList[nid].size() == 0 && endList[nid].size() == 0) {
            continue;
        }

        // CORRECTED: Use the free function pam::to_sequence
        parlay::sequence<Edge> current_cover_list = pam::to_sequence(coverList[nid]);
        parlay::sequence<Edge> current_end_list = pam::to_sequence(endList[nid]);
        
        // ... [Intersection logic remains the same] ...
        for (size_t i = 0; i < current_cover_list.size(); ++i) {
            for (size_t j = i + 1; j < current_cover_list.size(); ++j) {
                // Perform geometric intersection test for cover_list[i] vs cover_list[j]
                // On intersection, lock and push to result vectors.
            }
        }
        
        for (const auto& cover_edge : current_cover_list) {
            for (const auto& end_edge : current_end_list) {
                // Perform geometric intersection test for cover_edge vs end_edge
                // On intersection, lock and push to result vectors.
            }
        }
    }

    omp_destroy_lock(&write_lock);
}

