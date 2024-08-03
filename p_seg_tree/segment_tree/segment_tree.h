#ifndef SEGTREE_H
#define SEGTREE_H

#include <iostream>
#include <algorithm> // std::min
#include <set>
#include <map>
#include <utility>
#include <math.h> /* log2 */
#include <omp.h>

#include "node.h"
#include "../util/util.h"

typedef struct dsintersection
{
    int type;
    int bPolId;
    int cPolId;
} dsintersection;

class SegTree
{
private:
    void recQueryTree(REAL qx, Node *root, int nodeIdx, vector<Edge> *edges);
    // void recQueryTreeIterative(REAL qx, Node *root, int nodeIdx, vector<Edge> *edges);

    void insertToLeafLevelEdgeList(Edge *edges);

    map<REAL, pair<Node *, Node *>> *preparePointsToNodeMap();

    void buildInternalEndLists(int node, int start, int end);
    void buildInternalEndLists();

public:
    int intervalCount;
    int treeSize;
    int numLeaves;
    int numEdges;
    int elementaryPointCount;
    // int elementaryArrayCount;
    int elementaryEdgeCount;
    int intervalCountUpperBound;
    int POLY_TYPE_COUNT;

    REAL minInterval, maxInterval;

    vector<REAL> *m_elementaryPoints;
    vector<Edge> *elementaryVec;
    // this holds the skeleton of the segment-tree
    vector<Node *> *treeNode;
    omp_lock_t *lock;

    // array version of the tree end lists
    REAL *elementaryArray;
    // Edge *elementaryEdgeArray;
    REAL *elementaryEdgeArray;
    Edge *treeNodeEndListArray;
    Edge *treeNodeCoverListArray;
    Edge *intervals2;


    // duplicate array representation of cover lists in nodes using 2-D array
    REAL *coverArray;
    // following arrays store coverArray meta data
    unsigned long sumCoverListSize; // #nodes*maxsize(coverlist)
    int *coverListSizes;            // sizes of the cover lists at each node
    // int *coverListSizes2; //sizes of the cover lists at each node
    int *coverListCounts;
    REAL *nodeIntervalArray; // intervals of each node

    // duplicate array representation of end lists in nodes using 2-D array
    REAL *endArray;
    // following arrays store coverArray meta data
    unsigned long sumEndListSize; // #nodes*maxsize(coverlist)
    int *endListSizes;            // sizes of the cover lists at each node
    int *endListCounts;
    REAL delta;

    /*void init2(vector<Edge> *aVec)
    {
        elementaryVec = aVec;
        numLeaves = aVec->size();
        treeSize = 2 * aVec->size();

        treeNode = new vector<Node *>(treeSize);

        // coverListSizes = new int[treeSize+1];  //+1 is used for prefix sum later
        // nodeIntervalArray = new REAL[treeSize*2];
        // endListSizes = new int[treeSize+1];   //+1 is used for prefix sum later
        // endListCounts = new int[treeSize];

        coverListSizes = (int *)malloc((treeSize + 1) * sizeof(int));
        // coverListSizes2=(int *)malloc((treeSize+1)*sizeof(int));
        coverListCounts = (int *)malloc((treeSize) * sizeof(int));

        endListSizes = (int *)malloc((treeSize + 1) * sizeof(int));
        endListCounts = (int *)malloc((treeSize) * sizeof(int));
        nodeIntervalArray = (REAL *)malloc((treeSize * 2) * sizeof(REAL));

        // int ts=treeSize;
        // #pragma acc enter data copyin(this)
        // #pragma acc enter data create(coverListSizes[0:ts])
        // #pragma acc enter data create(nodeIntervalArray[0:ts*2])
        // #pragma acc enter data create(endListSizes[0:ts])
        // for(Edge x: *elementaryVec){
        // cout<<"ele Edge "<<x.start<<" "<<x.end<<endl;
        // }

        for (int i = 0; i < treeSize; i++)
        {
            treeNode->at(i) = new Node();
            //    coverListSizes[i]=0;
            //    coverListSizes2[i]=0;
            //    coverListCounts[i]=0;
            endListSizes[i] = 0;
            endListCounts[i] = 0;
        }

        // if(DEBUG_MODE) cout<<"Tree size "<<treeNode->size()<<endl;
    }*/

    // void init(vector<Edge> *aVec){
    void init()
    {
        // elementaryVec = aVec;
        // numLeaves = aVec->size();
        // treeSize = 2*aVec->size();

        treeSize = 2 * numLeaves;
        treeNode = new vector<Node *>(treeSize);

        // coverListSizes = new int[treeSize+1];  //+1 is used for prefix sum later
        // nodeIntervalArray = new REAL[treeSize*2];
        // endListSizes = new int[treeSize+1];   //+1 is used for prefix sum later
        // endListCounts = new int[treeSize];

        coverListSizes = (int *)malloc((treeSize + 1) * sizeof(int));
        // coverListSizes2=(int *)malloc((treeSize+1)*sizeof(int));
        coverListCounts = (int *)malloc((treeSize) * sizeof(int));

        endListSizes = (int *)malloc((treeSize + 1) * sizeof(int));
        // endListSizes = (long *)malloc((treeSize + 1) * sizeof(long));
        endListCounts = (int *)malloc((treeSize) * sizeof(int));
        nodeIntervalArray = (REAL *)malloc((treeSize * 2) * sizeof(REAL));

        // int lts=1;
        // int lts=log2(treeSize);
        int lts = treeSize;
        lock = (omp_lock_t *)malloc((lts) * sizeof(omp_lock_t));

        // #pragma omp parallel for
        // for (int i = 0; i < lts; ++i)
        // {
        //     omp_init_lock(&lock[i]);
        // }

        // int ts=treeSize;
        // // #pragma acc enter data copyin(this)
        // #pragma acc enter data create(coverListSizes[0:ts])
        // #pragma acc enter data create(nodeIntervalArray[0:ts*2])
        // #pragma acc enter data create(endListSizes[0:ts])
        // for(Edge x: *elementaryVec){
        // cout<<"ele Edge "<<x.start<<" "<<x.end<<endl;
        // }
        #pragma omp parallel for
        for (int i = 0; i < treeSize; i++)
        {
            // omp_init_lock(&lock[i]);

            // treeNode->at(i) = new Node();
            coverListSizes[i] = 0;
            //    coverListSizes2[i]=0;
            coverListCounts[i] = 0;
            endListSizes[i] = 0;
            endListCounts[i] = 0;
            // if(i==treeSize)  endListSizes[i+1] = 0;
        }

        // if(DEBUG_MODE) cout<<"Tree size "<<treeNode->size()<<endl;
    }

    void initHybrid()
    {
        treeSize = 2 * numLeaves;
        treeNode = new vector<Node *>(treeSize);
        coverListSizes = (int *)malloc((treeSize + 1) * sizeof(int));
        endListSizes = (int *)malloc((treeSize + 1) * sizeof(int));
        endListCounts = (int *)malloc((treeSize) * sizeof(int));
        nodeIntervalArray = (REAL *)malloc((treeSize * 2) * sizeof(REAL));
        int ts=treeSize;

        #pragma acc enter data create(endListSizes[0:ts+1], endListCounts[0:ts])
        #pragma acc parallel loop
        for (int i=0; i<ts; i++){
            endListSizes[i]=0;
            endListCounts[i]=0;
        }
    }
    SegTree()
    {
        numEdges = 0;
    }

    SegTree(int aNumEdges, int ic, REAL minVal, REAL maxVal)
    {
        numEdges = aNumEdges;
        intervalCount = ic;
        minInterval = minVal;
        maxInterval = maxVal;
    }

    ~SegTree()
    {
        free(treeNode);
        free(coverListSizes);
        free(coverListCounts);
        free(endListSizes);
        free(endListCounts);
        free(nodeIntervalArray);
    }

    Edge seg_union(Edge a, Edge b)
    {
        double min = std::min(a.start, b.start);
        double max = std::max(a.end, b.end);

        return Edge(min, max);
    }

    void intervalConstructionParallel(REAL *bPoly, REAL *cPoly, int bSize, int cSize,
                                      Edge **intervals);

    vector<REAL> *getPointsFromEdges(vector<Edge> *intervals);

    vector<Edge> *getElementaryIntervals(Edge *intervals);
    vector<REAL> *getPointsFromEdges2(Edge *intervals);

    void buildSkeleton(int node, int start, int end);

    void buildSkeleton2(int node, int start, int end);

    // https://stackoverflow.com/questions/20922755/print-heap-array-in-tree-format
    void display();

    void copyInputDataToGPU(REAL *bPoly, REAL *cPoly, int bSize, int cSize, int cPolysSize, int cPolyCount, int ic, Edge *edges);

    void insertParallel(Edge e);

    void insertEdgeParallel(Edge x, int nodeIdx);

    void insertAllParallel(Edge *edges);

    void insertAllParallel2(Edge *edges);

    void constructCoverListGPU();

    void constructCoverListMulticore(Edge *edges);

    void constructCoverListSerial(Edge *edges);

    void constructCoverEndLists(REAL *bPoly, REAL *cPoly, int bSize, int cPolysSize);

    void insert(Edge e);

    void insertEdge(Edge x, Node *root, int nodeIdx);

    void insertAll(Edge *edges);

    vector<Edge> *querySegTree(REAL qx);

    void insertAllToEdgeList(Edge *edges);

    void copyToEleArray(vector<Edge> *elementaryIntervals);

    void constructElementaryIntervalsArrayParallel(Edge *intervals);

    void constructElementaryIntervalsArrayGPU();

    void constructElementaryIntervalsArrayMulticore(Edge *intervals);

    void constructElementaryIntervalsArrayGPUParallel(Edge *intervals);

    void constructElementaryIntervalsArrayCMBR(Edge *intervals, REAL *cmbr);

    void constructTreeIntervalArray();

    void constructTreeIntervalArrayMulticore();

    void preparePointsToNodeMapParallel(REAL *boaders, LONG *intBoader, int *left, int *right);

    void preparePointsToNodeMapMulticore(REAL *boaders, LONG *intBoaders, int *left, int *right);

    void copyLeavesFromVector();

    // void insertEdgeParallel(Edge x);

    void insertToLeafLevelEdgeListParallel(REAL *boaders, LONG *intBoaders, int *left, int *right, int bSize);

    void insertToLeafLevelEdgeListMulticore(Edge *edges, REAL *boaders, LONG *intBoaders, int *left, int *right, int bSize);

    void buildInternalEndListsParallel(int bSize);

    void buildInternalEndListsGPU(int bSize);


    void buildInternalEndListsMulticore(int bSize);

    void insertAllToEdgeListParallel(int bSize);

    void insertAllToEdgeListMulticore(Edge *edges, int bSize);

    void displayEndLists();

    void displayEndSets();

    void saveIntersectionsIds(REAL *bPoly, REAL *cPoly, int bSize, int cSize,
                              int bType, int cType,
                              vector<int> &bPolLineIds, vector<int> &cPolLineIds, vector<int> &intersectTypesDup);

    void saveIntersectionsIdsMulticore(REAL *bPoly, REAL *cPoly, int bSize, int cSize,
                                       int bType, int cType,
                                       vector<int> &bPolLineIds, vector<int> &cPolLineIds, vector<int> &intersectTypesDup);
    // void saveIntersectionsIdsMulticore
    //   (REAL *bPoly, REAL *cPoly, int bSize, int cSize,
    //   int bType, int cType,
    //   list<int> &bPolLineIds, list<int> &cPolLineIds, list<int> &intersectTypesDup);

    void saveIntersectionsIdsMulticore2(REAL *bPoly, REAL *cPoly, int bSize, int cSize,
                                                  int bType, int cType,
                                                  vector<int> &bPolLineIds, vector<int> &cPolLineIds, vector<int> &intersectTypesDup);

    void saveIntersections(
        REAL *bPoly, REAL *cPoly,
        int *bPolyLineIds, int *cPolyLineIds, int polySize,
        REAL *intersections, int *intersectTypes,
        int *intersectAlphaValuesB, int *intersectAlphaValuesC);

    void sortIntersectionsAtEdges(int *intersectAlphaValues, int *sortedIndicies, int cStart,
                                  int *allIntersectCounts, int *nonDegenIntersectCounts, int *polLineUniqueIds, int uniqueIntersectCount, int intersectCount,
                                  int *sortedAlphaIndicies);

    void copyIntersections(
        REAL *bPoly, REAL *cPoly, int bSize, int cSize, int cStart,
        int *bPolyLineIds, int *cPolyLineIds, int polySize,
        int *bPolLineIdsSortedIndex, int *cPolLineIdsSortedIndex,
        int *sortedIndiciesB, int *sortedIndiciesC,
        REAL *intersections, int *intersectTypes,
        int *bAllIntersectCounts, int *cAllIntersectCounts,
        int *bNonDegenIntersectCountPS, int *cNonDegenIntersectCountPS,
        int *intersectAlphaValuesB, int *intersectAlphaValuesC,
        REAL *intersectedB, REAL *intersectedC,
        int *neighborsB, int *neighborsC,
        int *alphaValuesB, int *alphaValuesC);

    void calculateInitLabel(int bSize, int *bNonDegenIntersectCounts,
                            REAL *intersectedB, REAL *intersectedC, int *alphaValuesB,
                            int *neighborsB,
                            int bIntersectedCount, int cIntersectedCount, int *bInitLabels, int *cInitLabels);

    void copyNodeIntervalsToGPU();

    void copyNodeCoverListToGPU(REAL *bPoly, REAL *cPoly, int bSize, int cSize);

    void countEndlistCoverListSizes();

    void countEndlistCoverListSizesMulticore();

    // void copyNodeEndListToGPU();

    void copyNodeEndEdgeSetToGPU();

    void copyNodeCoverListMulticore();

    // void recQueryArrayParallel(REAL qx);

    void countIntersectionsParallel(REAL *bPoly, REAL *cPoly, int bSize, int cSize, int cPolyCount,
                                    int bType, int cType,
                                    int *iCount);

    void countIntersectionsMulticore(REAL *bPoly, REAL *cPoly, int bSize, int cSize,
                                     int bType, int cType,
                                     int *iCount);

    void polyClipArrayParallel(REAL *bPoly, REAL *sPoly, int bSize, int cSize,
                               int bType, int cType,
                               int *iCount, int *bPolLineIds, int *cPolLineIds, int *intersectTypes);

    void polyClipArrayGPUParallel(REAL *bPoly, REAL *cPoly, int bSize, int cSize,
                                  int bType, int cType,
                                  int *iCount, int *bPolLineIds, int *cPolLineIds, int *intersectTypes);

    void countLineSegIntersectsAtNodes(int *iCount);

    void saveIntersectionsParallel(
        REAL *bPol, REAL *cPol,
        int *bPolyLineIds, int *cPolyLineIds, int polySize,
        REAL *intersections, int *intersectTypes,
        int *intersectAlphaValuesB, int *intersectAlphaValuesC);

    void sortIntersectionsAtEdgesParallel(int *intersectAlphaValues, int *sortedIndicies, int cStart,
                                          int *allIntersectCounts, int *nonDegenIntersectCounts, int *bPolLineUniqueIds, int uniqueIntersectCount, int intersectCount,
                                          int *sortedAlphaIndicies);

    void copyIntersectionsParallel(REAL *bPoly, REAL *cPoly, int bSize, int cSize, int cStart,
                                   int *bPolyLineIds, int *cPolyLineIds, int polySize,
                                   int *bPolLineIdsSortedIndex, int *cPolLineIdsSortedIndex,
                                   int *sortedBPolyIdIndicies, int *sortedCPolyIdIndicies,
                                   REAL *intersections, int *intersectTypes,
                                   int *bAllIntersectCounts, int *cAllIntersectCounts,
                                   int *bNonDegenIntersectCountPS, int *cNonDegenIntersectCountPS,
                                   int *intersectAlphaValuesB, int *intersectAlphaValuesC,
                                   REAL *intersectedB, REAL *intersectedC,
                                   int *neighborsB, int *neighborsC,
                                   int *alphaValuesB, int *alphaValuesC);

    void calculateInitLabelParallel(int bSize, int *bNonDegenIntersectCounts,
                                    REAL *intersectedB, REAL *intersectedC, int *alphaValuesB,
                                    int *neighborsB,
                                    int bIntersectedCount, int cIntersectedCount, int *bInitLabels, int *cInitLabels);

    void neighborMapParallel(
        REAL **bPol, REAL **cPol, int bSize, int cSize,
        int *bPolyLineIds, int *cPolyLineIds, int polySize,
        int *bPS1, int *cPS1, int *cPS2,
        int *neighborMapQ);

    void saveIntersectionsParallelGH(
        REAL **bPol, REAL **cPol, int bSize, int cSize,
        int *bPolyLineIds, int *cPolyLineIds, int polySize,
        int *bPS1, int *cPS1, int *bPS2, int *cPS2,
        int *neighborMapQ,
        REAL *intersectBPol, REAL *intersectCPol,
        int *alphaValuesB, int *alphaValuesC,
        int *neighborB, int *neighborC);

    // void recQueryTree(REAL qx, Node root, int nodeIdx, vector<Edge> *edges);

    // void insertAllToEdgeList2(vector<Edge> *edges);
    // void insertEdgeToEdegSet(Edge x, Node *root, int nodeIdx);
    // void insertToEdgeSet(Edge x);
};
#endif