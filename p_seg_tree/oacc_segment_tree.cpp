#include <omp.h>
#include "segment_tree/segment_tree.h"
#include <thrust/scan.h>

#include <iomanip>
// #include <iostream>
#include <chrono>

#include "util/thrust_func.h"
#include <thrust/sort.h>

using namespace std::chrono;

// compile: g++ -std=c++0x -o seg segment_tree.cpp

// original source link:
// https://www.hackerearth.com/practice/data-structures/advanced-data-structures/segment-trees/tutorial/
void SegTree ::displayEndLists()
{
   // int firstLeafIndex = treeSize/2;
   int lastLeafIndex = treeSize - 1;
   // for(int nodeIdx = firstLeafIndex; nodeIdx <= lastLeafIndex; nodeIdx++)
   for (int nodeIdx = 1; nodeIdx <= lastLeafIndex; nodeIdx++)
   {
      Node *leaf = treeNode->at(nodeIdx);
      cout << " [" << nodeIdx << "] Interval: " << leaf->interval.toString() << " size: " << leaf->endList->size() << ", ";
      // for(Edge e: *leaf->endList)
      for (Edge e : *leaf->endList)
      {
         cout << e.toString() << ", ";
      }
      cout << endl;
   }
}

void SegTree ::displayEndSets()
{
   // int firstLeafIndex = treeSize/2;
   int lastLeafIndex = treeSize - 1;
   // for(int nodeIdx = firstLeafIndex; nodeIdx <= lastLeafIndex; nodeIdx++)
   for (int nodeIdx = 1; nodeIdx <= lastLeafIndex; nodeIdx++)
   {
      Node *leaf = treeNode->at(nodeIdx);
      cout << "[" << nodeIdx << "] Interval: " << leaf->interval.toString() << " size=" << leaf->endEdgeSet->size() << endl;
      // for(Edge e: *leaf->endList)
      // for(int j=0; j<leaf->endEdgeSet->size(); ++j)
      for (Edge e : *leaf->endEdgeSet)
      {
         cout << e.toString() << ", ";
         // cout<<leaf->endEdgeSet.si<<", ";
      }
      cout << endl;
      cout << endl;
   }
}

// displays interval and cover list for each node
void SegTree ::display()
{
   float numLevels = log2(treeSize);
   int nLevels = static_cast<int>(numLevels);
   cout << "float_#Levels " << numLevels << " int_#Levels = " << nLevels << endl;

   // for (int l = 1; l < treeSize; l++)
   // {
   //    // cout<<l<<" "<<treeNode[l].interval.toString()<<endl;
   //    cout<<"["<<l<<"] Interval: "<<treeNode->at(l)->interval.toString()<< " size="<<treeNode->at(l)->coverList->size()<<endl;
   //    // vector<Edge> *coverList = treeNode->at(l)->coverList;
   //    vector<Edge> *coverList = treeNode->at(l)->coverList;

   //    for (Edge e : *coverList)
   //       cout<<e.toString()<<", ";
   //    cout << endl;
   //    cout << endl;
   // }
   for (int l = 1; l < treeSize; l++)
   {
      // cout<<l<<" "<<treeNode[l].interval.toString()<<endl;
      cout << "[" << l << "] Interval: " << treeNode->at(l)->interval.toString() << " size=" << treeNode->at(l)->coverList->size() << endl;
      // vector<Edge> *coverList = treeNode->at(l)->coverList;
      for (auto &e : *treeNode->at(l)->coverList)
         cout << e.toString() << ", ";
      cout << endl;
      cout << endl;
   }
}

// OpenMP Quick sort
// https://github.com/koszio/QuickSort-with-OpenMP/blob/master/quickSortOpenMP.c
// ----------------------------------------------------
void swap(int *a, int *b)
{
   int t = *a;
   *a = *b;
   *b = t;
}

/* This function takes last element as pivot, places
the pivot element at its correct position in sorted
   array, and places all smaller (smaller than pivot)
to left of pivot and all greater elements to right
of pivot */
int partition(int arr[], int low, int high)
{
   int pivot = arr[high]; // pivot
   int i = (low - 1);     // Index of smaller element
   for (int j = low; j <= high - 1; j++)
   {
      // If current element is smaller than or
      // equal to pivot
      if (arr[j] <= pivot)
      {
         i++; // increment index of smaller element
         swap(&arr[i], &arr[j]);
      }
   }
   swap(&arr[i + 1], &arr[high]);
   return (i + 1);
}

/* The main function that implements QuickSort
arr[] --> Array to be sorted,
low --> Starting index,
high --> Ending index */
void quickSort(int arr[], int low, int high)
{
   if (low < high)
   {
      /* pi is partitioning index, arr[p] is now
      at right place */
      int pi = partition(arr, low, high);

      // As I capture each task,  its going to look at the value of  arr , low, pi at that point and
      //  it will create a private values for that variables and its going to initialize it the current values.
      // So each task will be execute its own copy.
#pragma omp task firstprivate(arr, low, pi)
      {
         quickSort(arr, low, pi - 1);
      }

      // #pragma omp task firstprivate(arr, high,pi)
      {
         quickSort(arr, pi + 1, high);
      }
   }
}

void swapGetIndexArray(REAL *a, REAL *b, int *ia, int *ib)
{
   REAL t = *a;
   int it = *ia;
   *a = *b;
   *ia = *ib;
   *b = t;
   *ib = it;
}
int partitionGetIndexArray(REAL arr[], int *sortedIndex, int low, int high)
{
   int pivot = arr[high]; // pivot
   int i = (low - 1);     // Index of smaller element
   for (int j = low; j <= high - 1; j++)
   {
      // If current element is smaller than or
      // equal to pivot
      if (arr[j] <= pivot)
      {
         i++; // increment index of smaller element
         swapGetIndexArray(&arr[i], &arr[j], &sortedIndex[i], &sortedIndex[j]);
      }
   }
   swapGetIndexArray(&arr[i + 1], &arr[high], &sortedIndex[i + 1], &sortedIndex[high]);
   return (i + 1);
}

/* The main function that implements QuickSort
arr[] --> Array to be sorted,
low --> Starting index,
high --> Ending index */
void quickSortGetIndexArray(REAL arr[], int *sortedIndex, int low, int high)
{
   if (low < high)
   {
      /* pi is partitioning index, arr[p] is now
      at right place */
      int pi = partitionGetIndexArray(arr, sortedIndex, low, high);
      // As I capture each task,  its going to look at the value of  arr , low, pi at that point and
      //  it will create a private values for that variables and its going to initialize it the current values.
      // So each task will be execute its own copy.
#pragma omp task firstprivate(arr, low, pi)
      {
         quickSortGetIndexArray(arr, sortedIndex, low, pi - 1);
      }

      // #pragma omp task firstprivate(arr, high,pi)
      {
         quickSortGetIndexArray(arr, sortedIndex, pi + 1, high);
      }
   }
}
// -----------------------------------------------------

// #pragma acc routine seq
LONG realToInt(REAL val, int positions)
{
   return (int)pow(10, positions) * val;
}

void SegTree ::buildSkeleton(int node, int start, int end)
{
   int nl = numLeaves;
   // cout<<"hereeee"<<endl;
   // #pragma acc update self(elementaryEdgeArray[0:nl])
   if (start == end)
   {
      // Leaf node will have a single element
      // treeNode->at(node)->interval=elementaryVec->at(start);
      treeNode->at(node)->interval = elementaryEdgeArray[start];
   }
   else
   {
      int mid = (start + end) / 2;
      // Recurse on the left child
      buildSkeleton(2 * node, start, mid);
      // Recurse on the right child
      buildSkeleton(2 * node + 1, mid + 1, end);
      // Internal node will have the sum of both of its children
      treeNode->at(node)->interval = seg_union(treeNode->at(2 * node)->interval, treeNode->at(2 * node + 1)->interval);
   }
}

void SegTree ::buildSkeleton2(int node, int start, int end)
{
   if (start == end)
   {
      // Leaf node will have a single element
      treeNode->at(node)->interval = elementaryVec->at(start);
   }
   else
   {
      int mid = (start + end) / 2;
      // Recurse on the left child
      buildSkeleton2(2 * node, start, mid);
      // Recurse on the right child
      buildSkeleton2(2 * node + 1, mid + 1, end);
      // Internal node will have the sum of both of its children
      treeNode->at(node)->interval = seg_union(treeNode->at(2 * node)->interval, treeNode->at(2 * node + 1)->interval);
   }
}

int getDigitCount(REAL val)
{
   int digitCount = 1;
   // int intVal=realToInt(val, MAX_DECIMAL_COUNT);
   int intVal = trunc(val);
   // cout<<"Max "<<val<<" intVal="<<intVal<<endl;
   while (intVal / 10 > 0)
   {
      intVal /= 10;
      digitCount++;
   }
   return digitCount;
}

// Buddhi
map<REAL, pair<Node *, Node *>> *SegTree ::preparePointsToNodeMap()
{
   if (m_elementaryPoints == NULL && false == m_elementaryPoints->empty())
   {
      cout << "Error: m_elementaryPoints in insertToLeafLevelEndList" << endl;
      return NULL;
   }
   map<REAL, pair<Node *, Node *>> *ptToAdjNodeMap = new map<REAL, pair<Node *, Node *>>();
   int firstLeafIndex = treeSize / 2;

   int leafIndex = firstLeafIndex;
   int i;
   for (i = 1; i < m_elementaryPoints->size() - 2; i++, leafIndex++)
   {
      Node *leftLeaf = treeNode->at(leafIndex);
      Node *rightLeaf = treeNode->at(leafIndex + 1);

      // cout<<i<<" m_elemetry points "<<m_elementaryPoints->at(i)<<" left: "<<leafIndex<<" right:"<<leafIndex+1<<" "<<m_elementaryPoints->at(i)<<endl;

      if (i != 1 && m_elementaryPoints->at(i) == m_elementaryPoints->at(i + 1))
      {
         rightLeaf = treeNode->at(leafIndex + 2);
         // cout<<i<<" **m_elemetry points "<<m_elementaryPoints->at(i)<<" left: "<<leafIndex<<" right:"<<leafIndex+2<<" "<<m_elementaryPoints->at(i)<<endl;
      }

      pair<Node *, Node *> ithNodePair = make_pair(leftLeaf, rightLeaf);

      REAL ithPt = m_elementaryPoints->at(i);
      pair<REAL, pair<Node *, Node *>> kv = make_pair(ithPt, ithNodePair);
      ptToAdjNodeMap->insert(kv);

      // cout<<"*** "<<i<<" "<<ithPt<<endl;
   }
   return ptToAdjNodeMap;
}

// Buddhi
void SegTree ::insertToLeafLevelEdgeList(Edge *edges)
{
   map<REAL, pair<Node *, Node *>> *ptToNodeMap = preparePointsToNodeMap();

   if (ptToNodeMap == NULL)
   {
      cout << "ERROR: map empty" << endl;
      exit(1);
   }
   else
   {         
      #ifdef DG
      cout << "Map Size " << ptToNodeMap->size() << endl;
      #endif
   }

   std::map<REAL, pair<Node *, Node *>>::iterator it;
   // cout<<"before for loop all edges "<<endl;
   Edge e;
   for (int eid = 0; eid < intervalCount; ++eid)
   {
      e = edges[eid];
      // end point
      it = ptToNodeMap->find(e.end);
      if (it != ptToNodeMap->end())
      {
         pair<Node *, Node *> nodePair = ptToNodeMap->at(e.end);
         // pair<Node *, Node *> nodePair = it->second;
         Node *after = nodePair.second;
         if (after->interval.start != after->interval.end)
         {
            // cout<<"** "<<e.toString()<<" end: "<<it->first<<" after "<<after->interval.toString()<<endl;
            after->addToEndList(e);
            after->addToEndEdgeSet(e);
         }
      }

      // start point
      it = ptToNodeMap->find(e.start);
      if (it != ptToNodeMap->end())
      {
         pair<Node *, Node *> nodePair = ptToNodeMap->at(e.start);
         // pair<Node *, Node *> nodePair = it->second;
         Node *before = nodePair.first;
         if (before->interval.start != before->interval.end)
         {
            // cout<<" start: "<<it->first<<" before "<<before->interval.toString()<<endl;
            before->addToEndList(e);
            before->addToEndEdgeSet(e);
         }
      }
   }
}

void SegTree ::buildInternalEndLists()
{
   int firstLeafIndex = treeSize / 2;
   int lastLeafIndex = treeSize - 1;

   buildInternalEndLists(1, firstLeafIndex, lastLeafIndex);
}

void SegTree ::buildInternalEndLists(int node, int start, int end)
{
   // cout<<"node "<<node<<" start: "<<start<<" end: "<<end<<endl;
   // cout<<"Ele size "<<elementaryVec->size()<<endl;
   if (start == end)
   {
      // Leaf node will have a single element
      // treeNode->at(node)->interval = elementaryVec->at(start);
      // set<Edge,edge_compare> *endSet = treeNode->at(node)->endEdgeSet;
      // vector<Edge> *lEndList = new vector<Edge>( endSet->begin(), endSet->end());
      // cout<<"test node "<<start<<endl<<endl;
      // cout<<elementaryVec->at(start).start<<endl;
   }
   else
   {
      int mid = (start + end) / 2;
      // Recurse on the left child
      buildInternalEndLists(2 * node, start, mid);
      // Recurse on the right child
      buildInternalEndLists(2 * node + 1, mid + 1, end);

      // cout<<"start="<<start<<" end="<<end<<" node="<<node<<endl;
      // Internal node will have the concatenation of both of its children
      set<Edge, edge_compare> *left = treeNode->at(2 * node)->endEdgeSet;
      set<Edge, edge_compare> *right = treeNode->at(2 * node + 1)->endEdgeSet;

      set<Edge, edge_compare> *root = treeNode->at(node)->endEdgeSet;
      // cout<<"root size before: "<<root->size()<<endl;
      root->insert(left->begin(), left->end());
      root->insert(right->begin(), right->end());
   }
}

// original
void SegTree ::insertAllToEdgeList(Edge *edges)
{
   #ifdef DG
   high_resolution_clock::time_point cp1, cp2;
   cp1 = high_resolution_clock::now();
   #endif
   insertToLeafLevelEdgeList(edges);
   #ifdef DG
   cp2 = high_resolution_clock::now();
   #endif

   #ifdef DG
   auto d1 = duration_cast<microseconds>(cp2 - cp1);
   cout<<"leaf build time: "<<d1.count()<<endl;
   #endif
   buildInternalEndLists();
}

// -----------------------------------
// https://stackoverflow.com/questions/6443569/implementation-of-c-lower-bound#:~:text=lower_bound%20is%20almost%20like%20doing%20a%20usual%20binary,return%20a%20pointer%2Fiterator%20to%20the%20first%20matching%20element.
// binary search that returns the first occurrance
// #pragma acc routine seq
extern "C" __device__ int findBoarder(REAL *boaders, int *sortedIndex, int epc, REAL query)
{
   int l = 0, h = epc;  // is epc-1??
   while (l < h)
   {
      int m = l + (h - l) / 2;
      // Check if x is present at mid
      if (boaders[sortedIndex[m]] == query)
      {
         // if(boaders[m]==query){
         // return m;
         // when there is a match, we should keep on searching for the next same element. If the same element is not  found, m is considered as the answer
         if (boaders[sortedIndex[m - 1]] != query)
         {
            // if(boaders[m-1]!=query){
            h = m;
            l = h;
         }
         else
         {
            h = m - 1;
         }
      }
      else if (boaders[sortedIndex[m]] < query)
      { // If x greater, ignore left half
         // }else if(boaders[m]<query){ // If x greater, ignore left half
         l = m + 1;
      }
      else
      { // If x is smaller, ignore right half
         h = m - 1;
      }
   }
   // If we reach here, then element was not present
   if (boaders[sortedIndex[h]] != query)
   {
      // if(boaders[h]!=query){
      return -1;
   }
   return h;
}
#pragma acc routine seq
extern "C" int findBoarder(REAL *boaders, int *sortedIndex, int epc, REAL query);

extern "C" __device__ int findBoarderSorted(REAL *arr, int size, REAL query){
   int l=0, h=size;
   while (l < h){
      int m=l+(h-l)/2;
      // Check if x is present at mid
      if (arr[m]==query){
         // when there is a match, we should keep on searching for the next same element. If the same element is not  found, m is considered as the answer
         if (arr[m-1]!=query){
            h=m;
            l=h;
         }else{
            h=m-1;
         }
      }
      else if (arr[m]<query){ // If x greater, ignore left half
         // }else if(boaders[m]<query){ // If x greater, ignore left half
         l=m+1;
      }else{ // If x is smaller, ignore right half
         h=m-1;
      }
   }
   // If we reach here, then element was not present
   if (arr[h]!=query){
      return -1;
   }
   return h;
}
#pragma acc routine seq
extern "C" int findBoarderSorted(REAL *arr, int size, REAL query);

/*
finds the last occurrance element of the target. 
If not, returns the next largest
*/
extern "C" __device__ int lastOccurrenceBinarySearch(REAL *arr, int *sortedIndex, int size, REAL query) {
   int low=0;
   int high=size-1;
   int result=-1;

   while (low<=high) {
      int mid=low+(high-low)/2;
      if (arr[sortedIndex[mid]]==query) {
      // if (abs(arr[sortedIndex[mid]]-query)<=EPSILON) {
         result=mid;
         // result=sortedIndex[mid];
         low=mid+1; // Look in the right half
      }else if (arr[sortedIndex[mid]]<query) {
      // }else if ((query-arr[sortedIndex[mid]])>EPSILON) {
         low=mid+1;
      } else{
         result=mid;
         high=mid-1;
      }
   }

    return result;
}
#pragma acc routine seq
extern "C" int lastOccurrenceBinarySearch(REAL *arr, int *sortedIndex, int size, REAL query);

/*
finds the first occurrance element of the target. 
If not, returns the next largest
*/
extern "C" __device__ int firstOccurrenceBinarySearch(REAL *arr, int *sortedIndex, int size, REAL query) {
    int low = 0, high = size - 1;
    int result = size; // Initialize result to size of array

    while (low <= high) {
        int mid = low + (high - low) / 2;

        if (arr[sortedIndex[mid]] == query) {
            result = mid; // Update result, but continue searching on the left side
            high = mid - 1;
        } else if (arr[sortedIndex[mid]] > query) {
            result = mid; // Update result to the current element
            high = mid - 1;
        } else {
            low = mid + 1;
        }
    }

    return result;
}
#pragma acc routine seq
extern "C" int firstOccurrenceBinarySearch(REAL *arr, int *sortedIndex, int size, REAL query);

int serialFindBoader(REAL *boaders, int epc, int query)
{
   for (int i = 0; i < epc; ++i)
   {
      if (boaders[i] == query)
      {
         return i;
      }
   }
   return -1;
}
// ------------------------------------------------------------------
/*
parallel methods to prepare end lists at each edge
*/
void SegTree::preparePointsToNodeMapParallel(REAL *boaders, LONG *intBoaders, int *left, int *right)
{

   int firstLeafIndex = treeSize / 2;
   int leafIndex = firstLeafIndex;
   int i;
   // int ic=elementaryArrayCount;
   int osize = intervalCount;
   int icupper = intervalCountUpperBound;

   REAL absMinEndPoint = abs(minInterval);
   int maxDigitCount = getDigitCount(maxInterval + absMinEndPoint);
   // cout<<"preparePointsToNodeMapParallel epc="<<epc<<endl;
   // -----------------
   // #pragma acc enter data copyin(elementaryArray[0 : epc]) ///*********
   int epc = elementaryPointCount;
   // #pragma acc data copyin(this)
   #pragma acc data present(elementaryArray[0:epc])
   #pragma acc parallel loop
   for (i = 1; i < epc; i++)
   {
      int leftLeaf = leafIndex + i - 1;
      int rightLeaf = leafIndex + i;
      // handling point intervals
      // ---------------------------------------------
      if (i != 1 && elementaryArray[i] == elementaryArray[i + 1])
      {
         rightLeaf = leafIndex + i + 1;
      }
      // ---------------------------------------------
      left[i - 1] = leftLeaf;
      right[i - 1] = rightLeaf;
      boaders[i - 1] = elementaryArray[i];
      if (elementaryArray[0] < 0)
      {
         intBoaders[i - 1] = realToInt(elementaryArray[i] + absMinEndPoint, (MAX_INTERVAL_DIGIT_COUNT - maxDigitCount));
      }
      else
      {
         intBoaders[i - 1] = realToInt(elementaryArray[i], (MAX_INTERVAL_DIGIT_COUNT - maxDigitCount));
      }
   }
}

/*
Radix Sort - V2
Sorts one array and output sorted index array
indexArr: first round has the sorted order of the polygon ids
*/
// ------------------------------------------
void countSort(LONG *input, int *indexArr, int start, int end)
{
   int n = end - start;
   LONG *output;     // output array
   int *lsVals;      // save least significant digits of the input
   int *outputIndex; // current sorted index order
   output = (LONG *)malloc((n) * sizeof(LONG));
   lsVals = (int *)malloc((n) * sizeof(int));
   outputIndex = (int *)malloc((n) * sizeof(int));

   int i, j, count[10] = {0};

   // Store count of occurrences in count[]
   for (i = 0, j = start; i < n; i++, j++)
   {
      lsVals[i] = input[j] % 10;
      count[lsVals[i]]++;
      // currOutputIndex[i]=outputIndex[i];
   }
   // Change count[i] so that count[i] now contains actual
   //  position of this digit in output[]
   for (i = 1; i < 10; i++)
      count[i] += count[i - 1];

   // Build the output array
   for (i = n - 1, j = end - 1; i >= 0; i--, j--)
   {
      output[count[lsVals[i]] - 1] = input[j];
      outputIndex[count[lsVals[i]] - 1] = indexArr[j]; ///////////////Error********
      count[lsVals[i]]--;

      // cout<<"== "<<lsVals[i]<<endl;
   }

   // Copy the output array to arr[], so that arr[] now
   // contains sorted numbers according to current digit
   for (i = start, j = 0; i < end; i++, j++)
   {
      input[i] = output[j] / 10;
      indexArr[i] = outputIndex[j];
      // cout<<j<<" i*** "<<i<<" "<<indexArr[i]<<" "<<input[i]<<" out:"<<output[j]<<endl;
   }
   free(output);
   free(lsVals);
   free(outputIndex);
}

// The main function to that sorts arr[] of size n using
// Radix Sort
void radixsort(LONG *data, int *outputIndex, int start, int end, int size, int m)
{
   // copy input data to use as the output array
   // cout<<"== "<<start<<" "<<end<<" "<<end-start<<endl;
   LONG *input;
   input = (LONG *)malloc((size) * sizeof(LONG));
   for (int i = start; i < end; ++i)
   {
      input[i] = data[i];
      //  cout<<i<<" sort "<<data[i]<<", ";
   }
   // cout<<"Data copied for sorting"<<endl;
   // Do counting sort for every digit. Note that instead
   // of passing digit number, exp is passed. exp is 10^i
   // where i is current digit number
   // maxType=1 to use max # of digits
   // maxType=2 to use max value
   for (int exp = 1, i = 0; m > i; exp *= 10, ++i)
   {
      // cout<<"round "<<exp<<"\n\n"<<endl;
      // countSort(data, outputIndex, start, end);
      countSort(input, outputIndex, start, end);
   }
   free(input);
}
// ------------------------------------------

/*
copy tree nodes to GPU
*/
/*void SegTree::copyNodeIntervalsToGPU(){
   int nis=treeSize*2;

   for(int nid=0, nid2=0; nid<treeSize; ++nid, nid2+=2){
      nodeIntervalArray[nid2]=treeNode->at(nid)->interval.start;
      nodeIntervalArray[nid2+1]=treeNode->at(nid)->interval.end;
      // cout<<nid<<" node intervals "<<nodeIntervalArray[nid2]<<", "<<nodeIntervalArray[nid2+1]<<endl;
   }
   #pragma acc enter data copyin(this, nodeIntervalArray[0:nis])
}*/

#if GPU!=0
void SegTree::insertToLeafLevelEdgeListParallel(REAL *boaders, LONG *intBoaders, int *left, int *right, int bSize)
{
   int epc = elementaryPointCount;
   int *boadersSortedIndex;
   boadersSortedIndex = (int *)malloc((epc) * sizeof(int));

   // for (int i = 0; i < epc; ++i)
   //    boadersSortedIndex[i] = i;

   preparePointsToNodeMapParallel(boaders, intBoaders, left, right);
   // REad MAX_INTERVAL_DIGIT_COUNT from bSize digits

   // #pragma acc data present(boaders[0 : epc], intBoaders[0 : epc], left[0 : epc], right[0 : epc])

   // #pragma acc update self(intBoaders[0 : epc])
   // radixsort(intBoaders, boadersSortedIndex, 0, epc, epc, MAX_INTERVAL_DIGIT_COUNT);

   #pragma acc enter data create(boadersSortedIndex[0:epc])
   #pragma acc parallel loop
   for (int i = 0; i < epc; ++i)
      boadersSortedIndex[i] = i;
   
   #pragma acc host_data use_device(intBoaders, boadersSortedIndex)
   thrustDeviceStableSortByKey(intBoaders, boadersSortedIndex, epc);

   // if(ptToNodeMap == NULL){
   //    cout << "ERROR: map empty" << endl;
   //    exit(1);
   // }else{
   //    if(DEBUG) cout << "Map Size " << ptToNodeMap->size() << endl;
   // }
   int ic = intervalCount;
   // #pragma acc update self(boaders[0:epc], left[0:epc], right[0:epc])
   int *enids, *stids;
   enids = (int *)malloc((ic) * sizeof(int));
   stids = (int *)malloc((ic) * sizeof(int));

   // cout<<"done part 3321"<<endl;
   // #pragma acc enter data copyin(this)
   #pragma acc enter data create(enids[0 : ic], stids[0 : ic]) /*copyin(edges[0 : ic])*/
   // #pragma acc enter data create(enids[0:ic], stids[0:ic]) copyin(this, edges[0:ic], boadersSortedIndex[0:epc])
   #pragma acc data present(boaders[0 : epc], boadersSortedIndex[0:epc], intervals2[0:ic])

   // use seperate loop to save foud index. Then work those values later
   #pragma acc parallel loop
   for (int eid = 0; eid < ic; ++eid)
   {
      enids[eid] = findBoarder(boaders, boadersSortedIndex, epc, intervals2[eid].end);
      stids[eid] = findBoarder(boaders, boadersSortedIndex, epc, intervals2[eid].start);
   }

   int ts = treeSize;
   // #pragma acc enter data copyin(this)
   // #pragma acc enter data create(endListSizes[0 : ts + 1])
   // #pragma acc enter data copyin(endListSizes[0 : ts + 1])

   // #pragma acc enter data copyin(this, nodeIntervalArray[0:ts2])
   // #pragma acc data present(enids[0:ic], stids[0:ic], left[0:epc], right[0:epc], boadersSortedIndex[0:epc])

   #pragma acc data present(nodeIntervalArray[0:ts*2], enids[0:ic], stids[0:ic], left[0:epc], right[0:epc], boadersSortedIndex[0:epc], endListSizes[0:ts+1])
   #pragma acc parallel loop
   for (int eid = 0; eid < ic; ++eid)
   {
      int after, before;

      // end point
      if (enids[eid] >= 0)
      {
         after = right[boadersSortedIndex[enids[eid]]];
         // after=right[enids[eid]];
         if (nodeIntervalArray[after * 2] != nodeIntervalArray[after * 2 + 1])
         {
            #pragma acc atomic update
            endListSizes[after] = endListSizes[after] + 1;
         }
      }
      // start point
      if (stids[eid] >= 0)
      {
         before = left[boadersSortedIndex[stids[eid]]];
         // before=left[stids[eid]];
         if (nodeIntervalArray[before * 2] != nodeIntervalArray[before * 2 + 1])
         {
            #pragma acc atomic update
            endListSizes[before] = endListSizes[before] + 1;
         }
      }
   }
   // #pragma acc update self(endListSizes[0 : ts + 1])
   // thrust::exclusive_scan(endListSizes, endListSizes + treeSize + 1, endListSizes);

   #pragma acc host_data use_device(endListSizes)
   thrustDeviceExclusiveScan(endListSizes, ts+1);
   // #pragma acc update self(endListSizes[ts])   // only copies the last element to CPU
   #pragma acc update self(endListSizes[ts])   //  copies all elements to CPU

   // for(int nid=treeSize/2; nid<treeSize; ++nid){
   //    cout<<"nid="<<nid<<" "<<endListSizes[nid]<<" "<<treeNode->at(nid)->endEdgeSet->size()<<endl;
   // }

   sumEndListSize = endListSizes[ts];
   int sels = sumEndListSize * log2(ts);
   treeNodeEndListArray = (Edge *)malloc((sels) * sizeof(Edge));
   int leafstart = sels - sumEndListSize;

   // cout<<"done part 3333331"<<endl;

   // serial test --------------
   // #pragma acc update self(enids[0:ic], stids[0:ic], left[0:epc], right[0:epc], nodeIntervalArray[0:ts*2])
   // cout<<treeSize<<" "<<endListSizes[treeSize]<<" ++++++ "<<sels<<endl;
   // --------------

   // #pragma acc update device(endListSizes[0 : ts + 1]) *****
   // #pragma acc enter data copyin(this)
   #pragma acc enter data create(treeNodeEndListArray[0:sels]) /*copyin(endListCounts[0 : ts])*/
   #pragma acc data present(enids[0:ic], stids[0:ic], left[0:epc], right[0:epc], boadersSortedIndex[0:epc], endListSizes[0:ts+1], intervals2[0:ic], endListCounts[0:ts])
   #pragma acc parallel loop
   for (int eid = 0; eid < ic; ++eid)
   {
      int after, before, localid;

      // end point
      if (enids[eid] >= 0)
      {
         after = right[boadersSortedIndex[enids[eid]]];
         // after=right[enids[eid]];
         if (nodeIntervalArray[after * 2] != nodeIntervalArray[after * 2 + 1])
         {
            #pragma acc atomic
            {
               treeNodeEndListArray[leafstart + endListSizes[after] + (endListCounts[after] += 1)] = intervals2[eid];
               // treeNodeEndListArray[leafstart+endListSizes[after]+(endListCounts[after])]=intervals2[eid];
               // endListCounts[after]=endListCounts[after]+1;
            }
         }
      }
      // start point
      if (stids[eid] >= 0)
      {
         before = left[boadersSortedIndex[stids[eid]]];
         // before=left[stids[eid]];
         if (nodeIntervalArray[before * 2] != nodeIntervalArray[before * 2 + 1])
         {
            #pragma acc atomic update
            {
               treeNodeEndListArray[leafstart + endListSizes[before] + (endListCounts[before] += 1)] = intervals2[eid];
               // treeNodeEndListArray[leafstart+endListSizes[before]+(endListCounts[before])]=intervals2[eid];
               // endListCounts[before]=endListCounts[before]+1;
            }
         }
      }
   }
}
#endif
// ----------------------------------------
// OpenAcc radix sort
void countSort(Edge *input, Edge *output, int *lsVals, int start, int n, int exp, int byType)
{
   // Edge *output; // output array
   // int *lsVals;  //save least significant digits of the input

   // output=(Edge *)malloc((n)*sizeof(Edge));
   // lsVals=(int *)malloc((n)*sizeof(int));

   int j, i, count[10] = {0};

   // Store count of occurrences in count[]
   for (i = 0, j = start; i < n; i++, j++)
   {
      if (byType == 1)
         lsVals[j] = (input[j].type / exp) % 10;
      else
      {
         lsVals[j] = (input[j].id / exp) % 10;
      }
      count[lsVals[j]]++;
      // currOutputIndex[i]=outputIndex[i];
   }

   // Change count[i] so that count[i] now contains actual
   //  position of this digit in output[]
   for (i = 1; i < 10; i++)
      count[i] += count[i - 1];

   // Build the output array
   for (i = n - 1, j = start + i; i >= 0; i--, j--)
   {
      output[start + count[lsVals[j]] - 1] = input[j];
      count[lsVals[j]]--;
   }

   // Copy the output array to arr[], so that arr[] now
   // contains sorted numbers according to current digit
   for (i = 0, j = start; i < n; i++, j++)
   {
      input[j] = output[j];
      // cout<<"i "<<i<<" "<<output[i]<<" "<<input[i]<<endl;
   }
}

// Radix Sort
// #pragma acc routine worker
void edgeRadixsort(Edge *data, Edge *output, int *lsVals, int start, int n, int max)
{
   // sort by type
   for (int exp = 1; MAX_POLY_TYPE / exp > 0; exp *= 10)
   {
      // cout<<"start="<<start<<" n="<<n<<endl;
      countSort(data, output, lsVals, start, n, exp, 1);
   }
   // sort by id
   for (int exp = 1; max / exp > 0; exp *= 10)
   {
      countSort(data, output, lsVals, start, n, exp, 2);
   }
}

// #pragma acc routine seq
// extern "C" void countSortDevice(Edge *input, int start, int n, int exp, int byType);

// Radix Sort
extern "C" __device__ void edgeRadixsortDevice(Edge *input, Edge *output, int *lsVals, int start, int n, int max)
{
   // sort by type
   for (int exp = 1; MAX_POLY_TYPE / exp > 0; exp *= 10)
   {
      // countSortDevice(input, start, n, exp, 1);

      int j, i, count[10];
      for (int i = 0; i < 10; ++i)
      {
         count[i] = 0;
      }

      // Store count of occurrences in count[]
      for (i = 0, j = start; i < n; i++, j++)
      {
         lsVals[j] = (input[j].type / exp) % 10;
         count[lsVals[j]]++;
         // currOutputIndex[i]=outputIndex[i];
      }

      // Change count[i] so that count[i] now contains actual
      //  position of this digit in output[]
      for (i = 1; i < 10; i++)
         count[i] += count[i - 1];

      // Build the output array
      for (i = n - 1, j = start + i; i >= 0; i--, j--)
      {
         output[start + count[lsVals[j]] - 1] = input[j];
         count[lsVals[j]]--;
      }

      // Copy the output array to arr[], so that arr[] now
      // contains sorted numbers according to current digit
      for (i = 0, j = start; i < n; i++, j++)
      {
         input[j] = output[j];
         // cout<<"i "<<i<<" "<<output[i]<<" "<<input[i]<<endl;
      }
   }
   // sort by id
   for (int exp = 1; max / exp > 0; exp *= 10)
   {
      // countSort(data, start, n, exp, 2);
      int j, i, count[10];

      for (int i = 0; i < 10; ++i)
      {
         count[i] = 0;
      }

      // Store count of occurrences in count[]
      for (i = 0, j = start; i < n; i++, j++)
      {
         lsVals[j] = (input[j].id / exp) % 10;
         count[lsVals[j]]++;
         // currOutputIndex[i]=outputIndex[i];
      }

      // Change count[i] so that count[i] now contains actual
      //  position of this digit in output[]
      for (i = 1; i < 10; i++)
         count[i] += count[i - 1];

      // Build the output array
      for (i = n - 1, j = start + i; i >= 0; i--, j--)
      {
         output[start + count[lsVals[j]] - 1] = input[j];
         count[lsVals[j]]--;
      }

      // Copy the output array to arr[], so that arr[] now
      // contains sorted numbers according to current digit
      for (i = 0, j = start; i < n; i++, j++)
      {
         input[j] = output[j];
         // cout<<"i "<<i<<" "<<output[i]<<" "<<input[i]<<endl;
      }
   }
}

#pragma acc routine seq
extern "C" void edgeRadixsortDevice(Edge *input, Edge *output, int *lsVals, int start, int n, int max);

// ----------------------------------------
// extern "C" __device__ void mergeTreeNodesInLevel(int nid, Edge *treeNodeEndListArray, int *endListSizes, int *endListCounts, int sumlist, int POLY_TYPE_COUNT){
   
// }


// #pragma acc routine seq
// extern "C" void mergeTreeNodesInLevel(int nid, Edge *treeNodeEndListArray, int *endListSizes, int *endListCounts, int sumlist, int POLY_TYPE_COUNT);

// End list construction: Hybrid version
void SegTree::buildInternalEndListsParallel(int bSize)
{

   int sumlist = sumEndListSize;
   int levelCount = log2(treeSize);
   int sels = sumlist * levelCount;
   int ts = treeSize;
   int leafstart = sels - sumlist;

   #ifdef TIME_L2
   high_resolution_clock::time_point cp1, cp2, cp3, cp4, cp5, cp6;
   cp1 = high_resolution_clock::now();
   #endif

   // cout<<leafstart<<" leaflevel preprocessing "<<ts<<endl;
   // pre-process leaf level
   #pragma acc parallel loop present(endListSizes[0 : ts+1])
   for (int nid =0; nid < ts/2; nid++)
   {
      endListSizes[nid+(ts/2)] = leafstart + endListSizes[nid+(ts/2)];
   }
   // cout<<"leaflevel preprocessing"<<endl;


   // #pragma acc update self(treeNodeEndListArray[0:sels], endListSizes[0:ts], endListCounts[0:ts])
   Edge *output;
   int *lsVals;
   output = (Edge *)malloc((sels) * sizeof(Edge));
   lsVals = (int *)malloc((sels) * sizeof(int));
     
   #ifdef DG
   cout << "POLY_TYPE_COUNT=" << POLY_TYPE_COUNT << endl;
   #endif

   // #pragma acc enter data copyin(this)
   #pragma acc enter data create(output[0 : sels], lsVals[0 : sels])
   #pragma acc parallel loop present(treeNodeEndListArray[0 : sels], endListSizes[0 : ts+1], endListCounts[0 : ts])
   for (int nid =0; nid < +(ts/2); nid++)
   {
      // edgeRadixsort(treeNodeEndListArray, output, lsVals, endListSizes[nid], endListCounts[nid], bSize*2);
      edgeRadixsortDevice(treeNodeEndListArray, output, lsVals, endListSizes[nid+(ts/2)], endListCounts[nid+(ts/2)], bSize * 2);
   }
   // cout<<"Edge radix sort"<<endl;

   // #pragma acc update device(treeNodeEndListArray[0:sels], endListSizes[0:ts+1], endListCounts[0:ts])

   #ifdef TIME_L2
      cp4 = high_resolution_clock::now();
   #endif
   #pragma acc update self(treeNodeEndListArray[0:sels], endListSizes[0:ts+1], endListCounts[0:ts])
   #ifdef TIME_L2
      cp5 = high_resolution_clock::now();
   #endif
   // #pragma acc update self(treeNodeEndListArray[0:sels], endListCounts[0:ts])

   // for(int i=0; i<sels; ++i){
   //    if (treeNodeEndListArray[i].id !=0)
   //       printf("%d->%d ", i, treeNodeEndListArray[i].id);
   // }
   // printf("complete\n");

   #ifdef TIME_L2
   cp3 = high_resolution_clock::now();
   #endif

   const int ptCount=POLY_TYPE_COUNT;

   // Edge *typeFilling;
   // typeFilling=(Edge *)malloc((ptCount)*sizeof(Edge));
   // #pragma enter data private(typeFilling[0:ptCount])
   // #pragma acc data present(treeNodeEndListArray[0:sels], endListSizes[0:ts+1], endListCounts[0:ts])

   // cout<<"Start End list building"<<endl;

   // int printc=0;
   // #pragma acc parallel loop
   for (int lid = 0, ldiv = 4, ldiv2 = 2; lid < levelCount - 1; ++lid, ldiv *= 2, ldiv2 *= 2)
   {
      // for a single level
      int startNode = ts / ldiv;
      int endNode = ts / ldiv2;
      // cout<<bSize<<" ***+++++++building startNode="<<startNode<<" endNode="<<endNode<<endl;

      // // Edge *typeFilling;
      // // typeFilling=(Edge *)malloc((ptCount)*sizeof(Edge));
      // // #pragma enter data private(typeFilling[0:ptCount])
      // #pragma acc enter data copyin(startNode, endNode, ptCount)

      // #pragma acc data present(treeNodeEndListArray[0:sels], endListSizes[0:ts+1], endListCounts[0:ts])
      // #pragma acc parallel loop
      #pragma omp parallel for
      for (int nnid=0; nnid < (endNode-startNode); nnid++)
      {
         int nid=nnid+startNode;
         int nid1 = nid * 2;
         int nid2 = nid1 + 1;
         
         // if(printc<6){
         //    cout<<lid<<" nid1="<<nid1<<" nid2="<<nid2<<endl;
         //    cout<<"s1="<<endListSizes[nid1]<<" s2="<<endListSizes[nid2]<<endl;
         //    printc++;
         // }

         // return;

         int merge = 0;
         int curr1, curr2, filling = 0, currEnd = 0, currid, previd;
         int id, id2;
         Edge nullEdge(-1);

         // Edge *typeFilling;
         // typeFilling = (Edge *)malloc((ptCount) * sizeof(Edge));
         Edge typeFilling[5];

         curr1 = endListSizes[nid1];
         curr2 = endListSizes[nid2];
         // update endlists of this level
         endListSizes[nid] = endListSizes[nid1] - sumlist;
         // initialize filling
         for (int fid = 0; fid < ptCount; ++fid)
         {
            typeFilling[fid] = nullEdge;
         }

         if (endListCounts[nid1] >= 1 && endListCounts[nid2] >= 1)
         {
            // traverse end lists of the child nodes of nid
            if (treeNodeEndListArray[curr1].id <= treeNodeEndListArray[curr2].id)
            {
               currid = treeNodeEndListArray[curr1].id;
               // previd=treeNodeEndListArray[curr1].id;
               previd = currid;
            }
            else
            {
               currid = treeNodeEndListArray[curr2].id;
               // previd=treeNodeEndListArray[curr2].id;
               previd = currid;
            }
            for (id = curr1, id2 = curr2; id < (curr1 + endListCounts[nid1]) && id2 < (curr2 + endListCounts[nid2]);)
            {
               // cout<<id<<" *** previd="<<previd<<" currid="<<currid<<" visiting id="<<treeNodeEndListArray[id].id<<" id2="<<treeNodeEndListArray[id2].id<<endl;
               // set smallest val for this step
               if (treeNodeEndListArray[id].id <= treeNodeEndListArray[id2].id)
               {
                  currid = treeNodeEndListArray[id].id;
               }
               else
               {
                  currid = treeNodeEndListArray[id2].id;
               }
               if (currid != previd)
               {
                  // append values in the filling arrays
                  for (int fid = 0; fid < ptCount; ++fid)
                  {
                     if (typeFilling[fid].id != -1)
                     {
                        treeNodeEndListArray[endListSizes[nid] + currEnd] = typeFilling[fid];
                        currEnd++;
                     }
                  }
                  // reset filling
                  for (int fid = 0; fid < ptCount; ++fid)
                  {
                     typeFilling[fid] = nullEdge;
                  }
                  previd = currid;
               }
               // cout<<id<<" *** previd="<<previd<<" currid="<<currid<<endl;
               if (treeNodeEndListArray[id].id < treeNodeEndListArray[id2].id)
               {
                  typeFilling[treeNodeEndListArray[id].type] = treeNodeEndListArray[id];
                  id++;
                  filling++;
               }
               else if (treeNodeEndListArray[id].id > treeNodeEndListArray[id2].id)
               {
                  typeFilling[treeNodeEndListArray[id2].type] = treeNodeEndListArray[id2];
                  id2++;
                  filling++;
               }
               else
               {
                  typeFilling[treeNodeEndListArray[id].type] = treeNodeEndListArray[id];
                  typeFilling[treeNodeEndListArray[id2].type] = treeNodeEndListArray[id2];
                  id++;
                  id2++;
               }
            }

            // cout<<"comparisons end "<<id<<" "<<(curr1+endListCounts[nid1])<<" "<<id2<<" "<<(curr2+endListCounts[nid2])<<endl;
            // copy remaining elements
            for (int rid = id; rid < (curr1 + endListCounts[nid1]); ++rid)
            {
               currid = treeNodeEndListArray[rid].id;
               if (currid != previd)
               {
                  // append values in the filling arrays
                  for (int fid = 0; fid < ptCount; ++fid)
                  {
                     if (typeFilling[fid].id != -1)
                     {
                        treeNodeEndListArray[endListSizes[nid] + currEnd] = typeFilling[fid];
                        currEnd++;
                     }
                  }
                  // reset filling
                  for (int fid = 0; fid < ptCount; ++fid)
                  {
                     typeFilling[fid] = nullEdge;
                  }
                  previd = currid;
               }
               typeFilling[treeNodeEndListArray[rid].type] = treeNodeEndListArray[rid];
               // treeNodeEndListArray[endListSizes[nid]+currEnd]=treeNodeEndListArray[rid];
               // currEnd++;
            }
            for (int rid2 = id2; rid2 < (curr2 + endListCounts[nid2]); ++rid2)
            {
               currid = treeNodeEndListArray[rid2].id;
               if (currid != previd)
               {
                  // append values in the filling arrays
                  for (int fid = 0; fid < ptCount; ++fid)
                  {
                     if (typeFilling[fid].id != -1)
                     {
                        treeNodeEndListArray[endListSizes[nid] + currEnd] = typeFilling[fid];
                        currEnd++;
                     }
                  }
                  // reset filling
                  for (int fid = 0; fid < ptCount; ++fid)
                  {
                     typeFilling[fid] = nullEdge;
                  }
                  previd = currid;
               }
               typeFilling[treeNodeEndListArray[rid2].type] = treeNodeEndListArray[rid2];
               // treeNodeEndListArray[endListSizes[nid]+currEnd]=treeNodeEndListArray[rid2];
               // currEnd++;
            }
            // append unsaved elelments from last step
            for (int fid = 0; fid < ptCount; ++fid)
            {
               if (typeFilling[fid].id != -1)
               {
                  treeNodeEndListArray[endListSizes[nid] + currEnd] = typeFilling[fid];
                  currEnd++;
               }
            }
         }
         else if (endListCounts[nid1] >= 1)
         {
            for (int id = curr1; id < (curr1 + endListCounts[nid1]); ++id)
            {
               treeNodeEndListArray[endListSizes[nid] + currEnd] = treeNodeEndListArray[id];
               currEnd++;
            }
         }
         else if (endListCounts[nid2] >= 1)
         {
            // cout<<"merge 3"<<endl;
            for (int id = curr2; id < (curr2 + endListCounts[nid2]); ++id)
            {
               treeNodeEndListArray[endListSizes[nid] + currEnd] = treeNodeEndListArray[id];
               currEnd++;
            }
            // cout<<"merge 3> "<<currEnd<<endl;
         }
         endListCounts[nid] = currEnd;
         // cout<<"element count "<<currEnd<<" "<< endListSizes[nid]<<endl;*/
      }
      // #pragma acc wait
      // #pragma omp barrier
      // break;
   }
   sumEndListSize *= levelCount;
   // cout<<"Start End list building complete"<<endl;

   // #pragma acc update self(treeNodeEndListArray[0:sels], endListSizes[0:ts], endListCounts[0:ts])
   // #pragma acc update self(treeNodeEndListArray[0 : sels], endListSizes[0 : ts+1], endListCounts[0 : ts])
   #ifdef TIME_L2
      cp6 = high_resolution_clock::now();
   #endif
   #pragma acc update device(treeNodeEndListArray[0:sels], endListSizes[0:ts+1], endListCounts[0:ts])

   // cout<<"Start End list copied"<<endl;
   #ifdef TIME_L2
      cp2 = high_resolution_clock::now();
   #endif

   #ifdef DT_INFO
      // Calculating total time taken by the program.
      auto dLeafBuild = duration_cast<microseconds>(cp4 - cp1);
      auto dBottomUpBuild = duration_cast<microseconds>(cp6 - cp3);
      auto dCpu2GPUMemcpy = duration_cast<microseconds>(cp2 - cp6);
      auto dGpu2CPUMemcpy = duration_cast<microseconds>(cp5 - cp4);

      cout << " |-endListConstruction: Leaf level build: " << fixed << dLeafBuild.count() << setprecision(10) << endl;
      cout << " |-endListConstruction: GPU to CPU memcpy: " << fixed << dGpu2CPUMemcpy.count() << setprecision(10) << endl;
      cout << " |-endListConstruction: bottom-up-build: " << fixed << dBottomUpBuild.count() << setprecision(10) << endl;
      cout << " |-endListConstruction: CPU to GPU memcpy: " << fixed << dCpu2GPUMemcpy.count() << setprecision(10) << endl;
   #elif TIME_L2
      // Calculating total time taken by the program.
      // auto dLeafBuild = duration_cast<microseconds>(cp4 - cp1);
      // auto dBottomUpBuild = duration_cast<microseconds>(cp6 - cp3);
      auto dCpu2GPUMemcpy = duration_cast<microseconds>(cp2 - cp6);
      auto dGpu2CPUMemcpy = duration_cast<microseconds>(cp5 - cp4);

      cout << fixed << dCpu2GPUMemcpy.count() << setprecision(10) << " " << fixed << dGpu2CPUMemcpy.count() << setprecision(10) << " ";
   #endif
}

// End list construction: GPU only version
void SegTree::buildInternalEndListsGPU(int bSize)
{

   int sumlist = sumEndListSize;
   int levelCount = log2(treeSize);
   int sels = sumlist * levelCount;
   int ts = treeSize;
   int leafstart = sels - sumlist;

   #ifdef TIME_L2
   high_resolution_clock::time_point cp1, cp2, cp3, cp4, cp5, cp6;
   cp1 = high_resolution_clock::now();
   #endif

   // cout<<leafstart<<" leaflevel preprocessing "<<ts<<endl;
   // pre-process leaf level
   #pragma acc parallel loop present(endListSizes[0 : ts+1])
   for (int nid =0; nid < ts/2; nid++)
   {
      endListSizes[nid+(ts/2)] = leafstart + endListSizes[nid+(ts/2)];
   }
   // cout<<"leaflevel preprocessing"<<endl;


   // #pragma acc update self(treeNodeEndListArray[0:sels], endListSizes[0:ts], endListCounts[0:ts])
   Edge *output;
   int *lsVals;
   output = (Edge *)malloc((sels) * sizeof(Edge));
   lsVals = (int *)malloc((sels) * sizeof(int));
     
   #ifdef DG
   cout << "POLY_TYPE_COUNT=" << POLY_TYPE_COUNT << endl;
   #endif

   // #pragma acc enter data copyin(this)
   #pragma acc enter data create(output[0 : sels], lsVals[0 : sels])
   #pragma acc parallel loop present(treeNodeEndListArray[0 : sels], endListSizes[0 : ts+1], endListCounts[0 : ts])
   for (int nid =0; nid < +(ts/2); nid++)
   {
      // edgeRadixsort(treeNodeEndListArray, output, lsVals, endListSizes[nid], endListCounts[nid], bSize*2);
      edgeRadixsortDevice(treeNodeEndListArray, output, lsVals, endListSizes[nid+(ts/2)], endListCounts[nid+(ts/2)], bSize * 2);
   }
   // cout<<"Edge radix sort"<<endl;

   // #pragma acc update device(treeNodeEndListArray[0:sels], endListSizes[0:ts+1], endListCounts[0:ts])

   #ifdef TIME_L2
      cp4 = high_resolution_clock::now();
   #endif
   #pragma acc update self(treeNodeEndListArray[0:sels], endListSizes[0:ts+1], endListCounts[0:ts])
   #ifdef TIME_L2
      cp5 = high_resolution_clock::now();
   #endif
   // #pragma acc update self(treeNodeEndListArray[0:sels], endListCounts[0:ts])

   // for(int i=0; i<sels; ++i){
   //    if (treeNodeEndListArray[i].id !=0)
   //       printf("%d->%d ", i, treeNodeEndListArray[i].id);
   // }
   // printf("complete\n");

   #ifdef TIME_L2
   cp3 = high_resolution_clock::now();
   #endif

   const int ptCount=POLY_TYPE_COUNT;
 
   // Edge *typeFilling;
   // typeFilling=(Edge *)malloc((ptCount)*sizeof(Edge));
   // #pragma enter data private(typeFilling[0:ptCount])
   // #pragma acc data present(treeNodeEndListArray[0:sels], endListSizes[0:ts+1], endListCounts[0:ts])

   // cout<<"Start End list building"<<endl;

   // int printc=0;
   // #pragma acc parallel loop
   for (int lid = 0, ldiv = 4, ldiv2 = 2; lid < levelCount - 1; ++lid, ldiv *= 2, ldiv2 *= 2)
   {
      // for a single level
      int startNode = ts / ldiv;
      int endNode = ts / ldiv2;
      // cout<<bSize<<" ***+++++++building startNode="<<startNode<<" endNode="<<endNode<<endl;

      // // Edge *typeFilling;
      // // typeFilling=(Edge *)malloc((ptCount)*sizeof(Edge));
      // // #pragma enter data private(typeFilling[0:ptCount])
      // #pragma acc enter data copyin(startNode, endNode, ptCount)

      // #pragma acc data present(treeNodeEndListArray[0:sels], endListSizes[0:ts+1], endListCounts[0:ts])
      // #pragma acc parallel loop
      #pragma omp parallel for
      for (int nnid=0; nnid < (endNode-startNode); nnid++)
      {
         int nid=nnid+startNode;
         int nid1 = nid * 2;
         int nid2 = nid1 + 1;
         
         // if(printc<6){
         //    cout<<lid<<" nid1="<<nid1<<" nid2="<<nid2<<endl;
         //    cout<<"s1="<<endListSizes[nid1]<<" s2="<<endListSizes[nid2]<<endl;
         //    printc++;
         // }

         // return;

         int merge = 0;
         int curr1, curr2, filling = 0, currEnd = 0, currid, previd;
         int id, id2;
         Edge nullEdge(-1);

         // Edge *typeFilling;
         // typeFilling = (Edge *)malloc((ptCount) * sizeof(Edge));
         Edge typeFilling[5];

         curr1 = endListSizes[nid1];
         curr2 = endListSizes[nid2];
         // update endlists of this level
         endListSizes[nid] = endListSizes[nid1] - sumlist;
         // initialize filling
         for (int fid = 0; fid < ptCount; ++fid)
         {
            typeFilling[fid] = nullEdge;
         }

         if (endListCounts[nid1] >= 1 && endListCounts[nid2] >= 1)
         {
            // traverse end lists of the child nodes of nid
            if (treeNodeEndListArray[curr1].id <= treeNodeEndListArray[curr2].id)
            {
               currid = treeNodeEndListArray[curr1].id;
               // previd=treeNodeEndListArray[curr1].id;
               previd = currid;
            }
            else
            {
               currid = treeNodeEndListArray[curr2].id;
               // previd=treeNodeEndListArray[curr2].id;
               previd = currid;
            }
            for (id = curr1, id2 = curr2; id < (curr1 + endListCounts[nid1]) && id2 < (curr2 + endListCounts[nid2]);)
            {
               // cout<<id<<" *** previd="<<previd<<" currid="<<currid<<" visiting id="<<treeNodeEndListArray[id].id<<" id2="<<treeNodeEndListArray[id2].id<<endl;
               // set smallest val for this step
               if (treeNodeEndListArray[id].id <= treeNodeEndListArray[id2].id)
               {
                  currid = treeNodeEndListArray[id].id;
               }
               else
               {
                  currid = treeNodeEndListArray[id2].id;
               }
               if (currid != previd)
               {
                  // append values in the filling arrays
                  for (int fid = 0; fid < ptCount; ++fid)
                  {
                     if (typeFilling[fid].id != -1)
                     {
                        treeNodeEndListArray[endListSizes[nid] + currEnd] = typeFilling[fid];
                        currEnd++;
                     }
                  }
                  // reset filling
                  for (int fid = 0; fid < ptCount; ++fid)
                  {
                     typeFilling[fid] = nullEdge;
                  }
                  previd = currid;
               }
               // cout<<id<<" *** previd="<<previd<<" currid="<<currid<<endl;
               if (treeNodeEndListArray[id].id < treeNodeEndListArray[id2].id)
               {
                  typeFilling[treeNodeEndListArray[id].type] = treeNodeEndListArray[id];
                  id++;
                  filling++;
               }
               else if (treeNodeEndListArray[id].id > treeNodeEndListArray[id2].id)
               {
                  typeFilling[treeNodeEndListArray[id2].type] = treeNodeEndListArray[id2];
                  id2++;
                  filling++;
               }
               else
               {
                  typeFilling[treeNodeEndListArray[id].type] = treeNodeEndListArray[id];
                  typeFilling[treeNodeEndListArray[id2].type] = treeNodeEndListArray[id2];
                  id++;
                  id2++;
               }
            }

            // cout<<"comparisons end "<<id<<" "<<(curr1+endListCounts[nid1])<<" "<<id2<<" "<<(curr2+endListCounts[nid2])<<endl;
            // copy remaining elements
            for (int rid = id; rid < (curr1 + endListCounts[nid1]); ++rid)
            {
               currid = treeNodeEndListArray[rid].id;
               if (currid != previd)
               {
                  // append values in the filling arrays
                  for (int fid = 0; fid < ptCount; ++fid)
                  {
                     if (typeFilling[fid].id != -1)
                     {
                        treeNodeEndListArray[endListSizes[nid] + currEnd] = typeFilling[fid];
                        currEnd++;
                     }
                  }
                  // reset filling
                  for (int fid = 0; fid < ptCount; ++fid)
                  {
                     typeFilling[fid] = nullEdge;
                  }
                  previd = currid;
               }
               typeFilling[treeNodeEndListArray[rid].type] = treeNodeEndListArray[rid];
               // treeNodeEndListArray[endListSizes[nid]+currEnd]=treeNodeEndListArray[rid];
               // currEnd++;
            }
            for (int rid2 = id2; rid2 < (curr2 + endListCounts[nid2]); ++rid2)
            {
               currid = treeNodeEndListArray[rid2].id;
               if (currid != previd)
               {
                  // append values in the filling arrays
                  for (int fid = 0; fid < ptCount; ++fid)
                  {
                     if (typeFilling[fid].id != -1)
                     {
                        treeNodeEndListArray[endListSizes[nid] + currEnd] = typeFilling[fid];
                        currEnd++;
                     }
                  }
                  // reset filling
                  for (int fid = 0; fid < ptCount; ++fid)
                  {
                     typeFilling[fid] = nullEdge;
                  }
                  previd = currid;
               }
               typeFilling[treeNodeEndListArray[rid2].type] = treeNodeEndListArray[rid2];
               // treeNodeEndListArray[endListSizes[nid]+currEnd]=treeNodeEndListArray[rid2];
               // currEnd++;
            }
            // append unsaved elelments from last step
            for (int fid = 0; fid < ptCount; ++fid)
            {
               if (typeFilling[fid].id != -1)
               {
                  treeNodeEndListArray[endListSizes[nid] + currEnd] = typeFilling[fid];
                  currEnd++;
               }
            }
         }
         else if (endListCounts[nid1] >= 1)
         {
            for (int id = curr1; id < (curr1 + endListCounts[nid1]); ++id)
            {
               treeNodeEndListArray[endListSizes[nid] + currEnd] = treeNodeEndListArray[id];
               currEnd++;
            }
         }
         else if (endListCounts[nid2] >= 1)
         {
            // cout<<"merge 3"<<endl;
            for (int id = curr2; id < (curr2 + endListCounts[nid2]); ++id)
            {
               treeNodeEndListArray[endListSizes[nid] + currEnd] = treeNodeEndListArray[id];
               currEnd++;
            }
            // cout<<"merge 3> "<<currEnd<<endl;
         }
         endListCounts[nid] = currEnd;
         // cout<<"element count "<<currEnd<<" "<< endListSizes[nid]<<endl;*/
      }
      // #pragma acc wait
      // #pragma omp barrier
      // break;
   }
   sumEndListSize *= levelCount;
   // cout<<"Start End list building complete"<<endl;

   // #pragma acc update self(treeNodeEndListArray[0:sels], endListSizes[0:ts], endListCounts[0:ts])
   // #pragma acc update self(treeNodeEndListArray[0 : sels], endListSizes[0 : ts+1], endListCounts[0 : ts])
   #ifdef TIME_L2
      cp6 = high_resolution_clock::now();
   #endif
   #pragma acc update device(treeNodeEndListArray[0:sels], endListSizes[0:ts+1], endListCounts[0:ts])

   // cout<<"Start End list copied"<<endl;
   #ifdef TIME_L2
      cp2 = high_resolution_clock::now();
   #endif

   #ifdef DT_INFO
      // Calculating total time taken by the program.
      auto dLeafBuild = duration_cast<microseconds>(cp4 - cp1);
      auto dBottomUpBuild = duration_cast<microseconds>(cp6 - cp3);
      auto dCpu2GPUMemcpy = duration_cast<microseconds>(cp2 - cp6);
      auto dGpu2CPUMemcpy = duration_cast<microseconds>(cp5 - cp4);

      cout << " |-endListConstruction: Leaf level build: " << fixed << dLeafBuild.count() << setprecision(10) << endl;
      cout << " |-endListConstruction: GPU to CPU memcpy: " << fixed << dGpu2CPUMemcpy.count() << setprecision(10) << endl;
      cout << " |-endListConstruction: bottom-up-build: " << fixed << dBottomUpBuild.count() << setprecision(10) << endl;
      cout << " |-endListConstruction: CPU to GPU memcpy: " << fixed << dCpu2GPUMemcpy.count() << setprecision(10) << endl;
   #elif TIME_L2
      // Calculating total time taken by the program.
      // auto dLeafBuild = duration_cast<microseconds>(cp4 - cp1);
      // auto dBottomUpBuild = duration_cast<microseconds>(cp6 - cp3);
      auto dCpu2GPUMemcpy = duration_cast<microseconds>(cp2 - cp6);
      auto dGpu2CPUMemcpy = duration_cast<microseconds>(cp5 - cp4);

      cout << fixed << dCpu2GPUMemcpy.count() << setprecision(10) << " " << fixed << dGpu2CPUMemcpy.count() << setprecision(10) << " ";
   #endif
}


void SegTree::copyLeavesFromVector()
{
   int j;
   int sumlist = sumEndListSize;
   int ts = treeSize;
   int levelCount = log2(treeSize);
   int sels = sumlist * levelCount;

   int leafstart = sels - sumEndListSize;

   treeNodeEndListArray = (Edge *)malloc((sels) * sizeof(Edge));

   // int leafstart=sels-sumlist;
   for (int i = treeSize / 2; i < treeSize; i++)
   {
      // for(int i=0; i<treeSize; i++){
      j = 0;
      for (Edge e : *treeNode->at(i)->endEdgeSet)
      {
         treeNodeEndListArray[leafstart + endListSizes[i] + j] = e;
         ++j;
      }
      endListCounts[i] = j;
   }
   // #pragma acc enter data copyin(treeNodeEndListArray[0:sels], endListCounts[0:ts])
}

#if GPU!=0

void SegTree ::insertAllToEdgeListParallel(int bSize)
{
   int *left, *right;
   LONG *intBoaders;
   REAL *boaders;
   int epc = elementaryPointCount;

   // int epc=m_elementaryPoints->size();
   // elementaryPointCount=epc;

   // cout<<"epc="<<epc<<endl;
   boaders = (REAL *)malloc((epc) * sizeof(REAL));
   intBoaders = (LONG *)malloc((epc) * sizeof(LONG));
   left = (int *)malloc((epc) * sizeof(int));
   right = (int *)malloc((epc) * sizeof(int));

   // #pragma acc enter data copyin(this)
   #pragma acc enter data create(boaders[0 : epc], intBoaders[0 : epc], left[0 : epc], right[0 : epc])

   insertToLeafLevelEdgeListParallel(boaders, intBoaders, left, right, bSize);
   // cout<<"done part 1"<<endl;
   // buildInternalEndListsParallel(bSize);
   buildInternalEndListsGPU(bSize);
   // cout<<"done part 2"<<endl;

}
#endif
// ------------------------------------------------------------------

/*
multicore end list build
*/

void SegTree ::buildInternalEndListsMulticore(int bSize)
{

   int sumlist = sumEndListSize;
   int levelCount = log2(treeSize);
   int sels = sumlist * levelCount;
   int ts = treeSize;
   int leafstart = sels - sumlist;


   // pre-process leaf level
   #pragma omp parallel for
   for (int nid = treeSize / 2; nid < treeSize; nid++)
   {
      endListSizes[nid] = leafstart + endListSizes[nid];
   }
   // cout<<"before edgeRadixSort"<<endl;
   // cout<<"sels="<<sels<<" sumlist="<<sumlist<<endl;

   #ifdef DG
   high_resolution_clock::time_point cp1, cp2, cp3;
   cp1 = high_resolution_clock::now();
   #endif

   // #pragma acc update self(treeNodeEndListArray[0:sels], endListSizes[0:ts], endListCounts[0:ts])
   Edge *output;
   int *lsVals;
   output = (Edge *)malloc((sels) * sizeof(Edge));
   lsVals = (int *)malloc((sels) * sizeof(int));
  
   // cout<<"sels="<<sels<<" sumlist="<<sumlist<<endl;

   #ifdef DG
   cout << "POLY_TYPE_COUNT=" << POLY_TYPE_COUNT << endl;
   #endif
   

   #pragma omp parallel for
   // #pragma omp parallel for schedule(dynamic, 20000)
   // #pragma omp parallel for schedule(auto)
   for (int nid = treeSize / 2; nid < treeSize; nid++)
   {
      edgeRadixsort(treeNodeEndListArray, output, lsVals, endListSizes[nid], endListCounts[nid], bSize * 2);
   }

   #ifdef DG 
      cp3 = high_resolution_clock::now();
   #endif

   int nid;
   // #pragma acc parallel loop
   for (int lid = 0, ldiv = 4, ldiv2 = 2; lid < levelCount - 1; ++lid, ldiv *= 2, ldiv2 *= 2)
   {
      // for a single level
      int startNode = treeSize / ldiv;
      int endNode = treeSize / ldiv2;
   // cout<<bSize<<" ***+++++++building startNode="<<startNode<<" endNode="<<endNode<<endl;

      // #pragma omp parallel for schedule(auto)
      #pragma omp parallel for
      for (nid = startNode; nid < endNode; nid++)
      {
         int nid1 = nid * 2;
         int nid2 = nid1 + 1;
         int merge = 0;
         int curr1, curr2, filling = 0, currEnd = 0, currid, previd;
         int id, id2;
         Edge nullEdge(-1);

         Edge *typeFilling;
         typeFilling = (Edge *)malloc((POLY_TYPE_COUNT) * sizeof(Edge));

         curr1 = endListSizes[nid1];
         curr2 = endListSizes[nid2];
         // update endlists of this level
         endListSizes[nid] = endListSizes[nid1] - sumlist;
         // initialize filling
         for (int fid = 0; fid < POLY_TYPE_COUNT; ++fid)
         {
            typeFilling[fid] = nullEdge;
         }

         if (endListCounts[nid1] >= 1 && endListCounts[nid2] >= 1)
         {
            // traverse end lists of the child nodes of nid
            if (treeNodeEndListArray[curr1].id <= treeNodeEndListArray[curr2].id)
            {
               currid = treeNodeEndListArray[curr1].id;
               // previd=treeNodeEndListArray[curr1].id;
               previd = currid;
            }
            else
            {
               currid = treeNodeEndListArray[curr2].id;
               // previd=treeNodeEndListArray[curr2].id;
               previd = currid;
            }
            for (id = curr1, id2 = curr2; id < (curr1 + endListCounts[nid1]) && id2 < (curr2 + endListCounts[nid2]);)
            {
               // cout<<id<<" *** previd="<<previd<<" currid="<<currid<<" visiting id="<<treeNodeEndListArray[id].id<<" id2="<<treeNodeEndListArray[id2].id<<endl;
               // set smallest val for this step
               if (treeNodeEndListArray[id].id <= treeNodeEndListArray[id2].id)
               {
                  currid = treeNodeEndListArray[id].id;
               }
               else
               {
                  currid = treeNodeEndListArray[id2].id;
               }
               if (currid != previd)
               {
                  // append values in the filling arrays
                  for (int fid = 0; fid < POLY_TYPE_COUNT; ++fid)
                  {
                     if (typeFilling[fid].id != -1)
                     {
                        treeNodeEndListArray[endListSizes[nid] + currEnd] = typeFilling[fid];
                        currEnd++;
                     }
                  }
                  // reset filling
                  for (int fid = 0; fid < POLY_TYPE_COUNT; ++fid)
                  {
                     typeFilling[fid] = nullEdge;
                  }
                  previd = currid;
               }
               // cout<<id<<" *** previd="<<previd<<" currid="<<currid<<endl;
               if (treeNodeEndListArray[id].id < treeNodeEndListArray[id2].id)
               {
                  typeFilling[treeNodeEndListArray[id].type] = treeNodeEndListArray[id];
                  id++;
                  filling++;
               }
               else if (treeNodeEndListArray[id].id > treeNodeEndListArray[id2].id)
               {
                  typeFilling[treeNodeEndListArray[id2].type] = treeNodeEndListArray[id2];
                  id2++;
                  filling++;
               }
               else
               {
                  typeFilling[treeNodeEndListArray[id].type] = treeNodeEndListArray[id];
                  typeFilling[treeNodeEndListArray[id2].type] = treeNodeEndListArray[id2];
                  id++;
                  id2++;
               }
            }

            // cout<<"comparisons end "<<id<<" "<<(curr1+endListCounts[nid1])<<" "<<id2<<" "<<(curr2+endListCounts[nid2])<<endl;
            // copy remaining elements
            for (int rid = id; rid < (curr1 + endListCounts[nid1]); ++rid)
            {
               currid = treeNodeEndListArray[rid].id;
               if (currid != previd)
               {
                  // append values in the filling arrays
                  for (int fid = 0; fid < POLY_TYPE_COUNT; ++fid)
                  {
                     if (typeFilling[fid].id != -1)
                     {
                        treeNodeEndListArray[endListSizes[nid] + currEnd] = typeFilling[fid];
                        currEnd++;
                     }
                  }
                  // reset filling
                  for (int fid = 0; fid < POLY_TYPE_COUNT; ++fid)
                  {
                     typeFilling[fid] = nullEdge;
                  }
                  previd = currid;
               }
               typeFilling[treeNodeEndListArray[rid].type] = treeNodeEndListArray[rid];
               // treeNodeEndListArray[endListSizes[nid]+currEnd]=treeNodeEndListArray[rid];
               // currEnd++;
            }
            for (int rid2 = id2; rid2 < (curr2 + endListCounts[nid2]); ++rid2)
            {
               currid = treeNodeEndListArray[rid2].id;
               if (currid != previd)
               {
                  // append values in the filling arrays
                  for (int fid = 0; fid < POLY_TYPE_COUNT; ++fid)
                  {
                     if (typeFilling[fid].id != -1)
                     {
                        treeNodeEndListArray[endListSizes[nid] + currEnd] = typeFilling[fid];
                        currEnd++;
                     }
                  }
                  // reset filling
                  for (int fid = 0; fid < POLY_TYPE_COUNT; ++fid)
                  {
                     typeFilling[fid] = nullEdge;
                  }
                  previd = currid;
               }
               typeFilling[treeNodeEndListArray[rid2].type] = treeNodeEndListArray[rid2];
               // treeNodeEndListArray[endListSizes[nid]+currEnd]=treeNodeEndListArray[rid2];
               // currEnd++;
            }
            // append unsaved elelments from last step
            for (int fid = 0; fid < POLY_TYPE_COUNT; ++fid)
            {
               if (typeFilling[fid].id != -1)
               {
                  treeNodeEndListArray[endListSizes[nid] + currEnd] = typeFilling[fid];
                  currEnd++;
               }
            }
         }
         else if (endListCounts[nid1] >= 1)
         {
            for (int id = curr1; id < (curr1 + endListCounts[nid1]); ++id)
            {
               treeNodeEndListArray[endListSizes[nid] + currEnd] = treeNodeEndListArray[id];
               currEnd++;
            }
         }
         else if (endListCounts[nid2] >= 1)
         {
            // cout<<"merge 3"<<endl;
            for (int id = curr2; id < (curr2 + endListCounts[nid2]); ++id)
            {
               treeNodeEndListArray[endListSizes[nid] + currEnd] = treeNodeEndListArray[id];
               currEnd++;
            }
            // cout<<"merge 3> "<<currEnd<<endl;
         }
         endListCounts[nid] = currEnd;
         // cout<<"element count "<<currEnd<<" "<< endListSizes[nid]<<endl;
      }
      // #pragma omp barrier
      // break;
   }
   sumEndListSize *= levelCount;

   free(output);
   free(lsVals);

   #ifdef DG 
   cp2 = high_resolution_clock::now();
   #endif
   #ifdef DG
      // Calculating total time taken by the program.
      auto d1 = duration_cast<microseconds>(cp2 - cp1);
      auto d2 = duration_cast<microseconds>(cp2 - cp3);

      cout << "endListConstruction: buildAll() " << fixed << d1.count() << setprecision(10) << endl;
      cout << "endListConstruction: build() " << fixed << d2.count() << setprecision(10) << endl;
   #endif
}

void SegTree::preparePointsToNodeMapMulticore(REAL *boaders, LONG *intBoaders, int *left, int *right)
{

   int firstLeafIndex = treeSize / 2;
   int leafIndex = firstLeafIndex;
   int i;
   // int ic=elementaryArrayCount;
   int osize = intervalCount;
   int epc = elementaryPointCount;
   int icupper = intervalCountUpperBound;

   REAL absMinEndPoint = abs(minInterval);
   int maxDigitCount = getDigitCount(maxInterval + absMinEndPoint);
// -----------------
#pragma omp parallel for
   for (i = 1; i < epc; i++)
   {
      int leftLeaf = leafIndex + i - 1;
      int rightLeaf = leafIndex + i;
      // handling point intervals
      // ---------------------------------------------
      if (i != 1 && elementaryArray[i] == elementaryArray[i + 1])
      {
         rightLeaf = leafIndex + i + 1;
      }
      // ---------------------------------------------
      left[i - 1] = leftLeaf;
      right[i - 1] = rightLeaf;
      boaders[i - 1] = elementaryArray[i];
      if (elementaryArray[0] < 0)
      {
         intBoaders[i - 1] = realToInt(elementaryArray[i] + absMinEndPoint, (MAX_INTERVAL_DIGIT_COUNT - maxDigitCount));
      }
      else
      {
         intBoaders[i - 1] = realToInt(elementaryArray[i], (MAX_INTERVAL_DIGIT_COUNT - maxDigitCount));
      }
   }
}

void SegTree::insertToLeafLevelEdgeListMulticore(Edge *edges, REAL *boaders, LONG *intBoaders, int *left, int *right, int bSize)
{
   int epc = elementaryPointCount;
   int *boadersSortedIndex;
   boadersSortedIndex = (int *)malloc((epc) * sizeof(int));

   for (int i = 0; i < epc; ++i)
      boadersSortedIndex[i] = i;

   preparePointsToNodeMapMulticore(boaders, intBoaders, left, right);
   radixsort(intBoaders, boadersSortedIndex, 0, epc, epc, MAX_INTERVAL_DIGIT_COUNT);

   int ic = intervalCount, ts = treeSize, ts2 = ts * 2;
   // #pragma acc update self(boaders[0:epc], left[0:epc], right[0:epc])
   int *enids, *stids;
   enids = (int *)malloc((ic) * sizeof(int));
   stids = (int *)malloc((ic) * sizeof(int));

   #pragma omp parallel for
   for (int eid = 0; eid < ic; ++eid)
   {
      enids[eid] = findBoarder(boaders, boadersSortedIndex, epc, edges[eid].end);
      stids[eid] = findBoarder(boaders, boadersSortedIndex, epc, edges[eid].start);
   }

   // use a linked list and copy**************************************

   #pragma omp parallel for
   for (int eid = 0; eid < ic; ++eid)
   {
      int after, before;

      // end point
      if (enids[eid] >= 0)
      {
         after = right[boadersSortedIndex[enids[eid]]];
         // after=right[enids[eid]];
         if (nodeIntervalArray[after * 2] != nodeIntervalArray[after * 2 + 1])
         {
            #pragma omp atomic update
            endListSizes[after] = endListSizes[after] + 1;
         }
      }
      // start point
      if (stids[eid] >= 0)
      {
         before = left[boadersSortedIndex[stids[eid]]];
         // before=left[stids[eid]];
         if (nodeIntervalArray[before * 2] != nodeIntervalArray[before * 2 + 1])
         {
            #pragma omp atomic update
            endListSizes[before] = endListSizes[before] + 1;
         }
      }
   }

   thrust::exclusive_scan(endListSizes, endListSizes + treeSize + 1, endListSizes);

   sumEndListSize = endListSizes[treeSize];
   // cout<<"Assigned sumEndListSize="<<sumEndListSize<<endl;
   int sels = sumEndListSize * log2(treeSize);
   treeNodeEndListArray = (Edge *)malloc((sels) * sizeof(Edge));
   int leafstart = sels - sumEndListSize;

   #pragma omp parallel for
   for (int eid = 0; eid < ic; ++eid)
   {
      int after, before, localid;

      // end point
      if (enids[eid] >= 0)
      {
         after = right[boadersSortedIndex[enids[eid]]];
         // after=right[enids[eid]];
         if (nodeIntervalArray[after * 2] != nodeIntervalArray[after * 2 + 1])
         {
            #pragma omp critical
            {
               // treeNodeEndListArray[leafstart+endListSizes[after]+(endListCounts[after]+=1)]=edges[eid];
               treeNodeEndListArray[leafstart + endListSizes[after] + (endListCounts[after])] = edges[eid];
               endListCounts[after] += 1;
            }
         }
      }
      // start point
      if (stids[eid] >= 0)
      {
         before = left[boadersSortedIndex[stids[eid]]];
         // before=left[stids[eid]];
         if (nodeIntervalArray[before * 2] != nodeIntervalArray[before * 2 + 1])
         {
            #pragma omp critical
            {
               // treeNodeEndListArray[leafstart+endListSizes[before]+(endListCounts[before]+=1)]=edges[eid];
               treeNodeEndListArray[leafstart + endListSizes[before] + (endListCounts[before])] = edges[eid];
               endListCounts[before] += 1;
            }
         }
      }
   }
   free(boadersSortedIndex);
   free(enids);
   free(stids);
}

void SegTree ::insertAllToEdgeListMulticore(Edge *edges, int bSize)
{
   int *left, *right;
   LONG *intBoaders;
   REAL *boaders;
   int epc = elementaryPointCount;

   boaders = (REAL *)malloc((epc) * sizeof(REAL));
   intBoaders = (LONG *)malloc((epc) * sizeof(LONG));
   left = (int *)malloc((epc) * sizeof(int));
   right = (int *)malloc((epc) * sizeof(int));

   insertToLeafLevelEdgeListMulticore(edges, boaders, intBoaders, left, right, bSize);
   
   free(boaders);
   free(intBoaders);
   free(left);
   free(right);
   // cout<<"done part 1"<<endl;
   buildInternalEndListsMulticore(bSize);
}

// ------------------------------------------------------------------
// extern "C" __device__ bool isIntersect(REAL start, REAL end, REAL es, REAL ee){
bool isIntersect(REAL start, REAL end, REAL es, REAL ee)
{
   // bool isIntersect(REAL start, REAL end, REAL es, REAL ee){
   // returns true if this edge contains arg. e or vice versa
   if ((start <= es && end >= ee) || (es <= start && ee >= end))
      return true;
   // partial overlap
   if ((start <= es && end >= es) || (ee <= end && ee >= start))
      return true;
   else
      return false;
}

extern "C" __device__ bool doesContain(REAL start, REAL end, REAL es, REAL ee){
// bool doesContain(REAL start, REAL end, REAL es, REAL ee){
   // bool isContain(REAL start, REAL end, REAL es, REAL ee){
   // if(start <= e.start && end >= e.end)
   // if(start <= es && end >= ee)
   if (start <= es && end >= ee)
      return true;
   // else if(start == e.start && end == e.end)
   else
      return false;
}
#pragma acc routine seq
extern "C" bool doesContain(REAL start, REAL end, REAL es, REAL ee);

// edge and node
extern "C" __device__ bool doesContainEndPoint(REAL start, REAL end, REAL es, REAL ee){  //********* check this
   if(es!=ee && ((start>=es && start<=ee) || (end>=es && end<=ee)))
      return true;
   else
      return false;
}
#pragma acc routine seq
extern "C" bool doesContainEndPoint(REAL start, REAL end, REAL es, REAL ee);
/*
#pragma acc routine seq
extern "C" bool isIntersect(REAL start, REAL end, REAL es, REAL ee);
#pragma acc routine seq
extern "C" bool isContain(REAL start, REAL end, REAL es, REAL ee);
*/
/*
extern "C" __device__ void insertEdgeCountParallel(Edge e, unsigned char *treeMap, REAL *nodeIntervalArray, int ts, int levelCount, int *coverListSizes2){
// extern "C" __device__ void insertEdgeCountParallel(Edge e, REAL *nodeIntervalArray, int ts, int levelCount, int *coverListSizes2){
// void insertEdgeCountParallel(Edge e, REAL *nodeIntervalArray, int ts, int levelCount, int *coverListSizes2){
   int ld=2, levelpass=0, go=1;
   // unsigned char *treeMap;

   int bitArraySize=(ts+8)/8;
   // cout<<"ts="<<ts<<" bitArraySize="<<bitArraySize<<endl;
   // treeMap=(unsigned char *)malloc((bitArraySize)*sizeof(unsigned char));
   // #pragma acc update self(nodeIntervalArray[0:ts*2])

   // set treemap
   for(int i=0; i<bitArraySize; ++i)
      treeMap[i]=0;
   // cout<<"allocated"<<endl;

   // #pragma acc enter data create(coverListSizes2[0:ts])
   // #pragma acc parallel loop present(nodeIntervalArray[0:ts*2])
   for(int lid=1; lid<levelCount; lid++){
      // if(lid>1 && levelpass!=1){
      //    lid=levelCount;
      //    go=0;
      // }
      // #pragma acc loop vector
      for(int nid=ld; nid<2*ld && go==1; ++nid){
      // for(int nid=ld; nid<2*ld; ++nid){
         int nid2=nid*2;
         // cout<<"--- nid="<<nid<<" lid="<<lid<<endl;
         // printf("nid=%d lid=%d\n", nid, lid);

         // if(treeMap[nid]==0){
         int val=(treeMap[nid/8]&(1<<(nid%8)))!=0;
         if(val==0){
            if(isContain(e.start, e.end, nodeIntervalArray[nid2], nodeIntervalArray[nid2+1])){
               // cout<<"[[[ "<<e.start<<" "<<e.end<<" -- "<<nodeIntervalArray[nid2]<<" "<<nodeIntervalArray[nid2+1]<<" /// "<<isIntersect(e.start, e.end, nodeIntervalArray[nid2], nodeIntervalArray[nid2+1])<<" +++ "<<isContain(e.start, e.end, nodeIntervalArray[nid2], nodeIntervalArray[nid2+1])<<endl;
               // cout<<"+++ "<<e.start<<" "<<e.end<<" -- "<<nodeIntervalArray[nid2]<<" "<<nodeIntervalArray[nid2+1]<<" "<<nid2<<endl;
               #pragma acc atomic
               coverListSizes2[nid]+=1;
               // printf("nid=%d size=%d\n", nid, coverListSizes2[nid]);
               if(lid<levelCount-1){
                  // treeMap[nid2]=1;   // mark not to visit childeren
                  // treeMap[nid2+1]=1;   // mark not to visit childeren
                  treeMap[nid2/8] |= (1 << (nid2%8));
                  treeMap[(nid2+1)/8] |= (1 << ((nid2+1)%8));
               }

            // }else if(isIntersect(e.start, e.end, nodeIntervalArray[nid2], nodeIntervalArray[nid2+1])){
            //    levelpass=1;
            // }
            }else if(lid<levelCount-1 && !isIntersect(e.start, e.end, nodeIntervalArray[nid2], nodeIntervalArray[nid2+1])){
               // treeMap[nid2]=1;   // mark not to visit childeren
               // treeMap[nid2+1]=1;   // mark not to visit childeren
               treeMap[nid2/8] |= (1 << (nid2%8));
               treeMap[(nid2+1)/8] |= (1 << ((nid2+1)%8));
            }
         }else if(lid<levelCount-1){
            // treeMap[nid2]=1;
            // treeMap[nid2+1]=1;
            treeMap[nid2/8] |= (1 << (nid2%8));
            treeMap[(nid2+1)/8] |= (1 << ((nid2+1)%8));
         }
      }
      // #pragma acc wait
      ld*=2;
   }
   // cout<<"end search"<<endl;
   // free(treeMap);
}


#pragma acc routine seq
extern "C" void insertEdgeCountParallel(Edge e, unsigned char *treeMap, REAL *nodeIntervalArray, int ts, int levelCount, int *coverListSizes2);
// #pragma acc routine seq
// extern "C" void insertEdgeCountParallel(Edge e, REAL *nodeIntervalArray, int ts, int levelCount, int *coverListSizes2);
// #pragma acc routine seq
// extern "C" void insertEdgeCountParallel(Edge e, REAL *nodeIntervalArray, int ts, int levelCount, int *coverListSizes2, Edge *treeNodeCoverListArray);
*/

/*
void SegTree::insertAllParallel2(Edge *edges){
   int ic=intervalCount, ts=treeSize;
   int lc=log2(ts);
   int bitArraySize=(ts+8)/8;

   unsigned char *treeMap;
   treeMap=(unsigned char *)malloc((bitArraySize)*sizeof(unsigned char));
   // cout<<"111"<<endl;
   // #pragma acc update self(nodeIntervalArray[0:ts*2])
   // cout <<"test11"<< endl;
   // #pragma acc update self(nodeIntervalArray[0:ts*2])


   // #pragma acc data present(nodeIntervalArray[0:ts*2])
   // #pragma acc enter data create(coverListSizes2[0:ts], treeMap[0:ts])
   // #pragma acc parallel loop present(nodeIntervalArray[0:ts*2]) private(treeMap[0:ts])

   #pragma acc enter data create(coverListSizes2[0:ts])
   #pragma acc parallel loop present(nodeIntervalArray[0:ts*2]) private(treeMap[0:bitArraySize])
   for(int eid=0; eid<ic; ++eid){
      if(isContain(edges[eid].start, edges[eid].end, nodeIntervalArray[2], nodeIntervalArray[3])){
         coverListSizes2[1]+=1;
      }else{
         // cout<<"&&&&&&** "<<edges[eid].start<<" "<<edges[eid].end<<endl;
         insertEdgeCountParallel(edges[eid], treeMap, nodeIntervalArray, ts, lc, coverListSizes2);
         // insertEdgeCountParallel(edges[eid], nodeIntervalArray, ts, lc, coverListSizes2);
         // break;
         // cout<<"** "<<edges[eid].start<<endl;
      }
   }

   // #pragma acc update self(coverListSizes2[0:ts])
   // thrust::exclusive_scan(coverListSizes2, coverListSizes2+ts, coverListSizes2);
   // #pragma acc update device(coverListSizes2[0:ts])

   // sumCoverListSize=coverListSizes2[ts];
   // treeNodeCoverListArray=(Edge *)malloc((sumCoverListSize)*sizeof(Edge));

   // #pragma acc enter data create(treeNodeCoverListArray[0:ts])
   // #pragma acc parallel loop present(nodeIntervalArray[0:ts*2], coverListSizes2[0:ts])
   // for(int eid=0; eid<ic; ++eid){
   //    if(isContain(edges[eid].start, edges[eid].end, nodeIntervalArray[2], nodeIntervalArray[3])){
   //       #pragma acc atomic
   //       treeNodeCoverListArray[coverListSizes2[1]+(coverListCounts[1]+=1)]=edges[eid];
   //    }else{
   //       // insertEdgeParallel(edges[eid], nodeIntervalArray, ts, lc, coverListSizes2, treeNodeCoverListArray);

   //       int ld=2, levelpass=0, go=1;
   //       Edge e=edges[eid];

   //       for(int lid=1; lid<lc-1; lid++){
   //          if(levelpass!=1){
   //             lid=lc;
   //             go=0;
   //          }
   //          // #pragma acc loop vector
   //          for(int nid=ld; nid<2*ld && go==1; ++nid){
   //             int nid2=nid*2;
   //             // cout<<"--- nid="<<nid<<" lid="<<lid<<endl;
   //             // cout<<"[[[ "<<e.start<<" "<<e.end<<" -- "<<nodeIntervalArray[nid2]<<" "<<nodeIntervalArray[nid2+1]<<" /// "<<isIntersect(e.start, e.end, nodeIntervalArray[nid2], nodeIntervalArray[nid2+1])<<" +++ "<<isContain(e.start, e.end, nodeIntervalArray[nid2], nodeIntervalArray[nid2+1])<<" "<<treeMap[nid]<<endl;

   //             // if(treeMap[nid]==0){
   //                if(isContain(e.start, e.end, nodeIntervalArray[nid2], nodeIntervalArray[nid2+1])){
   //                   // #pragma acc atomic
   //                   treeNodeCoverListArray[coverListSizes2[nid]+coverListCounts[nid]]=e;
   //                   coverListCounts[nid]+=1;
   //                   // coverListSizes2[nid]+=1;
   //                }else if(isIntersect(e.start, e.end, nodeIntervalArray[nid2], nodeIntervalArray[nid2+1])){
   //                   levelpass=1;
   //                }
   //          }
   //          // #pragma acc wait
   //          ld*=2;
   //       }

   //    }
   // }


//    cout<<"1114444444"<<endl;

//   cout<<"new "<<endl;
//    for(int i=0; i<treeSize; ++i){
//       cout<<i<<">"<<coverListSizes2[i]<<", ";
//       coverListSizes2[i]=0;    //reset
//    }
//   cout<<endl;
}
*/
// ------------------------------------------------------------------

/*void SegTree ::insertEdge(Edge x, Node *root, int nodeIdx){
   // cout<<"N# "<<nodeIdx<<endl;
   if (x.contains(root->interval)){
      // if (x.type==2) cout<<"inserting id="<<x.id<<" type="<<x.type<<" start="<<x.start<<" end="<<x.end<<endl;
      // cout<<" Root  "<<nodeIdx<<" "<<endl;
      root->addToCoverList(x);
   }else{
      // if (x.type==2) cout<<"inserting id="<<x.id<<" type="<<x.type<<" start="<<x.start<<" end="<<x.end<<endl;
      int leftIdx=2*nodeIdx;
      if (leftIdx<treeSize){
         Node *leftChild=treeNode->at(leftIdx);
         if (x.intersects(leftChild->interval)){
            insertEdge(x, leftChild, leftIdx);
         }
      }
      int rightIdx=2*nodeIdx+1;
      if (rightIdx < treeSize){
         Node *rightChild = treeNode->at(rightIdx);
         if (x.intersects(rightChild->interval)){
            insertEdge(x, rightChild, rightIdx);
         }
      }
   }
}

void SegTree ::insert(Edge x){
   insertEdge(x, treeNode->at(1), 1);
}

void SegTree ::insertAll(Edge *edges){
   // int counter = 0;
   for (int eid=0; eid<intervalCount; ++eid){
      insert(edges[eid]);
      // cout<<"inserted "<<counter<<" "<<endl;
      // counter++;
   }
   // cout <<"test"<< endl;
   // insert(Edge(1,4));
   //  insert(Edge(1,2));
}*/

// -------------------------------------
// void SegTree ::insertEdgeParallel(Edge x, Node *root, int nodeIdx, omp_lock_t &lock){

//    // cout<<"N# "<<nodeIdx<<endl;
//    // #pragma omp taskwait
//    if (x.contains(root->interval)){
//       // if (x.type==2) cout<<"inserting id="<<x.id<<" type="<<x.type<<" start="<<x.start<<" end="<<x.end<<endl;
//       // cout<<" Root  "<<nodeIdx<<" "<<endl;
//       omp_set_lock(&lock);
//       root->addToCoverList(x);
//       omp_unset_lock(&lock);
//       return;
//    }else{
//       // if (x.type==2) cout<<"inserting id="<<x.id<<" type="<<x.type<<" start="<<x.start<<" end="<<x.end<<endl;
//       int leftIdx=2*nodeIdx;
//       if (leftIdx<treeSize){
//          Node *leftChild=treeNode->at(leftIdx);
//          if (x.intersects(leftChild->interval)){
//             // #pragma omp task private(x, leftIdx) firstprivate(leftChild) shared(lock)
//             // #pragma omp task firstprivate(root) private(x, leftIdx) shared(lock)
//             insertEdgeParallel(x, leftChild, leftIdx, lock);
//          }
//       }
//       int rightIdx=2*nodeIdx+1;
//       if (rightIdx < treeSize){
//          Node *rightChild = treeNode->at(rightIdx);
//          if (x.intersects(rightChild->interval)){
//             // #pragma omp task private(x, rightIdx) firstprivate(rightChild) shared(lock)
//             // #pragma omp task firstprivate(root) private(x, rightIdx) shared(lock)
//             insertEdgeParallel(x, rightChild, rightIdx, lock);
//          }
//       }
//    }
// }

// void SegTree ::insertParallel(Edge x, omp_lock_t &lock){
//    // insertEdge(x, treeNode->at(1), 1);
//    // #pragma omp parallel
//    // #pragma omp task
//    {
//       // #pragma omp single nowait
//       // #pragma omp single
//       insertEdgeParallel(x, treeNode->at(1), 1, lock);
//    }
//    // int ts=treeSize;

//    // for(int nid=0; nid<ts; ++nid){
//    //    if(x.contains(treeNode->at(nid)->interval)){
//    //       omp_set_lock(&lock);
//    //       root->addToCoverList(x);
//    //       omp_unset_lock(&lock);
//    //    }

//    // }
//    // Stack to store the nodes
//    /*stack<pair(Node, Node)> nodes;
//    // push the current node onto the stack
//    nodes.push(treeNode->at(1));
//    // loop while the stack is not empty
//    while (!nodes.empty())
//       // store the current node and pop it from the stack
//       Node* curr = nodes.top();
//       nodes.pop();
//       // current node has been travarsed
//    if(curr != NULL){
//       cout << curr->key << " ";
//       if(x.contains(treeNode->at(nid)->interval)){
//          omp_set_lock(&lock);
//          root->addToCoverList(x);
//          omp_unset_lock(&lock);
//       }
//       // store all the childrent of current node from
//       // right to left.
//       vector<Node*>::iterator it = curr->child.end();

//       while (it != curr->child.begin()) {
//          it--;
//          nodes.push(*it);
//       }
//    }*/
// }

// void SegTree ::insertAllParallel(Edge *edges){
//    // int counter = 0;
//    omp_lock_t lock;
//    omp_init_lock(&lock);

//    // omp_set_dynamic(0);
//    // omp_set_num_threads(8);
//    // #pragma omp parallel for num_threads(8)
//    #pragma omp parallel for
//    for (int eid=0; eid<intervalCount; ++eid){
//       insertParallel(edges[eid], lock);
//       // cout<<"inserted "<<counter<<" "<<endl;
//       // counter++;
//    }
//    omp_destroy_lock(&lock);
//    // cout <<"test"<< endl;
//    // insert(Edge(1,4));
//    //  insert(Edge(1,2));
// }

// ----------------------------------------------------------
void SegTree ::insertEdgeParallel(Edge x, int nodeIdx)
{

   // cout<<"N# "<<nodeIdx<<endl;
   // #pragma omp taskwait
   // if (x.contains(root->interval)){
   if (x.contains(nodeIntervalArray[nodeIdx * 2], nodeIntervalArray[nodeIdx * 2 + 1]))
   // if(doesContain(x.start, x.end, nodeIntervalArray[nodeIdx*2], nodeIntervalArray[nodeIdx*2+1]))
   {
      // --------------------------------- test ----------------------
      /*int parent=nodeIdx;
      if(nodeIdx%2!=0){ //right child
         parent=nodeIdx-1;
         if(x.start<nodeIntervalArray[parent] || x.start>nodeIntervalArray[nodeIdx*2]){
            cout<<x.start<< " Mis match!! "<<nodeIntervalArray[parent]<<" "<<nodeIntervalArray[nodeIdx*2]<<endl;
         }
      }else{ //left child
         if(x.end>nodeIntervalArray[parent+1] || x.end<nodeIntervalArray[nodeIdx*2+1]){
         cout<<x.end<< " Mis match!! "<<nodeIntervalArray[parent]<<" "<<nodeIntervalArray[nodeIdx*2]<<endl;
         }
      }
      if(nodeIdx==19 || nodeIdx==15){
         cout<<"19-> "<<x.start<<" "<<x.end<<" parent="<<parent<<" "<<nodeIntervalArray[parent]<<" "<<nodeIntervalArray[nodeIdx*2]<<endl;
      }*/
      // --------------------------------- test ----------------------

      // if (x.type==2) cout<<"inserting id="<<x.id<<" type="<<x.type<<" start="<<x.start<<" end="<<x.end<<endl;
      // cout<<" Root  "<<nodeIdx<<" "<<endl;
      // int intleveId=1;
      // int intleveId=log2(nodeIdx);
      int intleveId = nodeIdx;
      omp_set_lock(&lock[intleveId]);
      // if node is null, construct it before inserting
      if(treeNode->at(nodeIdx)==nullptr)
         treeNode->at(nodeIdx) = new Node();
         
      treeNode->at(nodeIdx)->addToCoverList(x);
      omp_unset_lock(&lock[intleveId]);
      return;
   }
   else
   {
      // if (x.type==2) cout<<"inserting id="<<x.id<<" type="<<x.type<<" start="<<x.start<<" end="<<x.end<<endl;
      int leftIdx = 2 * nodeIdx;
      if (leftIdx < treeSize)
      {
         // Node *leftChild=treeNode->at(leftIdx);
         // if (x.intersects(leftChild->interval)){
         if (x.intersects(nodeIntervalArray[leftIdx * 2], nodeIntervalArray[leftIdx * 2 + 1]))
         {
            // #pragma omp task private(x, leftIdx) firstprivate(leftChild) shared(lock)
            // #pragma omp task firstprivate(root) private(x, leftIdx) shared(lock)
            insertEdgeParallel(x, leftIdx);
         }
      }
      int rightIdx = 2 * nodeIdx + 1;
      if (rightIdx < treeSize)
      {
         // Node *rightChild = treeNode->at(rightIdx);
         // if (x.intersects(rightChild->interval)){
         if (x.intersects(nodeIntervalArray[rightIdx * 2], nodeIntervalArray[rightIdx * 2 + 1]))
         {
            // #pragma omp task private(x, rightIdx) firstprivate(rightChild) shared(lock)
            // #pragma omp task firstprivate(root) private(x, rightIdx) shared(lock)
            insertEdgeParallel(x, rightIdx);
         }
      }
   }
}

void SegTree ::insertParallel(Edge x)
{
   insertEdgeParallel(x, 1);
}

void SegTree ::insertAllParallel(Edge *edges)
{
   // int counter = 0;
   // node level lock *************************************
   // omp_lock_t *lock;
   // int lts=log2(treeSize);
   // lock=(omp_lock_t)malloc((lts)*sizeof(omp_lock_t));

   // for(int i=0; i<lts; ++i){
   //    omp_init_lock(&lock[i]);
   // }

   // int ts=treeSize;

   // #pragma acc update self(nodeIntervalArray[0:ts*2])

   // omp_set_dynamic(0);
   // omp_set_num_threads(8);
   // #pragma omp parallel for num_threads(8)
   #pragma omp parallel for
   for (int eid = 0; eid < intervalCount; ++eid)
   {
      insertParallel(edges[eid]);
      // cout<<"inserted "<<counter<<" "<<endl;
      // counter++;
   }

   // // int lts=1;
   // // int lts=log2(treeSize);
   // int lts=treeSize;
   // for(int i=0; i<lts; ++i){
   //    omp_destroy_lock(&lock[i]);
   // }
   // cout <<"test"<< endl;
   // insert(Edge(1,4));
   //  insert(Edge(1,2));
}

// ----------------------------------------------------------
// GPU friendly cover list construction
// 1. find and save node ids of for each edge
// 2. merge and sort ids 
/*// ----------------------------------------------------------
// void SegTree ::saveNodeIds(Edge x) {
//    if (x.contains(nodeIntervalArray[nodeIdx * 2], nodeIntervalArray[nodeIdx * 2 + 1]))    {       
//       treeNode->at(nodeIdx)->addToCoverList(x);
//       return;
//    } else {
//       int leftIdx = 2 * nodeIdx;
//       if (leftIdx < treeSize) {
//          if (x.intersects(nodeIntervalArray[leftIdx * 2], nodeIntervalArray[leftIdx * 2 + 1])) {
//             saveNodeIds(x, leftIdx);
//          }
//       }
//       int rightIdx = 2 * nodeIdx + 1;
//       if (rightIdx < treeSize) {
//          if (x.intersects(nodeIntervalArray[rightIdx * 2], nodeIntervalArray[rightIdx * 2 + 1]))
//          {
//             saveNodeIds(x, rightIdx);
//          }
//       }
//    }

//    int ts=treeSize;
//    for(int nodeIdx=1; nodeIdx<ts*2; ){
//       if (x.contains(nodeIntervalArray[nodeIdx * 2], nodeIntervalArray[nodeIdx * 2 + 1]))    {       
//          treeNode->at(nodeIdx)->addToCoverList(x);
//          return;
//       } else {
//          int leftIdx = 2 * nodeIdx;
//          if (leftIdx < ts) {
//             if (x.intersects(nodeIntervalArray[leftIdx * 2], nodeIntervalArray[leftIdx * 2 + 1])) {
//                nodeIdx=leftIdx;
//             }
//          }
//          int rightIdx = 2 * nodeIdx + 1;
//          if (rightIdx < ts) {
//             if (x.intersects(nodeIntervalArray[rightIdx * 2], nodeIntervalArray[rightIdx * 2 + 1])) {
//                nodeIdx=rightIdx;
//             }
//          }
//       }
//    }
// }

// void SegTree::saveNodeIdsForCoverList(Edge *edges){
//    int ic=intervalCount;
//    // #pragma acc for 
//    for(int eid=0; eid<ic; ++eid){
//       saveNodeIds(edges[eid]);
//    }
// }
// ----------------------------------------------------------
*/

#if GPU!=0
/*
finds the intervals in the GPU given the polygons
*/
void SegTree::copyInputDataToGPU(REAL *bPoly, REAL *cPoly, int bSize, int cSize, int cPolysSize, int cPolyCount, int ic, Edge *edges){
   int bpolys=(bSize+1)*2, cpolys=(cPolysSize)*2;
   intervals2 = (Edge *)malloc(ic*sizeof(Edge));
   high_resolution_clock::time_point cp1, cp2, cp3;

   // cout<<"bpolys="<<bpolys<<" cpolys="<<cpolys<<" ic="<<ic<<endl;


   /*#pragma acc enter data copyin(bPoly[0:bpolys], cPoly[0:cpolys], this)
   #pragma acc enter data create(intervals2[0:ic])
   #pragma acc parallel loop
   for(int i=0; i<bSize; ++i){
      int id=i*2;
      Edge curr;

      curr.getData(id, 0, bPoly[id], bPoly[id+2]);
      intervals2[i]=curr;
   }

   #pragma acc data present(intervals2[0:ic])
   #pragma acc parallel loop
   for(int i=0; i<cSize; ++i){
      int id=i*2;
      Edge curr;

      curr.getData(id, 0, cPoly[id], cPoly[id+2]);
      intervals2[i]=curr;
   }

   #pragma acc update self(intervals2[0:ic])
   for (int i = 0; i < ic; i++)
   {
      if(edges[i].id!=intervals2[i].id || edges[i].start!=intervals2[i].start || edges[i].end!=intervals2[i].end || edges[i].type!=intervals2[i].type){
         cout<<"mis match "<<i<<" "<<edges[i].id<<" "<<intervals2[i].id<<endl;
      }
   }*/

   #ifdef TIME_L2
   cp1 = high_resolution_clock::now();
   #endif   
   // copy intervals to class variable
   for(int i=0; i<ic; ++i){
      intervals2[i]=edges[i];
   }
   #ifdef TIME_L2
   cp2 = high_resolution_clock::now();
   #endif   

   // #pragma acc enter data copyin(this)
   #pragma acc enter data copyin(this, bPoly[0:bpolys], cPoly[0:cpolys], intervals2[0:ic])

   #ifdef TIME_L2
   cp3 = high_resolution_clock::now();
   #endif   
   
   #ifdef DT_INFO
      // Calculating total time taken by the program.
      auto d1 = duration_cast<microseconds>(cp2 - cp1);
      auto dmemcpy = duration_cast<microseconds>(cp3 - cp2);
      printf(" |-intervals copy to array: %d \n |-copy polys and intervals to GPU:  %d\n", static_cast<int>(d1.count()), static_cast<int>(dmemcpy.count()));
   #elif TIME_L2
      auto d1 = duration_cast<microseconds>(cp2 - cp1);
      auto dmemcpy = duration_cast<microseconds>(cp3 - cp2);
      printf("%d ", static_cast<int>(dmemcpy.count()));

   #endif  


}
#endif

#if GPU!=0
// ----------------------------------------------------------
// constructing cover list and end list in parallel at each node
// ----------------------------------------------------------
void SegTree::constructCoverEndLists(REAL *bPoly, REAL *cPoly, int bSize, int cPolysSize){
   int ts=treeSize, ic=intervalCount;
   // int sels = sumEndListSize;
   int *cvTestStart, *cvTestEnd, *edTestStart, *edTestEnd;

   cvTestStart=(int*)malloc(ts*sizeof(int));
   cvTestEnd=(int*)malloc(ts*sizeof(int));
   edTestStart=(int*)malloc(ts*sizeof(int));
   edTestEnd=(int*)malloc(ts*sizeof(int));

   REAL *intervalsStart,  *intervalsEnd, *intervalsStart2,  *intervalsEnd2;
   int *intervalsStartSortedIndexArray, *intervalsEndSortedIndexArray;
   intervalsStart=(REAL*)malloc(ic*sizeof(REAL));
   intervalsEnd=(REAL*)malloc(ic*sizeof(REAL));
   intervalsStart2=(REAL*)malloc(ic*sizeof(REAL));
   intervalsEnd2=(REAL*)malloc(ic*sizeof(REAL));
   intervalsStartSortedIndexArray=(int*)malloc(ic*sizeof(int));
   intervalsEndSortedIndexArray=(int*)malloc(ic*sizeof(int)); 


   #pragma acc data present(intervals2[0:ic])
   #pragma acc enter data create(intervalsStart[0:ic], intervalsEnd[0:ic], intervalsStartSortedIndexArray[0:ic], intervalsEndSortedIndexArray[0:ic], intervalsStart2[0:ic], intervalsEnd2[0:ic])
   #pragma acc parallel loop
   // #pragma omp parallel for
   for(int intid=0; intid<ic; ++intid){
      REAL st=intervals2[intid].start;
      REAL ed=intervals2[intid].end;

      intervalsStart[intid]=st;
      intervalsEnd[intid]=ed;
      intervalsStart2[intid]=st;
      intervalsEnd2[intid]=ed;
      intervalsStartSortedIndexArray[intid]=intid;
      intervalsEndSortedIndexArray[intid]=intid;
   }

   // cout<<"ic="<<ic<<endl;
   // thrustHostStableSortByKey(intervalsStart2, intervalsStartSortedIndexArray, ic);
   // thrustHostStableSortByKey(intervalsEnd2, intervalsEndSortedIndexArray, ic);

   int *endListSizes2;
   endListSizes2=(int*)malloc((ts+1)*sizeof(int));


   #pragma acc host_data use_device(intervalsStart2, intervalsStartSortedIndexArray)
   thrustDeviceStableSortByKey(intervalsStart2, intervalsStartSortedIndexArray, ic);  
   #pragma acc host_data use_device(intervalsEnd2, intervalsEndSortedIndexArray)
   thrustDeviceStableSortByKey(intervalsEnd2, intervalsEndSortedIndexArray, ic);  

   #pragma acc enter data copyin(this)
   #pragma acc enter data create(coverListSizes[0:ts+1], cvTestStart[0:ts], cvTestEnd[0:ts], endListSizes2[0:ts+1], edTestStart[0:ts], edTestEnd[0:ts])

   #pragma acc data present(nodeIntervalArray[0:ts*2], intervalsStart[0:ic], intervalsEnd[0:ic], intervalsStartSortedIndexArray[0:ic], intervalsEndSortedIndexArray[0:ic], intervals2[0:ic])
   #pragma acc parallel loop
   // #pragma omp parallel for
   for(int nodeIdx=0; nodeIdx<ts; ++nodeIdx){
      int clCount=0, clst=0, cled=-1, parent=nodeIdx;
      int edCount=0, edst=0, eded=-1;
      int eidStart=0, eidEnd=ic-1,  edidStart=0, edidEnd=ic-1;
      if(nodeIdx%2!=0) parent=nodeIdx-1; //right child

      int tmpCount=0;
      
      if(nodeIdx==0){
         coverListSizes[nodeIdx]=clCount;
         endListSizes2[nodeIdx]=edCount;
         continue;
      }else if(nodeIdx==1){
         eidStart=0;
         eidEnd=ic-1;
      }
      else if(nodeIdx%2!=0){   //right child

         eidStart=firstOccurrenceBinarySearch(intervalsStart, intervalsStartSortedIndexArray, ic, nodeIntervalArray[parent]);
         eidEnd=lastOccurrenceBinarySearch(intervalsStart, intervalsStartSortedIndexArray, ic, nodeIntervalArray[nodeIdx*2]);
      }else{  //left child
         eidStart=firstOccurrenceBinarySearch(intervalsEnd, intervalsEndSortedIndexArray, ic, nodeIntervalArray[nodeIdx*2+1]);
         eidEnd=lastOccurrenceBinarySearch(intervalsEnd, intervalsEndSortedIndexArray, ic, nodeIntervalArray[parent+1]);
      }

      for(int eeid=eidStart; eeid<=eidEnd; ++eeid){
         int eid=intervalsEndSortedIndexArray[eeid];
         if(nodeIdx%2!=0) eid=intervalsStartSortedIndexArray[eeid];

         if(doesContain(intervals2[eid].start, intervals2[eid].end, nodeIntervalArray[nodeIdx*2], nodeIntervalArray[nodeIdx*2+1]) && (nodeIdx==1 || !doesContain(intervals2[eid].start, intervals2[eid].end, nodeIntervalArray[parent], nodeIntervalArray[parent+1]))){
            clCount+=1;
            cled=eeid;
            if(clCount==1) clst=eeid;
         }

      }
      coverListSizes[nodeIdx]=clCount;
      cvTestStart[nodeIdx]=clst;
      cvTestEnd[nodeIdx]=cled+1;

      for(int eeid=edidStart; eeid<=edidEnd; ++eeid){ 
         int eid=intervalsEndSortedIndexArray[eeid];
         if(nodeIdx%2!=0) eid=intervalsStartSortedIndexArray[eeid];

         if(doesContain(intervals2[eid].start, intervals2[eid].end, nodeIntervalArray[nodeIdx*2], nodeIntervalArray[nodeIdx*2+1]) && (nodeIdx==1 || !doesContain(intervals2[eid].start, intervals2[eid].end, nodeIntervalArray[parent], nodeIntervalArray[parent+1]))){
            tmpCount+=1;
         }

         if(!(doesContain(intervals2[eid].start, intervals2[eid].end, nodeIntervalArray[nodeIdx*2], nodeIntervalArray[nodeIdx*2+1]) && (nodeIdx==1 || !doesContain(intervals2[eid].start, intervals2[eid].end, nodeIntervalArray[parent], nodeIntervalArray[parent+1]))) && doesContainEndPoint(intervals2[eid].start, intervals2[eid].end, nodeIntervalArray[nodeIdx*2], nodeIntervalArray[nodeIdx*2+1])){
            edCount+=1;
            eded=eeid;
            if(edCount==1) edst=eeid;
         }

      }
      endListSizes2[nodeIdx]=edCount;
      edTestStart[nodeIdx]=edst;
      edTestEnd[nodeIdx]=eded+1;
      if(tmpCount!=clCount){
         printf("*********** ERROR------ %d %d %d", tmpCount, clCount, nodeIdx);
      }
   }

   // thrust::exclusive_scan(coverListSizes, coverListSizes + (ts+1), coverListSizes);
   #pragma acc host_data use_device(coverListSizes)
   thrustDeviceExclusiveScan(coverListSizes, (ts+1));  
   #pragma acc update self(coverListSizes[ts])

   #pragma acc host_data use_device(endListSizes2)
   thrustDeviceExclusiveScan(endListSizes2, (ts+1));  
   // #pragma acc update self(endListSizes2[ts])

      // -------------------------------------
   // test endlist sizes
   #pragma acc update self(endListSizes2[0:ts+1], endListSizes[0:ts+1])
   for(int i=1; i<ts; ++i){
      if(endListSizes[i]!=endListSizes2[i]){
         cout<<"i="<<i<<" "<<endListSizes[i]<<" "<<endListSizes2[i]<<endl;
      }
   }
   // -------------------------------------
   
   int scl = coverListSizes[ts]*4;
   // cout<<coverListSizes[ts]<<" scl="<<scl<<endl;
   coverArray = (REAL*)malloc(scl*sizeof(REAL));

   // int sels = endListSizes[ts];
   // cout<<endListSizes[ts]<<" sels="<<sels<<endl;
   // treeNodeEndListArray = (Edge *)malloc((sels) * sizeof(Edge));

   #pragma acc enter data copyin(this)
   #pragma acc enter data create(coverArray[0:scl]/*, treeNodeEndListArray[0:sels]*/)
   #pragma acc data present(nodeIntervalArray[0:ts*2], coverListSizes[0:ts+1], /*edges[0:ic],*/ cvTestStart[0:ts], cvTestEnd[0:ts], /*endListSizes[0:ts+1],*/ edTestStart[0:ts], edTestEnd[0:ts], intervalsStartSortedIndexArray[0:ic], intervalsEndSortedIndexArray[0:ic], intervals2[0:ic])
   #pragma acc parallel loop
   // #pragma omp parallel for  schedule(dynamic, 100)
   for(int nodeIdx=0; nodeIdx<ts; ++nodeIdx){
      if (coverListSizes[nodeIdx+1]!=coverListSizes[nodeIdx]){      
         int clid2=coverListSizes[nodeIdx]*4, parent=nodeIdx;

         if(nodeIdx==0) continue;

         if(nodeIdx%2!=0) parent=nodeIdx-1;

         // cout<<"clid2="<<clid2<<" coverArray[nodeIdx]="<<coverArray[nodeIdx]<<endl;
         // for (int eid=0; eid<ic; ++eid){
         for (int eeid=cvTestStart[nodeIdx]; eeid<cvTestEnd[nodeIdx]; ++eeid){
            int eid=intervalsEndSortedIndexArray[eeid];
            if(nodeIdx%2!=0) eid=intervalsStartSortedIndexArray[eeid];

            if(doesContain(intervals2[eid].start, intervals2[eid].end, nodeIntervalArray[nodeIdx*2], nodeIntervalArray[nodeIdx*2+1]) && (nodeIdx==1 || !doesContain(intervals2[eid].start, intervals2[eid].end, nodeIntervalArray[parent], nodeIntervalArray[parent+1]))){            
               coverArray[clid2] = intervals2[eid].id;
               coverArray[clid2 + 1] = intervals2[eid].type;
               coverArray[clid2 + 2] = intervals2[eid].start;
               coverArray[clid2 + 3] = intervals2[eid].end;
               clid2+=4;
            }
         }
      }

      // if (endListSizes[nodeIdx+1]!=endListSizes[nodeIdx]){      
      //    int edid2=endListSizes[nodeIdx], parent=nodeIdx;

      //    if(nodeIdx==0) continue;

      //    if(nodeIdx%2!=0) parent=nodeIdx-1;

      //    // cout<<"clid2="<<clid2<<" coverArray[nodeIdx]="<<coverArray[nodeIdx]<<endl;
      //    // for (int eid=0; eid<ic; ++eid){
      //    for (int eeid=edTestStart[nodeIdx]; eeid<edTestEnd[nodeIdx]; ++eeid){
      //       int eid=intervalsEndSortedIndexArray[eeid];
      //       if(nodeIdx%2!=0) eid=intervalsStartSortedIndexArray[eeid];

      //       if(!(doesContain(intervals2[eid].start, intervals2[eid].end, nodeIntervalArray[nodeIdx*2], nodeIntervalArray[nodeIdx*2+1]) && (nodeIdx==1 || !doesContain(intervals2[eid].start, intervals2[eid].end, nodeIntervalArray[parent], nodeIntervalArray[parent+1]))) && doesContainEndPoint(intervals2[eid].start, intervals2[eid].end, nodeIntervalArray[nodeIdx*2], nodeIntervalArray[nodeIdx*2+1])){            
      //          treeNodeEndListArray[edid2] = intervals2[eid];
      //          edid2+=1;
      //       }
      //    }
      // }
   }

   sumCoverListSize = coverListSizes[ts];
   // sumEndListSize = endListSizes[ts];


   int bpolys = (bSize + 1) * 2, cpolys = (cPolysSize) * 2;
   // #pragma acc update device(treeNodeEndListArray[0 : sels], endListSizes[0 : ts], endListCounts[0 : ts])
   // #pragma acc enter data copyin(coverArray[0 : scl], bPoly[0 : bpolys], cPoly[0 : cpolys], coverListSizes[0:ts+1])

   // #pragma acc update self(coverArray[0:scl])
   #pragma acc enter data copyin(bPoly[0 : bpolys], cPoly[0 : cpolys])

}
#endif

#if GPU!=0
// ----------------------------------------------------------
// constructing cover list in parallel at each node
// ----------------------------------------------------------
void SegTree::constructCoverListGPU(){
   int ts=treeSize, ic=intervalCount;
   int sels = sumEndListSize;
   int *cvTestStart, *cvTestEnd;

   cvTestStart=(int*)malloc(ts*sizeof(int));
   cvTestEnd=(int*)malloc(ts*sizeof(int));

   REAL *intervalsStart,  *intervalsEnd, *intervalsStart2,  *intervalsEnd2;
   int *intervalsStartSortedIndexArray, *intervalsEndSortedIndexArray;
   intervalsStart=(REAL*)malloc(ic*sizeof(REAL));
   intervalsEnd=(REAL*)malloc(ic*sizeof(REAL));
   intervalsStart2=(REAL*)malloc(ic*sizeof(REAL));
   intervalsEnd2=(REAL*)malloc(ic*sizeof(REAL));
   intervalsStartSortedIndexArray=(int*)malloc(ic*sizeof(int));
   intervalsEndSortedIndexArray=(int*)malloc(ic*sizeof(int)); 

   #ifdef TIME_L2 
   high_resolution_clock::time_point cp1, cp2;
   #endif


   #pragma acc data present(intervals2[0:ic])
   #pragma acc enter data create(intervalsStart[0:ic], intervalsEnd[0:ic], intervalsStartSortedIndexArray[0:ic], intervalsEndSortedIndexArray[0:ic], intervalsStart2[0:ic], intervalsEnd2[0:ic])
   #pragma acc parallel loop
   // #pragma omp parallel for
   for(int intid=0; intid<ic; ++intid){
      REAL st=intervals2[intid].start;
      REAL ed=intervals2[intid].end;

      intervalsStart[intid]=st;
      intervalsEnd[intid]=ed;
      intervalsStart2[intid]=st;
      intervalsEnd2[intid]=ed;
      intervalsStartSortedIndexArray[intid]=intid;
      intervalsEndSortedIndexArray[intid]=intid;
   }

   // cout<<"ic="<<ic<<endl;
   // thrustHostStableSortByKey(intervalsStart2, intervalsStartSortedIndexArray, ic);
   // thrustHostStableSortByKey(intervalsEnd2, intervalsEndSortedIndexArray, ic);

   #pragma acc host_data use_device(intervalsStart2, intervalsStartSortedIndexArray)
   thrustDeviceStableSortByKey(intervalsStart2, intervalsStartSortedIndexArray, ic);  
   #pragma acc host_data use_device(intervalsEnd2, intervalsEndSortedIndexArray)
   thrustDeviceStableSortByKey(intervalsEnd2, intervalsEndSortedIndexArray, ic);  

   #pragma acc enter data copyin(this)
   #pragma acc enter data create(coverListSizes[0:ts+1], cvTestStart[0:ts], cvTestEnd[0:ts])
   #pragma acc data present(nodeIntervalArray[0:ts*2], intervalsStart[0:ic], intervalsEnd[0:ic], intervalsStartSortedIndexArray[0:ic], intervalsEndSortedIndexArray[0:ic], intervals2[0:ic])
   #pragma acc parallel loop
   // #pragma omp parallel for
   for(int nodeIdx=0; nodeIdx<ts; ++nodeIdx){
      int clCount=0, clst=0, cled=-1, parent=nodeIdx;
      int eidStart=0, eidEnd=ic-1;
      if(nodeIdx%2!=0) parent=nodeIdx-1; //right child
      
      if(nodeIdx==0){
          coverListSizes[nodeIdx]=clCount;
         continue;
      }else if(nodeIdx==1){
         eidStart=0;
         eidEnd=ic-1;
      }
      else if(nodeIdx%2!=0){   //right child

         eidStart=firstOccurrenceBinarySearch(intervalsStart, intervalsStartSortedIndexArray, ic, nodeIntervalArray[parent]);
         eidEnd=lastOccurrenceBinarySearch(intervalsStart, intervalsStartSortedIndexArray, ic, nodeIntervalArray[nodeIdx*2]);
      }else{  //left child
         eidStart=firstOccurrenceBinarySearch(intervalsEnd, intervalsEndSortedIndexArray, ic, nodeIntervalArray[nodeIdx*2+1]);
         eidEnd=lastOccurrenceBinarySearch(intervalsEnd, intervalsEndSortedIndexArray, ic, nodeIntervalArray[parent+1]);
      }

      for(int eeid=eidStart; eeid<=eidEnd; ++eeid){
         int eid=intervalsEndSortedIndexArray[eeid];
         if(nodeIdx%2!=0) eid=intervalsStartSortedIndexArray[eeid];

         if(doesContain(intervals2[eid].start, intervals2[eid].end, nodeIntervalArray[nodeIdx*2], nodeIntervalArray[nodeIdx*2+1]) && (nodeIdx==1 || !doesContain(intervals2[eid].start, intervals2[eid].end, nodeIntervalArray[parent], nodeIntervalArray[parent+1]))){
            clCount+=1;
            cled=eeid;
            if(clCount==1) clst=eeid;
         }
      }
      coverListSizes[nodeIdx]=clCount;
      cvTestStart[nodeIdx]=clst;
      cvTestEnd[nodeIdx]=cled+1;
   }

   // test with old DS
   // #pragma acc update self(coverListSizes[0:ts+1])

   /*for (int nid = 1; nid < treeSize; ++nid)
   {
      // cout<<"nid="<<nid<<" "<<endl;
      // cout<<"nid="<<nid<<" "<<(treeNode->at(nid-1)==nullptr)<<endl;
      if (treeNode->at(nid)!=nullptr){
         // cout<<"nid="<<nid<<" "<<endl;
         // cout<<"nid="<<nid<<" "<<treeNode->at(nid-1)->coverList->size()<<endl;
         if(coverListSizes[nid]!=treeNode->at(nid)->coverList->size()){
            cout<<"error nid="<<nid<<" "<<treeNode->at(nid)->coverList->size()<<" "<<coverListSizes[nid]<<endl;
         }
      }else{
         if(coverListSizes[nid]!=0){
            cout<<"error nid="<<nid<<" "<<coverListSizes[nid]<<endl;
         }
      }
   }
   cout<<"testing complete "<<endl;*/

   // thrust::exclusive_scan(coverListSizes, coverListSizes + (ts+1), coverListSizes);
   #pragma acc host_data use_device(coverListSizes)
   thrustDeviceExclusiveScan(coverListSizes, (ts+1));  

   #ifdef TIME_L2 
   cp1 = high_resolution_clock::now();
   #endif
   #pragma acc update self(coverListSizes[ts])
   #ifdef TIME_L2 
   cp2 = high_resolution_clock::now();
   #endif

   int scl = coverListSizes[ts]*4;
   // cout<<coverListSizes[ts]<<" scl="<<scl<<endl;
   coverArray = (REAL*)malloc(scl*sizeof(REAL));

   #pragma acc enter data copyin(this)
   #pragma acc enter data create(coverArray[0:scl])
   #pragma acc data present(nodeIntervalArray[0:ts*2], coverListSizes[0:ts+1], /*edges[0:ic],*/ cvTestStart[0:ts], cvTestEnd[0:ts], intervalsStartSortedIndexArray[0:ic], intervalsEndSortedIndexArray[0:ic], intervals2[0:ic])
   #pragma acc parallel loop
   // #pragma omp parallel for  schedule(dynamic, 100)
   for(int nodeIdx=0; nodeIdx<ts; ++nodeIdx){
      if (coverListSizes[nodeIdx+1]!=coverListSizes[nodeIdx]){      
         int clid2=coverListSizes[nodeIdx]*4, parent=nodeIdx;

         if(nodeIdx==0) continue;

         if(nodeIdx%2!=0) parent=nodeIdx-1;

         // cout<<"clid2="<<clid2<<" coverArray[nodeIdx]="<<coverArray[nodeIdx]<<endl;
         // for (int eid=0; eid<ic; ++eid){
         for (int eeid=cvTestStart[nodeIdx]; eeid<cvTestEnd[nodeIdx]; ++eeid){
            int eid=intervalsEndSortedIndexArray[eeid];
            if(nodeIdx%2!=0) eid=intervalsStartSortedIndexArray[eeid];

            if(doesContain(intervals2[eid].start, intervals2[eid].end, nodeIntervalArray[nodeIdx*2], nodeIntervalArray[nodeIdx*2+1]) && (nodeIdx==1 || !doesContain(intervals2[eid].start, intervals2[eid].end, nodeIntervalArray[parent], nodeIntervalArray[parent+1]))){            
               coverArray[clid2] = intervals2[eid].id;
               coverArray[clid2 + 1] = intervals2[eid].type;
               coverArray[clid2 + 2] = intervals2[eid].start;
               coverArray[clid2 + 3] = intervals2[eid].end;
               clid2+=4;
            }
         }
      }
   }

   sumCoverListSize = coverListSizes[ts];

   // #pragma acc update device(treeNodeEndListArray[0 : sels], endListSizes[0 : ts], endListCounts[0 : ts])
   // #pragma acc enter data copyin(coverArray[0 : scl], bPoly[0 : bpolys], cPoly[0 : cpolys], coverListSizes[0:ts+1])

   // #pragma acc update self(coverArray[0:scl])

   #ifdef DT_INFO 
    auto d1 = duration_cast<microseconds>(cp2 - cp1);
    cout<<" |-Time: cover list sizes max memcpy to CPU: "<<d1.count() << setprecision(10) << endl;
   #elif TIME_L2
      auto d1 = duration_cast<microseconds>(cp2 - cp1);
      cout<<d1.count() << setprecision(10)<<" ";
   #endif

}
#endif

// ----------------------------------------------------------
// constructing cover list in parallel at each node: Multi-core version
// ----------------------------------------------------------
void SegTree::constructCoverListMulticore(Edge *edges){
   int ts=treeSize, ic=intervalCount;
   int sels = sumEndListSize;
   int *cvTestStart, *cvTestEnd;

   cvTestStart=(int*)malloc(ts*sizeof(int));
   cvTestEnd=(int*)malloc(ts*sizeof(int));

   REAL *intervalsStart,  *intervalsEnd, *intervalsStart2,  *intervalsEnd2;
   int *intervalsStartSortedIndexArray, *intervalsEndSortedIndexArray;
   intervalsStart=(REAL*)malloc(ic*sizeof(REAL));
   intervalsEnd=(REAL*)malloc(ic*sizeof(REAL));
   intervalsStart2=(REAL*)malloc(ic*sizeof(REAL));
   intervalsEnd2=(REAL*)malloc(ic*sizeof(REAL));
   intervalsStartSortedIndexArray=(int*)malloc(ic*sizeof(int));
   intervalsEndSortedIndexArray=(int*)malloc(ic*sizeof(int)); 

   #pragma acc update self(nodeIntervalArray[0:ts*2])

   #pragma omp parallel for
   for(int intid=0; intid<ic; ++intid){
      REAL st=edges[intid].start;
      REAL ed=edges[intid].end;

      intervalsStart[intid]=st;
      intervalsEnd[intid]=ed;
      intervalsStart2[intid]=st;
      intervalsEnd2[intid]=ed;
      intervalsStartSortedIndexArray[intid]=intid;
      intervalsEndSortedIndexArray[intid]=intid;
   }

   // cout<<"ic="<<ic<<endl;
   thrust::stable_sort_by_key(intervalsStart2, intervalsStart2+ic, intervalsStartSortedIndexArray);
   thrust::stable_sort_by_key(intervalsEnd2, intervalsEnd2+ic, intervalsEndSortedIndexArray);

   #pragma omp parallel for schedule(dynamic, 100)
   for(int nodeIdx=0; nodeIdx<ts; ++nodeIdx){
      int clCount=0, clst=0, cled=-1, parent=nodeIdx;
      int eidStart=0, eidEnd=ic-1;
      if(nodeIdx%2!=0) parent=nodeIdx-1; //right child
      
      if(nodeIdx==0){
          coverListSizes[nodeIdx]=clCount;
         continue;
      }else if(nodeIdx==1){
         eidStart=0;
         eidEnd=ic-1;
      }
      else if(nodeIdx%2!=0){   //right child

         eidStart=firstOccurrenceBinarySearch(intervalsStart, intervalsStartSortedIndexArray, ic, nodeIntervalArray[parent]);
         eidEnd=lastOccurrenceBinarySearch(intervalsStart, intervalsStartSortedIndexArray, ic, nodeIntervalArray[nodeIdx*2]);
      }else{  //left child
         eidStart=firstOccurrenceBinarySearch(intervalsEnd, intervalsEndSortedIndexArray, ic, nodeIntervalArray[nodeIdx*2+1]);
         eidEnd=lastOccurrenceBinarySearch(intervalsEnd, intervalsEndSortedIndexArray, ic, nodeIntervalArray[parent+1]);
      }

      for(int eeid=eidStart; eeid<=eidEnd; ++eeid){
         int eid=intervalsEndSortedIndexArray[eeid];
         if(nodeIdx%2!=0) eid=intervalsStartSortedIndexArray[eeid];

         if(doesContain(edges[eid].start, edges[eid].end, nodeIntervalArray[nodeIdx*2], nodeIntervalArray[nodeIdx*2+1]) && (nodeIdx==1 || !doesContain(edges[eid].start, edges[eid].end, nodeIntervalArray[parent], nodeIntervalArray[parent+1]))){
            clCount+=1;
            cled=eeid;
            if(clCount==1) clst=eeid;
         }
      }
      coverListSizes[nodeIdx]=clCount;
      cvTestStart[nodeIdx]=clst;
      cvTestEnd[nodeIdx]=cled+1;
   }

   thrust::exclusive_scan(coverListSizes, coverListSizes + (ts+1), coverListSizes);
   
   int scl = coverListSizes[ts]*4;
   // cout<<coverListSizes[ts]<<" scl="<<scl<<endl;
   coverArray = (REAL*)malloc(scl*sizeof(REAL));

   #pragma omp parallel for  schedule(dynamic, 100)
   for(int nodeIdx=0; nodeIdx<ts; ++nodeIdx){
      if (coverListSizes[nodeIdx+1]!=coverListSizes[nodeIdx]){      
         int clid2=coverListSizes[nodeIdx]*4, parent=nodeIdx;

         if(nodeIdx==0) continue;

         if(nodeIdx%2!=0) parent=nodeIdx-1;

         // cout<<"clid2="<<clid2<<" coverArray[nodeIdx]="<<coverArray[nodeIdx]<<endl;
         // for (int eid=0; eid<ic; ++eid){
         for (int eeid=cvTestStart[nodeIdx]; eeid<cvTestEnd[nodeIdx]; ++eeid){
            int eid=intervalsEndSortedIndexArray[eeid];
            if(nodeIdx%2!=0) eid=intervalsStartSortedIndexArray[eeid];

            if(doesContain(edges[eid].start, edges[eid].end, nodeIntervalArray[nodeIdx*2], nodeIntervalArray[nodeIdx*2+1]) && (nodeIdx==1 || !doesContain(edges[eid].start, edges[eid].end, nodeIntervalArray[parent], nodeIntervalArray[parent+1]))){            
               coverArray[clid2] = edges[eid].id;
               coverArray[clid2 + 1] = edges[eid].type;
               coverArray[clid2 + 2] = edges[eid].start;
               coverArray[clid2 + 3] = edges[eid].end;
               clid2+=4;
            }
         }
      }
   }

   sumCoverListSize = coverListSizes[ts];

   #pragma acc enter data copyin(coverArray[0 : scl], coverListSizes[0:ts+1])

}

// ----------------------------------------------------------
// constructing cover list in parallel at each node: Serial version
// ----------------------------------------------------------
void SegTree::constructCoverListSerial(Edge *edges){
   int ts=treeSize, ic=intervalCount;
   int sels = sumEndListSize;
   int *cvTestStart, *cvTestEnd;

   cvTestStart=(int*)malloc(ts*sizeof(int));
   cvTestEnd=(int*)malloc(ts*sizeof(int));

   REAL *intervalsStart,  *intervalsEnd, *intervalsStart2,  *intervalsEnd2;
   int *intervalsStartSortedIndexArray, *intervalsEndSortedIndexArray;
   intervalsStart=(REAL*)malloc(ic*sizeof(REAL));
   intervalsEnd=(REAL*)malloc(ic*sizeof(REAL));
   intervalsStart2=(REAL*)malloc(ic*sizeof(REAL));
   intervalsEnd2=(REAL*)malloc(ic*sizeof(REAL));
   intervalsStartSortedIndexArray=(int*)malloc(ic*sizeof(int));
   intervalsEndSortedIndexArray=(int*)malloc(ic*sizeof(int)); 

   #pragma acc update self(nodeIntervalArray[0:ts*2])

   for(int intid=0; intid<ic; ++intid){
      REAL st=edges[intid].start;
      REAL ed=edges[intid].end;

      intervalsStart[intid]=st;
      intervalsEnd[intid]=ed;
      intervalsStart2[intid]=st;
      intervalsEnd2[intid]=ed;
      intervalsStartSortedIndexArray[intid]=intid;
      intervalsEndSortedIndexArray[intid]=intid;
   }

   // cout<<"ic="<<ic<<endl;
   thrust::stable_sort_by_key(intervalsStart2, intervalsStart2+ic, intervalsStartSortedIndexArray);
   thrust::stable_sort_by_key(intervalsEnd2, intervalsEnd2+ic, intervalsEndSortedIndexArray);

   for(int nodeIdx=0; nodeIdx<ts; ++nodeIdx){
      int clCount=0, clst=0, cled=-1, parent=nodeIdx;
      int eidStart=0, eidEnd=ic-1;
      if(nodeIdx%2!=0) parent=nodeIdx-1; //right child
      
      if(nodeIdx==0){
          coverListSizes[nodeIdx]=clCount;
         continue;
      }else if(nodeIdx==1){
         eidStart=0;
         eidEnd=ic-1;
      }
      else if(nodeIdx%2!=0){   //right child

         eidStart=firstOccurrenceBinarySearch(intervalsStart, intervalsStartSortedIndexArray, ic, nodeIntervalArray[parent]);
         eidEnd=lastOccurrenceBinarySearch(intervalsStart, intervalsStartSortedIndexArray, ic, nodeIntervalArray[nodeIdx*2]);
      }else{  //left child
         eidStart=firstOccurrenceBinarySearch(intervalsEnd, intervalsEndSortedIndexArray, ic, nodeIntervalArray[nodeIdx*2+1]);
         eidEnd=lastOccurrenceBinarySearch(intervalsEnd, intervalsEndSortedIndexArray, ic, nodeIntervalArray[parent+1]);
      }

      for(int eeid=eidStart; eeid<=eidEnd; ++eeid){
         int eid=intervalsEndSortedIndexArray[eeid];
         if(nodeIdx%2!=0) eid=intervalsStartSortedIndexArray[eeid];

         if(doesContain(edges[eid].start, edges[eid].end, nodeIntervalArray[nodeIdx*2], nodeIntervalArray[nodeIdx*2+1]) && (nodeIdx==1 || !doesContain(edges[eid].start, edges[eid].end, nodeIntervalArray[parent], nodeIntervalArray[parent+1]))){
            clCount+=1;
            cled=eeid;
            if(clCount==1) clst=eeid;
         }
      }
      coverListSizes[nodeIdx]=clCount;
      cvTestStart[nodeIdx]=clst;
      cvTestEnd[nodeIdx]=cled+1;
   }

   thrust::exclusive_scan(coverListSizes, coverListSizes + (ts+1), coverListSizes);
   
   int scl = coverListSizes[ts]*4;
   // cout<<coverListSizes[ts]<<" scl="<<scl<<endl;
   coverArray = (REAL*)malloc(scl*sizeof(REAL));

   for(int nodeIdx=0; nodeIdx<ts; ++nodeIdx){
      if (coverListSizes[nodeIdx+1]!=coverListSizes[nodeIdx]){      
         int clid2=coverListSizes[nodeIdx]*4, parent=nodeIdx;

         if(nodeIdx==0) continue;

         if(nodeIdx%2!=0) parent=nodeIdx-1;

         // cout<<"clid2="<<clid2<<" coverArray[nodeIdx]="<<coverArray[nodeIdx]<<endl;
         // for (int eid=0; eid<ic; ++eid){
         for (int eeid=cvTestStart[nodeIdx]; eeid<cvTestEnd[nodeIdx]; ++eeid){
            int eid=intervalsEndSortedIndexArray[eeid];
            if(nodeIdx%2!=0) eid=intervalsStartSortedIndexArray[eeid];

            if(doesContain(edges[eid].start, edges[eid].end, nodeIntervalArray[nodeIdx*2], nodeIntervalArray[nodeIdx*2+1]) && (nodeIdx==1 || !doesContain(edges[eid].start, edges[eid].end, nodeIntervalArray[parent], nodeIntervalArray[parent+1]))){            
               coverArray[clid2] = edges[eid].id;
               coverArray[clid2 + 1] = edges[eid].type;
               coverArray[clid2 + 2] = edges[eid].start;
               coverArray[clid2 + 3] = edges[eid].end;
               clid2+=4;
            }
         }
      }
   }

   sumCoverListSize = coverListSizes[ts];

   #pragma acc enter data copyin(coverArray[0 : scl], coverListSizes[0:ts+1])

}


/*void SegTree::constructCoverList(Edge *edges, REAL *bPoly, REAL *cPoly, int bSize, int cPolysSize){
   int ts=treeSize, ic=intervalCount;
   int sels = sumEndListSize;
   int *cvTestStart, *cvTestEnd;

   cvTestStart=(int*)malloc(ts*sizeof(int));
   cvTestEnd=(int*)malloc(ts*sizeof(int));


   #pragma acc enter data copyin(this)
   #pragma acc enter data create(coverListSizes[0:ts+1], cvTestStart[0:ts], cvTestEnd[0:ts]) copyin(edges[0:ic])
   #pragma acc data present(nodeIntervalArray[0:ts*2])
   #pragma acc parallel loop
   // #pragma omp parallel for
   for(int nodeIdx=0; nodeIdx<ts; ++nodeIdx){
      int clCount=0, clst=0, cled=-1, parent=nodeIdx;
      if(nodeIdx%2!=0) parent=nodeIdx-1;

      for (int eid=0; nodeIdx!=0 && eid<ic; ++eid){
         if(nodeIntervalArray[nodeIdx*2]==nodeIntervalArray[nodeIdx*2+1] && edges[eid].start==edges[eid].end && nodeIntervalArray[nodeIdx*2]==edges[eid].start){
            clCount+=1;
            cled=eid;
            if(clCount==1) clst=eid;
         }else if(doesContain(edges[eid].start, edges[eid].end, nodeIntervalArray[nodeIdx*2], nodeIntervalArray[nodeIdx*2+1]) && (nodeIdx==1 || !doesContain(edges[eid].start, edges[eid].end, nodeIntervalArray[parent], nodeIntervalArray[parent+1]))){
            clCount+=1;
            cled=eid;
            if(clCount==1) clst=eid;
         }
      }
      coverListSizes[nodeIdx]=clCount;
      cvTestStart[nodeIdx]=clst;
      cvTestEnd[nodeIdx]=cled+1;
   }

   // thrust::exclusive_scan(coverListSizes, coverListSizes + (ts+1), coverListSizes);
   #pragma acc host_data use_device(coverListSizes)
   thrustDeviceExclusiveScan(coverListSizes, (ts+1));  
   #pragma acc update self(coverListSizes[0:ts+1])

   
   int scl = coverListSizes[ts]*4;
   cout<<"scl="<<scl<<endl;
   coverArray = (REAL*)malloc(scl*sizeof(REAL));
   cout<<"scl="<<scl<<endl;

   // #pragma acc update device(treeNodeEndListArray[0 : sels],  endListSizes[0 : ts], endListCounts[0 : ts])

   #pragma acc enter data copyin(this)
   #pragma acc enter data create(coverArray[0:scl])
   #pragma acc data present(nodeIntervalArray[0:ts*2], coverListSizes[0:ts+1], edges[0:ic], cvTestStart[0:ts], cvTestEnd[0:ts])
   #pragma acc parallel loop
   // #pragma omp parallel for  schedule(dynamic, 100)
   for(int nodeIdx=0; nodeIdx<ts; ++nodeIdx){
      if (coverListSizes[nodeIdx+1]!=coverListSizes[nodeIdx]){      
         int clid2=coverListSizes[nodeIdx]*4, parent=nodeIdx;
         if(nodeIdx%2!=0) parent=nodeIdx-1;

         // cout<<"clid2="<<clid2<<" coverArray[nodeIdx]="<<coverArray[nodeIdx]<<endl;
         // for (int eid=0; nodeIdx!=0 && eid<ic; ++eid){
         for (int eid=cvTestStart[nodeIdx]; eid<cvTestEnd[nodeIdx]; ++eid){
         // for (int eid=0; eid<cvTestEnd[nodeIdx]; ++eid){
            if(nodeIntervalArray[nodeIdx*2]==nodeIntervalArray[nodeIdx*2+1] && edges[eid].start==edges[eid].end && nodeIntervalArray[nodeIdx*2]==edges[eid].start){
               coverArray[clid2] = edges[eid].id;
               coverArray[clid2 + 1] = edges[eid].type;
               coverArray[clid2 + 2] = edges[eid].start;
               coverArray[clid2 + 3] = edges[eid].end;
               clid2+=4;
            }else if(doesContain(edges[eid].start, edges[eid].end, nodeIntervalArray[nodeIdx*2], nodeIntervalArray[nodeIdx*2+1]) && (nodeIdx==1 || !doesContain(edges[eid].start, edges[eid].end, nodeIntervalArray[parent], nodeIntervalArray[parent+1]))){            
               coverArray[clid2] = edges[eid].id;
               coverArray[clid2 + 1] = edges[eid].type;
               coverArray[clid2 + 2] = edges[eid].start;
               coverArray[clid2 + 3] = edges[eid].end;
               clid2+=4;
            }
         }
      }
   }

   sumCoverListSize = coverListSizes[ts];
   // cout<<"countEndlistCoverListSizes: counted"<<endl;
   // #pragma acc enter data copyin(coverListSizes[0 : (ts + 1)])

   int bpolys = (bSize + 1) * 2, cpolys = (cPolysSize) * 2;
   // #pragma acc update device(treeNodeEndListArray[0 : sels], endListSizes[0 : ts], endListCounts[0 : ts])
   // #pragma acc enter data copyin(coverArray[0 : scl], bPoly[0 : bpolys], cPoly[0 : cpolys])
   
   // #pragma acc update self(endListSizes[0 : ts], endListCounts[0 : ts])

   #pragma acc update self(coverArray[0:scl])
   #pragma acc update device(treeNodeEndListArray[0 : sels], endListSizes[0 : ts], endListCounts[0 : ts])
   #pragma acc enter data copyin(bPoly[0 : bpolys], cPoly[0 : cpolys])

}*/
// ----------------------------------------------------------


void SegTree ::insertEdge(Edge x, Node *root, int nodeIdx)
{
   if (x.contains(root->interval))
   {
      root->addToCoverList(x);
      return;
   }
   else
   {
      int leftIdx = 2 * nodeIdx;
      if (leftIdx < treeSize)
      {
         Node *leftChild = treeNode->at(leftIdx);
         if (x.intersects(leftChild->interval))
         {
            insertEdge(x, leftChild, leftIdx);
         }
      }
      int rightIdx = 2 * nodeIdx + 1;
      if (rightIdx < treeSize)
      {
         Node *rightChild = treeNode->at(rightIdx);
         if (x.intersects(rightChild->interval))
         {
            insertEdge(x, rightChild, rightIdx);
         }
      }
   }
}

void SegTree ::insert(Edge x)
{
   insertEdge(x, treeNode->at(1), 1);
}

void SegTree ::insertAll(Edge *edges)
{
   for (int eid = 0; eid < intervalCount; ++eid)
   {
      insert(edges[eid]);
   }
}

// -------------------------------------

// original
vector<REAL> *SegTree ::getPointsFromEdges(vector<Edge> *intervals)
{
   set<REAL> *ptList = new set<REAL>();
   for (Edge e : *intervals)
   {
      ptList->insert(e.start);
      ptList->insert(e.end);
   }
   vector<REAL> *ptVect = new vector<REAL>(ptList->size());
   //**************** why copy is needed here? just return the ptList??? ****************
   std::copy(ptList->begin(), ptList->end(), ptVect->begin());

   //  sort(ptList->begin(), ptList->end());
   //  for(Edge e : *intervals)
   //    {
   //       ptList->push_back(e.start);
   //       ptList->push_back(e.end);
   //    }
   //
   //    sort(ptList->begin(), ptList->end());

   return ptVect;
}

// considers duplicate start and end points
// Buddhi
vector<REAL> *SegTree ::getPointsFromEdges2(Edge *intervals)
{
   REAL intLength, minIntLength;
   bool min = false;
   set<REAL> dupEdgeList;
   set<REAL> *ptList = new set<REAL>();
   Edge e;

   // for (Edge e : *intervals)
   for (int eid = 0; eid < intervalCount; ++eid)
   {
      e = intervals[eid];
      ptList->insert(e.start);
      ptList->insert(e.end);
      if (e.start == e.end)
      {
         dupEdgeList.insert(e.start);
      }
      if (!min)
      {
         if (e.start != e.end)
         {
            minIntLength = e.end - e.start;
            min = true;
         }
      }
      else
      {
         intLength = e.end - e.start;
         if (intLength != 0 && intLength < minIntLength)
         {
            minIntLength = intLength;
         }
      }
      // cout<<"test length "<<intLength<<endl;
      // cout<<"test min length "<<minIntLength<<endl;
   }
   delta = minIntLength / 10;
   vector<REAL> *ptVect = new vector<REAL>(ptList->size());

   std::copy(ptList->begin(), ptList->end(), ptVect->begin());
   for (REAL dupEle : dupEdgeList)
   {
      ptVect->push_back(dupEle);
      // cout<<"dup point "<<dupEle<<endl;
   }
   sort(ptVect->begin(), ptVect->end());

   return ptVect;
}

/*
copy elementary intervals to array and GPU
*/
void SegTree::copyToEleArray(vector<Edge> *elementaryIntervals)
{
   // int eleCount=elementaryIntervals->size();
   // int pCount=m_elementaryPoints->size();
   // elementaryPointCount=eleCount;

   // elementaryEdgeArray=(Edge *)malloc((eleCount)*sizeof(Edge));
   // elementaryArray=(REAL *)malloc((pCount)*sizeof(REAL));

   // // elementaryEdgeArray=&elementaryIntervals[0];
   // // #pragma omp parallel for
   // for(int i=0; i<eleCount; ++i){
   //    elementaryEdgeArray[i]=elementaryIntervals->at(i);
   //    elementaryArray[i]=m_elementaryPoints->at(i);
   // }
   // elementaryArray[eleCount]=m_elementaryPoints->at(eleCount);

   // #pragma acc enter data copyin(elementaryEdgeArray[0:eleCount], elementaryArray[0:pCount])
}

/*
paralllel interval construction
*/
/*
void SegTree::intervalConstructionParallel
   (REAL *bPoly, REAL *cPoly, int bSize, int cSize,
   Edge **intervals){

   int ic=bSize+cSize;
   int bpolys=(bSize+1)*2, cpolys=(cSize+1)*2;
   cout<<ic<<" === "<<intervalCount<<endl;
   // intervalCount=ic;
   *intervals=(Edge *)malloc((ic)*sizeof(Edge));

   #pragma acc enter data copyin(bPoly[0:bpolys], cPoly[0:cpolys])
   #pragma acc enter data create(intervals[0][0:ic])
   #pragma acc parallel loop async
   for(int i=0; i<bSize; ++i){
      int i2=i*2;
      Edge curr;
      curr.id=i2;
      curr.type=0;
      if(bPoly[i2]<=bPoly[i2+2]){
         curr.start=bPoly[i2];
         curr.end=bPoly[i2+2];
      }else{
         curr.end=bPoly[i2];
         curr.start=bPoly[i2+2];
      }
      // curr.getData(i2, 0, bPoly[i2], bPoly[i2+2]);
      *(*intervals+i)=curr;
      // cout<<bPoly[i2]<<" "<<bPoly[i2+2]<<endl;
   }
   #pragma acc parallel loop async
   for(int i=0; i<cSize; ++i){
      int i2=i*2;
      Edge curr;
      curr.id=i2;
      curr.type=1;
      if(cPoly[i2]<=cPoly[i2+2]){
         curr.start=cPoly[i2];
         curr.end=cPoly[i2+2];
      }else{
         curr.end=cPoly[i2];
         curr.start=cPoly[i2+2];
      }
      // curr.getData(i2, 1, cPoly[i2], cPoly[i2+2]);
      *(*intervals+i+bSize)=curr;
      // cout<<cPoly[i2]<<" "<<cPoly[i2+2]<<endl;
   }
   #pragma acc update self(intervals[0][0:ic])
   cout<<"Done"<<endl;
}
*/
/*
Uses CMBR filtering
*/
void SegTree::constructElementaryIntervalsArrayCMBR(Edge *intervals, REAL *cmbr)
{
   int eleCount = 0, edgeCount = 0, prevCount;
   REAL currEle, *eleArray;
   Edge *eleEdgeArray;
   LONG *eleIntArray;
   int *sortedIndex, lastVertical = -1, ic = intervalCount;
   int icupper = Util::getNPowOf2(ic);
   intervalCountUpperBound = icupper;
   vector<Edge> filteredIntervals;

   for (int i = 0; i < ic; ++i)
   {
      if (cmbr[0] <= intervals[i].end && cmbr[2] >= intervals[i].start)
      {
         filteredIntervals.push_back(intervals[i]);
      }
   }
   cout << "filtered size " << filteredIntervals.size();

   // cout<<intervalCountUpperBound<<" interval count "<<intervalCount<<endl;
   eleArray = (REAL *)malloc((ic * 2) * sizeof(REAL));
   eleEdgeArray = (Edge *)malloc((ic * 2) * sizeof(Edge));
   eleIntArray = (LONG *)malloc((ic * 2) * sizeof(LONG));
   sortedIndex = (int *)malloc((ic * 2) * sizeof(int));

   REAL absMinEndPoint = abs(minInterval);
   int maxDigitCount = getDigitCount(maxInterval + absMinEndPoint);

   #ifdef DG 
   high_resolution_clock::time_point cp1, cp2, cp3, cp4, cp5, cp6;
   cp1 = high_resolution_clock::now();
   #endif

   int j;
   // fill array with start vertex only
   // #pragma acc enter data copyin(intervals[0:ic])
   // // #pragma acc enter data copyin(this)
   // #pragma acc enter data create(eleArray[0:ic*2], sortedIndex[0:ic*2], eleIntArray[0:ic*2])
   // #pragma acc parallel loop
   // #pragma omp parallel for
   for (j = 0; j < ic; j++)
   {
      int i = j * 2;
      eleArray[i] = intervals[j].start;
      eleArray[i + 1] = intervals[j].end;
      sortedIndex[i] = i;
      sortedIndex[i + 1] = i + 1;
      eleIntArray[i] = realToInt(eleArray[i] + absMinEndPoint, (MAX_INTERVAL_DIGIT_COUNT - maxDigitCount));
      eleIntArray[i + 1] = realToInt(eleArray[i + 1] + absMinEndPoint, (MAX_INTERVAL_DIGIT_COUNT - maxDigitCount));
      // cout<<"--- "<<eleArray[i]<<" "<<eleArray[i+1]<<" "<<absMinEndPoint<<" "<<sortedIndex[i]<<" "<<sortedIndex[i+1]<<endl;
   }

   #ifdef DG 
   cp2 = high_resolution_clock::now();
   #endif

   // #pragma acc update self(eleArray[0:ic*2], sortedIndex[0:ic*2], eleIntArray[0:ic*2])
   // cout<<"111111111"<<endl;
   // sort elementary array
   radixsort(eleIntArray, sortedIndex, 0, ic * 2, ic * 2, MAX_INTERVAL_DIGIT_COUNT);

   #ifdef DG 
   cp3 = high_resolution_clock::now();
   #endif

   // for(int i=0; i<ic; ++i){
   // cout<<"++ "<<sortedIndex[i]<<" "<<eleIntArray[sortedIndex[i]]<<" "<<eleArray[sortedIndex[i]]<<endl;
   // }
   // cout<<"sort below"<<endl;

   REAL *elementaryArray2, *elementaryEdgeArray2;
   elementaryArray2 = (REAL *)malloc((icupper + 1) * sizeof(REAL));
   // elementaryEdgeArray=(Edge *)malloc((icupper)*sizeof(Edge));
   elementaryEdgeArray2 = (REAL *)malloc((icupper * 2) * sizeof(REAL));

   // remove duplicates and shrink elemnetary array
   currEle = eleArray[sortedIndex[0]];
   for (int i = 1; i < ic * 2; ++i)
   {
      prevCount = eleCount;
      if (currEle != eleArray[sortedIndex[i]])
      {
         elementaryArray2[eleCount] = currEle;
         eleCount += 1;
         currEle = eleArray[sortedIndex[i]];
         // cout<<i<<" ***"<<eleCount<<" "<<sortedIndex[i]<<" "<<endl;
      }
      else if (intervals[sortedIndex[i] / 2].start == intervals[sortedIndex[i] / 2].end && (lastVertical == -1 || intervals[lastVertical].start != intervals[sortedIndex[i] / 2].start))
      {
         elementaryArray2[eleCount] = intervals[sortedIndex[i] / 2].start;
         eleCount += 1;
         lastVertical = sortedIndex[i] / 2;
      }
      if (prevCount < eleCount && eleCount > 1)
      {
         // Edge e(elementaryArray2[eleCount-2], elementaryArray2[eleCount-1]);
         // elementaryEdgeArray[edgeCount]=e;
         elementaryEdgeArray2[edgeCount] = elementaryArray2[eleCount - 2];
         elementaryEdgeArray2[edgeCount + 1] = elementaryArray2[eleCount - 1];
         edgeCount += 2;
      }
   }

   #ifdef DG 
   cp4= high_resolution_clock::now();
   #endif

   // adding last element
   elementaryArray2[eleCount++] = eleArray[sortedIndex[ic * 2 - 1]];
   // Edge e(elementaryArray2[eleCount-2], elementaryArray2[eleCount-1]);
   // elementaryEdgeArray[edgeCount]=e;
   elementaryEdgeArray2[edgeCount] = elementaryArray2[eleCount - 2];
   elementaryEdgeArray2[edgeCount + 1] = elementaryArray2[eleCount - 1];
   edgeCount += 2;

   // high_resolution_clock::time_point cp1, cp2;
   // if(DEBUG) cp1=high_resolution_clock::now();

   int upperBound = Util::getNPowOf2(eleCount);
   // normalize : Make the number of edges power of 2 for convenience
   int last = eleCount - 1;
   // Edge efinal(elementaryArray2[last], elementaryArray2[last]);

   // #pragma omp parallel for
   for (j = eleCount; j <= upperBound; j++)
   {
      elementaryArray2[j] = elementaryArray2[last];
      // elementaryEdgeArray[edgeCount]=efinal;
      elementaryEdgeArray2[edgeCount] = elementaryArray2[last];
      elementaryEdgeArray2[edgeCount + 1] = elementaryArray2[last];
      // #pragma omp atomic
      edgeCount += 2;
   }

   #ifdef DG 
   cp5 = high_resolution_clock::now();
   #endif
   eleCount = j;
   elementaryPointCount = eleCount;

   elementaryArray = (REAL *)malloc((eleCount) * sizeof(REAL));
   elementaryEdgeArray = (REAL *)malloc((edgeCount) * sizeof(REAL));
   for (int ii = 0; ii < eleCount; ++ii)
   {
      elementaryArray[ii] = elementaryArray2[ii];
   }
   // #pragma acc enter data create(this, elementaryEdgeArray[0:edgeCount])
   for (int ii = 0; ii < edgeCount; ++ii)
   {
      elementaryEdgeArray[ii] = elementaryEdgeArray2[ii];
   }
   // #pragma acc enter data copyin(elementaryArray[0:eleCount])
   // #pragma acc enter data copyin(elementaryEdgeArray[0:edgeCount])

   numLeaves = edgeCount / 2;
   // #pragma acc enter data copyin(elementaryArray[0:eleCount], elementaryEdgeArray[0:edgeCount]) async

   #ifdef DG 
   cp6 = high_resolution_clock::now();
   #endif
   #ifdef DG
      // Calculating total time taken by the program.
      auto d1 = duration_cast<microseconds>(cp2 - cp1);
      auto d2 = duration_cast<microseconds>(cp3 - cp2);
      auto d3 = duration_cast<microseconds>(cp4 - cp3);
      auto d4 = duration_cast<microseconds>(cp5 - cp4);
      auto d5 = duration_cast<microseconds>(cp6 - cp5);

      cout << "CEIA: copy() " << fixed << d1.count() << setprecision(10) << endl;
      cout << "CEIA: radix sort() " << fixed << d2.count() << setprecision(10) << endl;
      cout << "CEIA: dup remove() " << fixed << d3.count() << setprecision(10) << endl;
      cout << "CEIA: last elements handling() " << fixed << d4.count() << setprecision(10) << endl;
      cout << "CEIA: to GPU copy() " << fixed << d5.count() << setprecision(10) << endl;
   #endif
}

/*
construct elementary intervals in the GPU
*/
void SegTree::constructElementaryIntervalsArrayGPUParallel(Edge *intervals)
{
   int eleCount = 0, edgeCount = 0, prevCount;
   REAL currEle, *eleArray;
   Edge *eleEdgeArray;
   LONG *eleIntArray;
   int *sortedIndex, lastVertical = -1, ic = intervalCount;
   int icupper = Util::getNPowOf2(ic);
   intervalCountUpperBound = icupper;

   // cout<<intervalCountUpperBound<<" interval count "<<intervalCount<<endl;
   eleArray = (REAL *)malloc((ic * 2) * sizeof(REAL));
   eleEdgeArray = (Edge *)malloc((ic * 2) * sizeof(Edge));
   eleIntArray = (LONG *)malloc((ic * 2) * sizeof(LONG));
   sortedIndex = (int *)malloc((ic * 2) * sizeof(int));

   REAL absMinEndPoint = abs(minInterval);
   int maxDigitCount = getDigitCount(maxInterval + absMinEndPoint);

   #ifdef DG
   high_resolution_clock::time_point cp1, cp2, cp3, cp4, cp5, cp6;
   cp1 = high_resolution_clock::now();
   #endif

   int j;
   // fill array with start vertex only
   // #pragma acc enter data copyin(intervals[0:ic])
   // // #pragma acc enter data copyin(this)
   // #pragma acc enter data create(eleArray[0:ic*2], sortedIndex[0:ic*2], eleIntArray[0:ic*2])
   // #pragma acc parallel loop
   // #pragma omp parallel for
   for (j = 0; j < ic; j++)
   {
      int i = j * 2;
      eleArray[i] = intervals[j].start;
      eleArray[i + 1] = intervals[j].end;
      sortedIndex[i] = i;
      sortedIndex[i + 1] = i + 1;
      eleIntArray[i] = realToInt(eleArray[i] + absMinEndPoint, (MAX_INTERVAL_DIGIT_COUNT - maxDigitCount));
      eleIntArray[i + 1] = realToInt(eleArray[i + 1] + absMinEndPoint, (MAX_INTERVAL_DIGIT_COUNT - maxDigitCount));
      // cout<<"--- "<<eleArray[i]<<" "<<eleArray[i+1]<<" "<<absMinEndPoint<<" "<<sortedIndex[i]<<" "<<sortedIndex[i+1]<<endl;
   }
   #ifdef DG
   cp2 = high_resolution_clock::now();
   #endif

   radixsort(eleIntArray, sortedIndex, 0, ic * 2, ic * 2, MAX_INTERVAL_DIGIT_COUNT);

   #ifdef DG
   cp3 = high_resolution_clock::now();
   #endif

   REAL *elementaryArray2, *elementaryEdgeArray2;
   elementaryArray2 = (REAL *)malloc((icupper + 1) * sizeof(REAL));
   // elementaryEdgeArray=(Edge *)malloc((icupper)*sizeof(Edge));
   elementaryEdgeArray2 = (REAL *)malloc((icupper * 2) * sizeof(REAL));

   // remove duplicates and shrink elemnetary array
   currEle = eleArray[sortedIndex[0]];
   for (int i = 1; i < ic * 2; ++i)
   {
      prevCount = eleCount;
      if (currEle != eleArray[sortedIndex[i]])
      {
         elementaryArray2[eleCount] = currEle;
         eleCount += 1;
         currEle = eleArray[sortedIndex[i]];
         // cout<<i<<" ***"<<eleCount<<" "<<sortedIndex[i]<<" "<<endl;
      }
      else if (intervals[sortedIndex[i] / 2].start == intervals[sortedIndex[i] / 2].end && (lastVertical == -1 || intervals[lastVertical].start != intervals[sortedIndex[i] / 2].start))
      {
         elementaryArray2[eleCount] = intervals[sortedIndex[i] / 2].start;
         eleCount += 1;
         lastVertical = sortedIndex[i] / 2;
      }
      if (prevCount < eleCount && eleCount > 1)
      {
         // Edge e(elementaryArray2[eleCount-2], elementaryArray2[eleCount-1]);
         // elementaryEdgeArray[edgeCount]=e;
         elementaryEdgeArray2[edgeCount] = elementaryArray2[eleCount - 2];
         elementaryEdgeArray2[edgeCount + 1] = elementaryArray2[eleCount - 1];
         edgeCount += 2;
      }
   }

   #ifdef DG
   cp4 = high_resolution_clock::now();
   #endif
   // adding last element
   elementaryArray2[eleCount++] = eleArray[sortedIndex[ic * 2 - 1]];
   // Edge e(elementaryArray2[eleCount-2], elementaryArray2[eleCount-1]);
   // elementaryEdgeArray[edgeCount]=e;
   elementaryEdgeArray2[edgeCount] = elementaryArray2[eleCount - 2];
   elementaryEdgeArray2[edgeCount + 1] = elementaryArray2[eleCount - 1];
   edgeCount += 2;

   // high_resolution_clock::time_point cp1, cp2;
   // if(DEBUG_TIME) cp1=high_resolution_clock::now();

   int upperBound = Util::getNPowOf2(eleCount);
   // normalize : Make the number of edges power of 2 for convenience
   int last = eleCount - 1;
   // Edge efinal(elementaryArray2[last], elementaryArray2[last]);

   // #pragma omp parallel for
   for (j = eleCount; j <= upperBound; j++)
   {
      elementaryArray2[j] = elementaryArray2[last];
      // elementaryEdgeArray[edgeCount]=efinal;
      elementaryEdgeArray2[edgeCount] = elementaryArray2[last];
      elementaryEdgeArray2[edgeCount + 1] = elementaryArray2[last];
      // #pragma omp atomic
      edgeCount += 2;
   }

   #ifdef DG
   cp5 = high_resolution_clock::now();
   #endif
   eleCount = j;
   elementaryPointCount = eleCount;

   elementaryArray = (REAL *)malloc((eleCount) * sizeof(REAL));
   elementaryEdgeArray = (REAL *)malloc((edgeCount) * sizeof(REAL));
   for (int ii = 0; ii < eleCount; ++ii)
   {
      elementaryArray[ii] = elementaryArray2[ii];
   }
   // #pragma acc enter data create(this, elementaryEdgeArray[0:edgeCount])
   for (int ii = 0; ii < edgeCount; ++ii)
   {
      elementaryEdgeArray[ii] = elementaryEdgeArray2[ii];
   }
   // #pragma acc enter data copyin(elementaryArray[0:eleCount])
   // #pragma acc enter data copyin(elementaryEdgeArray[0:edgeCount])

   numLeaves = edgeCount / 2;
   // #pragma acc enter data copyin(elementaryArray[0:eleCount], elementaryEdgeArray[0:edgeCount]) async

   #ifdef DG
      cp6 = high_resolution_clock::now();
   #endif
   #ifdef DT_INFO
      // Calculating total time taken by the program.
      auto d1 = duration_cast<microseconds>(cp2 - cp1);
      auto d2 = duration_cast<microseconds>(cp3 - cp2);
      auto d3 = duration_cast<microseconds>(cp4 - cp3);
      auto d4 = duration_cast<microseconds>(cp5 - cp4);
      auto d5 = duration_cast<microseconds>(cp6 - cp5);

      cout << "CEIA: copy() " << fixed << d1.count() << setprecision(10) << endl;
      cout << "CEIA: radix sort() " << fixed << d2.count() << setprecision(10) << endl;
      cout << "CEIA: dup remove() " << fixed << d3.count() << setprecision(10) << endl;
      cout << "CEIA: last elements handling() " << fixed << d4.count() << setprecision(10) << endl;
      cout << "CEIA: to GPU copy() " << fixed << d5.count() << setprecision(10) << endl;
   #endif
}

/*
method to save elementary intervals given the intervals in the GPU
*/
#if GPU!=0
void SegTree::constructElementaryIntervalsArrayGPU()
{
   int eleCount = 0, edgeCount = 0, prevCount;
   REAL currEle, *eleArray, *eleArray2, *eleIntArrayReal;
   Edge *eleEdgeArray;
   // LONG *eleIntArray;
   int *sortedIndex, lastVertical = -1, ic = intervalCount;
   int icupper = Util::getNPowOf2(ic);
   intervalCountUpperBound = icupper;

   // intervals2 = (Edge *)malloc(ic*sizeof(Edge));
   // /*int*/ ic=intervalCount;
   // for(int i=0; i<ic; ++i){
   //    intervals2[i]=intervals2[i];
   // }
   // int ic2=intervalCount;
   // #pragma acc enter data copyin(this)
   // #pragma acc enter data copyin(intervals2[0:ic2])
   // cout<<"intrevals2 in gpu"<<endl;


   // cout<<intervalCountUpperBound<<" interval count "<<intervalCount<<endl;
   eleArray = (REAL *)malloc((ic * 2) * sizeof(REAL));
   eleArray2 = (REAL *)malloc((ic * 2) * sizeof(REAL));
   eleEdgeArray = (Edge *)malloc((ic * 2) * sizeof(Edge));
   // eleIntArray = (LONG *)malloc((ic * 2) * sizeof(LONG));
   sortedIndex = (int *)malloc((ic * 2) * sizeof(int));

   // eleArray = (REAL *)malloc((ic) * sizeof(REAL));
   // eleArray2 = (REAL *)malloc((ic) * sizeof(REAL));
   // eleEdgeArray = (Edge *)malloc((ic*2) * sizeof(Edge));
   // sortedIndex = (int *)malloc((ic) * sizeof(int));
   // sortedIndex1 = (int *)malloc((ic) * sizeof(int));

   REAL absMinEndPoint = abs(minInterval);
   int maxDigitCount = getDigitCount(maxInterval + absMinEndPoint);


   #ifdef TIME_L2 
   high_resolution_clock::time_point cp1, cp2, cp3, cp4, cp5, cp6, cp7, cp8, cp9, cp10, cp11, cp12;
   #endif
   int j;

   // // first round sort by end
   // #pragma acc enter data copyin(intervals2[0:ic])
   // #pragma acc enter data create(sortedIndex1[0:ic],  eleArray1[0:ic])
   // #pragma acc parallel loop
   // for (j = 0; j < ic; j++)
   // {
   //    eleArray1[j] = intervals2[j].end;
   //    sortedIndex1[j] = j;
   // }

   // #pragma acc host_data use_device(eleArray1, sortedIndex1)
   // thrustDeviceStableSortByKey(eleArray1, sortedIndex1, ic);
   

   // fill array with start vertex only
   #ifdef TIME_L2
   cp1 = high_resolution_clock::now();
   #endif
   // #pragma acc enter data copyin(intervals[0:ic])
   #pragma acc enter data create(eleArray[0:ic*2], sortedIndex[0:ic*2], eleArray2[0:ic*2])
   // #pragma acc enter data create(eleArray[0:ic+1], sortedIndex[0:ic+1], eleArray2[0:ic+1])
   // #pragma acc data present(sortedIndex1[0:ic],  eleArray1[0:ic])
   #pragma acc data present(intervals2[0:ic])
   #pragma acc parallel loop
   for (j = 0; j < ic; j++)
   {
      int i = j * 2;
      eleArray[i] = intervals2[j].start;
      eleArray[i + 1] = intervals2[j].end;
      sortedIndex[i] = i;
      sortedIndex[i + 1] = i + 1;
      // eleIntArray[i] = realToInt(eleArray[i] + absMinEndPoint, (MAX_INTERVAL_DIGIT_COUNT - maxDigitCount));
      // eleIntArray[i + 1] = realToInt(eleArray[i + 1] + absMinEndPoint, (MAX_INTERVAL_DIGIT_COUNT - maxDigitCount));
      eleArray2[i] = eleArray[i];
      eleArray2[i + 1] = eleArray[i + 1];

      // eleArray[j] = intervals2[j].end;
      // sortedIndex[j] = j;
      // eleArray2[j] = eleArray[j];

      // if(j==ic-1){
      //    eleArray[j+1] = intervals2[j].end;
      //    sortedIndex[j+1] = j+1;
      //    eleArray2[j+1] = eleArray[j+1];
      // }
      // cout<<"--- "<<eleArray[i]<<" "<<eleArray[i+1]<<" "<<absMinEndPoint<<" "<<sortedIndex[i]<<" "<<sortedIndex[i+1]<<endl;
   }
   #ifdef TIME_L2 
   cp2 = high_resolution_clock::now();
   #endif
   // #pragma acc update self(eleArray[0:ic*2], sortedIndex[0:ic*2], eleIntArray[0:ic*2])
   // sort elementary array

   // quickSortGetIndexArray(eleIntArrayReal, sortedIndex, 0, ic*2);
   // radixsort(eleIntArray, sortedIndex, 0, ic * 2, ic * 2, MAX_INTERVAL_DIGIT_COUNT);

   // #pragma acc data present(eleIntArray[0:ic*2], sortedIndex[0:ic*2])
   #pragma acc host_data use_device(eleArray2, sortedIndex)
   thrustDeviceStableSortByKey(eleArray2, sortedIndex, ic*2);  
   // #pragma acc update self(eleArray[0:ic*2], sortedIndex[0:ic*2])
   // thrustDeviceStableSortByKey(eleArray2, sortedIndex, ic+1);
   // #pragma acc update self(eleArray[0:ic], sortedIndex[0:ic])
   #ifdef TIME_L2
   cp3 = high_resolution_clock::now();
   #endif

   int *elementaryArrayIndex;
   elementaryArrayIndex = (int *)malloc((ic*2+1)*sizeof(int));  //******** use calloc in CPU
   // elementaryArrayIndex = (int *)calloc(ic+2, sizeof(int));
   int *verticalLines;
   verticalLines = (int *)malloc((ic*2+1)*sizeof(int)); //******** use calloc in CPU

   // mark vertical lines
   #pragma acc data present(sortedIndex[0:ic*2], intervals2[0:ic])
   #pragma acc enter data create(verticalLines[0:ic*2+1])
   #pragma acc parallel loop
   for(int j=0; j<ic; ++j){
      int i = sortedIndex[j]/2;
      verticalLines[i]=0;
      if(i==ic-1) verticalLines[i+1]=0;
      if(intervals2[i].start == intervals2[i].end){
         verticalLines[i]=1;
      }
   }
   #ifdef TIME_L2
   cp4 = high_resolution_clock::now();
   #endif

   // thrust::exclusive_scan(verticalLines, verticalLines+(ic*2+1), verticalLines);
   #pragma acc host_data use_device(verticalLines)
   thrustDeviceExclusiveScan(verticalLines, (ic*2+1));  
   #pragma acc update self(verticalLines[ic*2])
   #ifdef TIME_L2
   cp5 = high_resolution_clock::now();
   #endif
   // compact vertical lines into new array
   int *verticalLineLoc, numVLines=verticalLines[ic*2];
   // cout<<"Num vertical lines="<<numVLines<<endl;
   verticalLineLoc = (int *)malloc(numVLines*sizeof(int));

   #pragma acc data present(verticalLines[0:ic*2+1])
   #pragma acc enter data create(verticalLineLoc[0:numVLines])
   #pragma acc parallel loop
   for(int i=0; i<ic; ++i){
      if(verticalLines[i] != verticalLines[i+1]){
         verticalLineLoc[verticalLines[i]]=i;
      }
   }

   #ifdef TIME_L2 
   cp6 = high_resolution_clock::now();
   #endif
   // for(int i=0; i<ic*2; ++i){
   //    cout<<"++ i="<<i<<" "<<sortedIndex[i]<<" "<<eleArray[sortedIndex[i]]<<endl;
   // }
   // cout<<"sort below"<<endl;
   

   REAL *elementaryArray2, *elementaryEdgeArray2;
   REAL *elementaryArray3;

   // for(int i=0; i<ic*2-1; ++i){
   //    if(intervals2[sortedIndex[i]/2].end!=intervals2[sortedIndex[i+1]/2].start){
   //       cout<<"inequal i="<<i<<" "<<intervals2[sortedIndex[i]/2].start<<" "<<intervals2[sortedIndex[i]/2].end<<" "<<intervals2[sortedIndex[i+1]/2].start<<" "<<intervals2[sortedIndex[i+1]/2].end<<endl;
   //       break;
   //    }
   // }


   // cout<<"ic="<<ic<<" dup removal icupper="<<icupper<<endl;

   // /*
   elementaryArray3 = (REAL *)malloc(ic*2 * sizeof(REAL));
   // elementaryArray3 = (REAL *)malloc(ic+1 * sizeof(REAL));
   // elementaryEdgeArray2 = (REAL *)malloc(ic*2 * sizeof(REAL));

   // #pragma acc update self(eleArray[0:ic*2], sortedIndex[0:ic*2], verticalLines[0:ic*2+1], verticalLineLoc[0:numVLines])

   // remove duplicates and shrink elemnetary array
   #ifdef TIME_L2
   cp7 = high_resolution_clock::now();
   #endif
   #pragma acc data present(sortedIndex[0:ic*2], intervals2[0:ic], verticalLines[0:ic*2+1], verticalLineLoc[0:numVLines], eleArray[0:ic*2])
   #pragma acc enter data create(elementaryArray3[0:ic*2], elementaryArrayIndex[0:ic*2+1])
   #pragma acc parallel loop
   for(int i=0; i<ic*2-1; ++i){
   // for(int i=0; i<ic; ++i){
      // cout<<i<<" "; 
      // int intvi = sortedIndex[sortedIndex1[i]];
      // int intvi2 = sortedIndex[sortedIndex1[i-1]];
      elementaryArrayIndex[i]=0;
      if(i==ic*2-2) {
         elementaryArrayIndex[i+1]=0;
         elementaryArrayIndex[i+2]=0;
      }

      int intvi = sortedIndex[i]/2;
      // int intvi2 = sortedIndex[i-1]/2;
      int comid = verticalLines[i/2];
      int lookup;
      if (comid>0) lookup = verticalLineLoc[comid-1];

      if(eleArray[sortedIndex[i]]!=eleArray[sortedIndex[i+1]]){
         // cout<<"case 1 i="<<i<<" "<<eleArray[sortedIndex[i]]<<endl;
         elementaryArray3[i]=eleArray[sortedIndex[i]];
         elementaryArrayIndex[i]=1;
      }
      // else if (intervals2[sortedIndex[i]/2].start == intervals2[sortedIndex[i]/2].end && (i==0 || (intervals2[sortedIndex[i-1]/2].start != intervals2[sortedIndex[i]/2].start && intervals2[sortedIndex[i-1]/2].end != intervals2[sortedIndex[i]/2].start))){
      else if (intervals2[intvi].start == intervals2[intvi].end && (comid==0 || (intervals2[lookup].start != intervals2[intvi].start && intervals2[lookup].end != intervals2[intvi].start))){
         // cout<<"case 2 i="<<i<<" "<<intervals2[sortedIndex[i]/2].start<<endl;
         elementaryArray3[i]=intervals2[intvi].start;
         // elementaryArray3[i]=intervals2[sortedIndex[i]/2].start;
         elementaryArrayIndex[i]=1;
      }

      if(i==ic*2-2) {
         if(eleArray[sortedIndex[ic*2-2]]!=eleArray[sortedIndex[ic*2-1]]){
            elementaryArray3[ic*2-1]=eleArray[sortedIndex[ic*2-2]];
            elementaryArrayIndex[ic*2-1]=1;
         }
      }
   }
   #ifdef TIME_L2
   cp8 = high_resolution_clock::now();
   #endif

   // #pragma acc update self(elementaryArrayIndex[0:ic*2+1], elementaryArray3[0:ic*2])
   // if(eleArray[sortedIndex[ic*2-2]]!=eleArray[sortedIndex[ic*2-1]]){
   //    elementaryArray3[ic*2-1]=eleArray[sortedIndex[ic*2-2]];
   //    elementaryArrayIndex[ic*2-1]=1;
   // }
   // if(eleArray[sortedIndex[ic-1]]!=eleArray[sortedIndex[ic]]){
   //    elementaryArray3[ic]=eleArray[sortedIndex[ic-1]];
   //    elementaryArrayIndex[ic]=1;
   // }


   // cout<<"before prefix"<<endl;

   // thrust::exclusive_scan(elementaryArrayIndex, elementaryArrayIndex+(ic*2+1), elementaryArrayIndex);
   #pragma acc host_data use_device(elementaryArrayIndex)
   thrustDeviceExclusiveScan(elementaryArrayIndex, (ic*2+1));  
   // thrust::exclusive_scan(elementaryArrayIndex, elementaryArrayIndex+(ic+2), elementaryArrayIndex);
   #ifdef TIME_L2
   cp9 = high_resolution_clock::now();
   #endif
   #pragma  acc update self(elementaryArrayIndex[ic*2])
   #ifdef TIME_L2
   cp10 = high_resolution_clock::now();
   #endif
   // #pragma  acc update self(elementaryArrayIndex[0:ic*2+1], elementaryArray3[0:ic*2])

   eleCount = elementaryArrayIndex[ic*2];
   // eleCount = elementaryArrayIndex[ic+1];
   // REAL *elementaryArray11, *elementaryEdgeArray11;
   // cout<<"***after prefix eleCount="<<eleCount<<endl;

   int eleCountRound = Util::getNPowOf2(eleCount)+1;
   elementaryArray = (REAL *)malloc((eleCountRound) * sizeof(REAL));
   #pragma acc enter data copyin(this)
   #pragma acc enter data create(elementaryArray[0:eleCountRound])

   #pragma acc data present(elementaryArrayIndex[0:ic*2+1], elementaryArray3[0:ic*2], elementaryArray[0:eleCountRound])
   #pragma acc parallel loop
   for(int i=0; i<ic*2; ++i){
   // for(int i=0; i<ic+1; ++i){
      if(elementaryArrayIndex[i]!=elementaryArrayIndex[i+1]){
         elementaryArray[elementaryArrayIndex[i]]=elementaryArray3[i];
         // cout<<"c1 i="<<i<<" at="<<elementaryArrayIndex[i]<<" val="<<elementaryArray3[i]<<endl;
      }
      // if(elementaryArrayIndex[i]>eleCountRound) cout<<"erro i="<<i<<" "<<elementaryArrayIndex[i]<<endl; 
   }

   // cout<<"init copy eleCountRound="<<eleCountRound<<endl;

   // #pragma acc update self(elementaryArray[0:eleCountRound])

   #pragma acc data present(elementaryArray[0:eleCountRound])
   #pragma acc parallel loop
   for(int i=0; i<(eleCountRound-eleCount); ++i){
      elementaryArray[i+eleCount]=elementaryArray[eleCount-1];
   }
   #ifdef TIME_L2
   cp11 = high_resolution_clock::now();
   #endif

   /*// --------------------------------------
   // int *elementaryArray2;
   elementaryArray2 = (REAL *)malloc((upperBound) * sizeof(REAL));
   // #pragma acc enter data copyin(this)
   // #pragma acc enter data create(elementaryArray2[0:upperBound])

   // #pragma acc data present(elementaryArrayIndex[0:ic*2+1], elementaryArray3[0:ic*2])
   // #pragma acc parallel loop
   for(int i=0; i<ic*2; ++i){
   // for(int i=0; i<ic+1; ++i){
      if(elementaryArrayIndex[i]!=elementaryArrayIndex[i+1]){
         elementaryArray2[elementaryArrayIndex[i]]=elementaryArray3[i];
         // cout<<"c1 i="<<i<<" at="<<elementaryArrayIndex[i]<<" val="<<elementaryArray3[i]<<endl;
      }
   }

   // cout<<"init copy upperBound="<<upperBound<<endl;

   // #pragma acc update self(elementaryArray[0:upperBound])

   // #pragma acc data present(elementaryArray[0:upperBound])
   // #pragma acc parallel loop
   for(int i=0; i<(upperBound-eleCount); ++i){
      elementaryArray2[i+eleCount]=elementaryArray2[eleCount-1];
   }
   // --------------------------------------*/


   // #pragma acc update self(elementaryArray[0:eleCountRound])
   /*// --------------------------------------

   // cout<<"start "<<elementaryArray[0]<<" "<<elementaryArray2[0]<<endl;
   // for (int i = 0; i < upperBound; i++)
   // {
   //    if(elementaryArray[i]!=elementaryArray2[i]){
   //       cout<<"i="<<i<<" val1="<<elementaryArray[i]<<" val2="<<elementaryArray2[i]<<endl;
   //    }
   // }
   // cout<<"end=== "<<endl;
   // --------------------------------------*/
   

   // #pragma acc enter data copyin(elementaryArray[0:eleCountRound])

   eleCount = eleCountRound;
   // cout<<"eleCount="<<eleCount<<" eleCountRound="<<eleCountRound<<endl;

   edgeCount = eleCountRound*2-2;
   elementaryEdgeArray = (REAL *)malloc(edgeCount * sizeof(REAL));



   #pragma acc enter data copyin(this)
   #pragma acc enter data create(elementaryEdgeArray[0:edgeCount])
   #pragma acc data present(elementaryArray[0:eleCountRound])
   #pragma acc parallel loop
   for(int i=0; i<eleCount-1; ++i){
      elementaryEdgeArray[i*2]=elementaryArray[i];
      elementaryEdgeArray[i*2+1]=elementaryArray[i+1];
      // cout<<"* "<<i<<" "<<elementaryEdgeArray[i*2]<<" "<<elementaryEdgeArray[i*2+1]<<endl;
   }
   #ifdef TIME_L2
   cp12 = high_resolution_clock::now();
   #endif
   // cout<<"cone"<<endl;

   // #pragma acc update self(elementaryArray[0:eleCountRound], elementaryEdgeArray[0:eleCountRound*2-2])

   elementaryPointCount = eleCount;
   numLeaves = eleCount-1;
   int ipc=numLeaves;
   // eleCount=0;
   // edgeCount=0;
   
   // */

   /*
   elementaryArray2 = (REAL *)malloc((icupper + 1) * sizeof(REAL));
   elementaryEdgeArray2 = (REAL *)malloc((icupper * 2) * sizeof(REAL));
   // remove duplicates and shrink elemnetary array
   currEle = eleArray[sortedIndex[0]];
   // for (int i = 1; i < ic * 2; ++i)
   for (int i = 1; i < ic; ++i)
   {
      prevCount = eleCount;
      if (currEle != eleArray[sortedIndex[i]])
      {
         elementaryArray2[eleCount] = currEle;
         eleCount += 1;
         currEle = eleArray[sortedIndex[i]];
         // cout<<i<<" ***"<<eleCount<<" "<<sortedIndex[i]<<" "<<endl;
      }
      // else if (intervals[sortedIndex[i] / 2].start == intervals[sortedIndex[i] / 2].end && (lastVertical == -1 || intervals[lastVertical].start != intervals[sortedIndex[i] / 2].start))
      else if (intervals[sortedIndex[i]].start == intervals[sortedIndex[i]].end && (lastVertical == -1 || intervals[lastVertical].start != intervals[sortedIndex[i]].start))
      // else if (intervals[sortedIndex[i] / 2].start == intervals[sortedIndex[i] / 2].end && (lastVertical == -1 || intervals[sortedIndex[i-1] / 2].start != intervals[sortedIndex[i] / 2].start)) // values fomr here in the txt
      // else if (intervals[sortedIndex[i] / 2].start == intervals[sortedIndex[i] / 2].end && (i == 1 || ((intervals[sortedIndex[i-1] / 2].start != intervals[sortedIndex[i] / 2].start && intervals[sortedIndex[i-1] / 2].start != intervals[sortedIndex[i-1] / 2].end))))
      // if (intervals[sortedIndex[i]].start == intervals[sortedIndex[i]].end && (i == 1 || (intervals[sortedIndex[i-1]].start != intervals[sortedIndex[i]].start && intervals[sortedIndex[i-1]].end != intervals[sortedIndex[i]].start)))
      {
         // elementaryArray2[eleCount] = intervals[sortedIndex[i] / 2].start;
         elementaryArray2[eleCount] = intervals[sortedIndex[i]].start;
         eleCount += 1;
         // lastVertical = sortedIndex[i] / 2;
         lastVertical = sortedIndex[i];
      }
      if (prevCount < eleCount && eleCount > 1)
      {
         // Edge e(elementaryArray2[eleCount-2], elementaryArray2[eleCount-1]);
         // elementaryEdgeArray[edgeCount]=e;
         elementaryEdgeArray2[edgeCount] = elementaryArray2[eleCount - 2];
         elementaryEdgeArray2[edgeCount + 1] = elementaryArray2[eleCount - 1];
         edgeCount += 2;
      }
   }
   cout<<"eleCout="<<eleCount<<endl;
   #ifdef DG 
   cp4 = high_resolution_clock::now();
   #endif
   // adding last element
   // elementaryArray2[eleCount++] = eleArray[sortedIndex[ic * 2 - 1]];
   elementaryArray2[eleCount++] = eleArray[sortedIndex[ic]];
   // Edge e(elementaryArray2[eleCount-2], elementaryArray2[eleCount-1]);
   // elementaryEdgeArray[edgeCount]=e;
   elementaryEdgeArray2[edgeCount] = elementaryArray2[eleCount - 2];
   elementaryEdgeArray2[edgeCount + 1] = elementaryArray2[eleCount - 1];
   edgeCount += 2;
   
   cout<<"eleCount="<<eleCount<<" edgeCount="<<edgeCount<<endl;

   // high_resolution_clock::time_point cp1, cp2;
   // if(DEBUG) cp1=high_resolution_clock::now();

   int upperBound = Util::getNPowOf2(eleCount);
   // normalize : Make the number of edges power of 2 for convenience
   int last = eleCount - 1;
   // Edge efinal(elementaryArray2[last], elementaryArray2[last]);

   // #pragma omp parallel for
   for (j = eleCount; j <= upperBound; j++)
   {
      elementaryArray2[j] = elementaryArray2[last];
      // elementaryEdgeArray[edgeCount]=efinal;
      elementaryEdgeArray2[edgeCount] = elementaryArray2[last];
      elementaryEdgeArray2[edgeCount + 1] = elementaryArray2[last];
      // #pragma omp atomic
      edgeCount += 2;
   }

   cout<<"eleCount="<<eleCount<<" edgeCount="<<edgeCount<<endl;

   #ifdef DG 
   cp5 = high_resolution_clock::now();
   #endif
   eleCount = j;
   elementaryPointCount = eleCount;

   cout<<"eleCount="<<eleCount<<" edgeCount="<<edgeCount<<endl;

   elementaryArray = (REAL *)malloc((eleCount) * sizeof(REAL));
   elementaryEdgeArray = (REAL *)malloc((edgeCount) * sizeof(REAL));
   for (int ii = 0; ii < eleCount; ++ii)
   {
      elementaryArray[ii] = elementaryArray2[ii];
   }
   // #pragma acc enter data create(this, elementaryEdgeArray[0:edgeCount])
   for (int ii = 0; ii < edgeCount; ++ii)
   {
      elementaryEdgeArray[ii] = elementaryEdgeArray2[ii];
   }
   // #pragma acc enter data copyin(elementaryArray[0:eleCount])
   // #pragma acc enter data copyin(elementaryEdgeArray[0:edgeCount])
   
   int jj=0;
   // for(int i=0; i<eleCount; ++i){
   //       // if(i%2==0 && elementaryArray[i]!=elementaryArray[i+1]){
   //          cout<<"arrays i="<<i<<" elementaryArray[i]="<<elementaryArray[i]<<" elementaryArray11[i]="<<elementaryArray11[i]<<endl;
   //       //    break;
   //       // }
   //       // if(i>20) break;
   //    // if(elementaryArray[i]!=elementaryArray11[i]){
   //    //    cout<<"Error i="<<i<<" elementaryArray[i]="<<elementaryArray[i]<<" elementaryArray11[i]="<<elementaryArray11[i]<<endl;
   //    //    if(jj>5) break;
   //    //    jj++;
   //    // }
   // }

   // for(int i=0; i<edgeCount; ++i){
   //    cout<<"arrays i="<<i<<" elementaryEdgeArray[i]="<<elementaryArray[i]<<" elementaryEdgeArray11[i]="<<elementaryArray11[i]<<endl;
   // }

   numLeaves = edgeCount / 2; 
   // */
   // #pragma acc enter data copyin(elementaryArray[0:eleCountRound], elementaryEdgeArray[0:eleCountRound*2-2]) 
   // #pragma acc enter data copyin(this)
   // #pragma acc enter data copyin(elementaryEdgeArray[0:ipc*2])

   #ifdef DT_INFO
      // Calculating total time taken by the program.
      auto d1 = duration_cast<microseconds>(cp2 - cp1);
      auto d2 = duration_cast<microseconds>(cp3 - cp2);
      auto d3 = duration_cast<microseconds>(cp4 - cp3);
      auto d4 = duration_cast<microseconds>(cp5 - cp4);
      auto d5 = duration_cast<microseconds>(cp6 - cp5);
      auto d7 = duration_cast<microseconds>(cp8 - cp7);
      auto d8 = duration_cast<microseconds>(cp9 - cp8);
      auto d9 = duration_cast<microseconds>(cp10 - cp9);
      auto d10 = duration_cast<microseconds>(cp11 - cp10);
      auto d11 = duration_cast<microseconds>(cp12 - cp11);

      cout << " |-CEIA: copy intervals to array " << fixed << d1.count() << setprecision(10) << endl;
      cout << " |-CEIA: thrust sort " << fixed << d2.count() << setprecision(10) << endl;
      cout << " |-CEIA: vertical line counting " << fixed << d3.count() << setprecision(10) << endl;
      cout << " |-CEIA: prefix sum: vertical line counts " << fixed << d4.count() << setprecision(10) << endl;
      cout << " |-CEIA: compact vertical lines " << fixed << d5.count() << setprecision(10) << endl;
      cout << " |-CEIA: duplicate removal " << fixed << d7.count() << setprecision(10) << endl;
      cout << " |-CEIA: prefix sum: elementaryArrayIndex " << fixed << d8.count() << setprecision(10) << endl;
      cout << " |-CEIA: copy max elementaryArrayIndex to CPU " << fixed << d9.count() << setprecision(10) << endl;
      cout << " |-CEIA: save elementaryArray values " << fixed << d10.count() << setprecision(10) << endl;
      cout << " |-CEIA: save elementaryEdgeArray values " << fixed << d11.count() << setprecision(10) << endl;
   #elif TIME_L2
      // Calculating total time taken by the program.
      auto d9 = duration_cast<microseconds>(cp10 - cp9);

      cout << fixed << d9.count() << setprecision(10) << " ";
   #endif
}
#endif

void SegTree::constructElementaryIntervalsArrayMulticore(Edge *intervals)
{
   int eleCount = 0, edgeCount = 0, prevCount;
   REAL currEle, *eleArray, *eleArray2, *eleIntArrayReal;
   Edge *eleEdgeArray;
   // LONG *eleIntArray;
   int *sortedIndex, lastVertical = -1, ic = intervalCount;
   int icupper = Util::getNPowOf2(ic);
   intervalCountUpperBound = icupper;

   eleArray = (REAL *)malloc((ic * 2) * sizeof(REAL));
   eleArray2 = (REAL *)malloc((ic * 2) * sizeof(REAL));
   eleEdgeArray = (Edge *)malloc((ic * 2) * sizeof(Edge));
   // eleIntArray = (LONG *)malloc((ic * 2) * sizeof(LONG));
   sortedIndex = (int *)malloc((ic * 2) * sizeof(int));


   REAL absMinEndPoint = abs(minInterval);
   int maxDigitCount = getDigitCount(maxInterval + absMinEndPoint);


   #ifdef DG 
   high_resolution_clock::time_point cp1, cp2, cp3, cp4, cp5, cp6;
   cp1 = high_resolution_clock::now();
   #endif
   int j;

   #pragma omp parallel for
   for (j = 0; j < ic; j++)
   {
      int i = j * 2;
      eleArray[i] = intervals[j].start;
      eleArray[i + 1] = intervals[j].end;
      sortedIndex[i] = i;
      sortedIndex[i + 1] = i + 1;

      eleArray2[i] = eleArray[i];
      eleArray2[i + 1] = eleArray[i + 1];
   }
   #ifdef DG 
   cp2 = high_resolution_clock::now();
   #endif

   thrust::stable_sort_by_key(eleArray2, eleArray2+ic*2, sortedIndex);
   // thrustDeviceStableSortByKey(eleArray2, sortedIndex, ic*2);
   // radixsort(eleArray2, sortedIndex, 0, ic * 2, ic * 2, MAX_INTERVAL_DIGIT_COUNT);
   

   int *elementaryArrayIndex;
   elementaryArrayIndex = (int *)malloc((ic*2+1)*sizeof(int));  //******** use calloc in CPU
   // elementaryArrayIndex = (int *)calloc(ic+2, sizeof(int));
   int *verticalLines;
   verticalLines = (int *)malloc((ic*2+1)*sizeof(int)); //******** use calloc in CPU

   // mark vertical lines
   #pragma omp parallel for
   for(int j=0; j<ic; ++j){
      int i = sortedIndex[j]/2;
      verticalLines[i]=0;
      if(i==ic-1) verticalLines[i+1]=0;
      if(intervals[i].start == intervals[i].end){
         verticalLines[i]=1;
      }
   }

   thrust::exclusive_scan(verticalLines, verticalLines+(ic*2+1), verticalLines);

   // compact vertical lines into new array
   int *verticalLineLoc, numVLines=verticalLines[ic*2];
   // cout<<"Num vertical lines="<<numVLines<<endl;
   verticalLineLoc = (int *)malloc(numVLines*sizeof(int));

   #pragma omp parallel for
   for(int i=0; i<ic; ++i){
      if(verticalLines[i] != verticalLines[i+1]){
         verticalLineLoc[verticalLines[i]]=i;
      }
   }

   #ifdef DG 
   cp3 = high_resolution_clock::now();
   #endif
   

   REAL *elementaryArray2, *elementaryEdgeArray2;
   REAL *elementaryArray3;

   elementaryArray3 = (REAL *)malloc(ic*2 * sizeof(REAL));

   // remove duplicates and shrink elemnetary array
   #pragma omp parallel for
   for(int i=0; i<ic*2-1; ++i){
      elementaryArrayIndex[i]=0;
      if(i==ic*2-2) {
         elementaryArrayIndex[i+1]=0;
         elementaryArrayIndex[i+2]=0;
      }

      int intvi = sortedIndex[i]/2;
      // int intvi2 = sortedIndex[i-1]/2;
      int comid = verticalLines[i/2];
      int lookup;
      if (comid>0) lookup = verticalLineLoc[comid-1];

      if(eleArray[sortedIndex[i]]!=eleArray[sortedIndex[i+1]]){
         // cout<<"case 1 i="<<i<<" "<<eleArray[sortedIndex[i]]<<endl;
         elementaryArray3[i]=eleArray[sortedIndex[i]];
         elementaryArrayIndex[i]=1;
      }
      // else if (intervals[sortedIndex[i]/2].start == intervals[sortedIndex[i]/2].end && (i==0 || (intervals[sortedIndex[i-1]/2].start != intervals[sortedIndex[i]/2].start && intervals[sortedIndex[i-1]/2].end != intervals[sortedIndex[i]/2].start))){
      else if (intervals[intvi].start == intervals[intvi].end && (comid==0 || (intervals[lookup].start != intervals[intvi].start && intervals[lookup].end != intervals[intvi].start))){
         // cout<<"case 2 i="<<i<<" "<<intervals[sortedIndex[i]/2].start<<endl;
         elementaryArray3[i]=intervals[intvi].start;
         // elementaryArray3[i]=intervals[sortedIndex[i]/2].start;
         elementaryArrayIndex[i]=1;
      }

      if(i==ic*2-2) {
         if(eleArray[sortedIndex[ic*2-2]]!=eleArray[sortedIndex[ic*2-1]]){
            elementaryArray3[ic*2-1]=eleArray[sortedIndex[ic*2-2]];
            elementaryArrayIndex[ic*2-1]=1;
         }
      }
   }

   thrust::exclusive_scan(elementaryArrayIndex, elementaryArrayIndex+(ic*2+1), elementaryArrayIndex);

   eleCount = elementaryArrayIndex[ic*2];

   int eleCountRound = Util::getNPowOf2(eleCount)+1;
   elementaryArray = (REAL *)malloc((eleCountRound) * sizeof(REAL));

   #pragma omp parallel for
   for(int i=0; i<ic*2; ++i){
   // for(int i=0; i<ic+1; ++i){
      if(elementaryArrayIndex[i]!=elementaryArrayIndex[i+1]){
         elementaryArray[elementaryArrayIndex[i]]=elementaryArray3[i];
         // cout<<"c1 i="<<i<<" at="<<elementaryArrayIndex[i]<<" val="<<elementaryArray3[i]<<endl;
      }
   }

   #pragma omp parallel for
   for(int i=0; i<(eleCountRound-eleCount); ++i){
      elementaryArray[i+eleCount]=elementaryArray[eleCount-1];
   }


   eleCount = eleCountRound;

   edgeCount = eleCountRound*2-2;
   elementaryEdgeArray = (REAL *)malloc(edgeCount * sizeof(REAL));

   #pragma omp parallel for
   for(int i=0; i<eleCount-1; ++i){
      elementaryEdgeArray[i*2]=elementaryArray[i];
      elementaryEdgeArray[i*2+1]=elementaryArray[i+1];
   }

   elementaryPointCount = eleCount;
   numLeaves = eleCount-1;
   int ipc=numLeaves;

   #ifdef DG 
   cp6 = high_resolution_clock::now();
   #endif
   #ifdef DG
      // Calculating total time taken by the program.
      auto d1 = duration_cast<microseconds>(cp2 - cp1);
      auto d2 = duration_cast<microseconds>(cp3 - cp2);
      auto d3 = duration_cast<microseconds>(cp4 - cp3);
      auto d4 = duration_cast<microseconds>(cp5 - cp4);
      auto d5 = duration_cast<microseconds>(cp6 - cp5);

      cout << "CEIA: copy() " << fixed << d1.count() << setprecision(10) << endl;
      cout << "CEIA: radix sort() " << fixed << d2.count() << setprecision(10) << endl;
      cout << "CEIA: dup remove() " << fixed << d3.count() << setprecision(10) << endl;
      cout << "CEIA: last elements handling() " << fixed << d4.count() << setprecision(10) << endl;
      cout << "CEIA: to GPU copy() " << fixed << d5.count() << setprecision(10) << endl;
   #endif
}

/*
method to save elementary intervals given the intervals
*/
void SegTree::constructElementaryIntervalsArrayParallel(Edge *intervals)
{
   int eleCount = 0, edgeCount = 0, prevCount;
   REAL currEle, *eleArray, *eleIntArrayReal;
   Edge *eleEdgeArray;
   LONG *eleIntArray;
   int *sortedIndex, lastVertical = -1, ic = intervalCount;
   int icupper = Util::getNPowOf2(ic);
   intervalCountUpperBound = icupper;

   // cout<<intervalCountUpperBound<<" interval count "<<intervalCount<<endl;
   eleArray = (REAL *)malloc((ic * 2) * sizeof(REAL));
   eleEdgeArray = (Edge *)malloc((ic * 2) * sizeof(Edge));
   eleIntArray = (LONG *)malloc((ic * 2) * sizeof(LONG));
   eleIntArrayReal = (REAL *)malloc((ic * 2) * sizeof(REAL));
   sortedIndex = (int *)malloc((ic * 2) * sizeof(int));

   REAL absMinEndPoint = abs(minInterval);
   int maxDigitCount = getDigitCount(maxInterval + absMinEndPoint);

   #ifdef DG 
   high_resolution_clock::time_point cp1, cp2, cp3, cp4, cp5, cp6;
   cp1 = high_resolution_clock::now();
   #endif
   int j;
   // fill array with start vertex only
   // #pragma acc enter data copyin(intervals[0:ic])
   // // #pragma acc enter data copyin(this)
   // #pragma acc enter data create(eleArray[0:ic*2], sortedIndex[0:ic*2], eleIntArray[0:ic*2])
   // #pragma acc parallel loop
   // #pragma omp parallel for
   for (j = 0; j < ic; j++)
   {
      int i = j * 2;
      eleArray[i] = intervals[j].start;
      eleArray[i + 1] = intervals[j].end;
      sortedIndex[i] = i;
      sortedIndex[i + 1] = i + 1;
      eleIntArrayReal[i] = eleArray[i];
      eleIntArrayReal[i + 1] = eleArray[i + 1];
      eleIntArray[i] = realToInt(eleArray[i] + absMinEndPoint, (MAX_INTERVAL_DIGIT_COUNT - maxDigitCount));
      eleIntArray[i + 1] = realToInt(eleArray[i + 1] + absMinEndPoint, (MAX_INTERVAL_DIGIT_COUNT - maxDigitCount));
      // cout<<"--- "<<eleArray[i]<<" "<<eleArray[i+1]<<" "<<absMinEndPoint<<" "<<sortedIndex[i]<<" "<<sortedIndex[i+1]<<endl;
   }
   #ifdef DG 
   cp2 = high_resolution_clock::now();
   #endif
   // #pragma acc update self(eleArray[0:ic*2], sortedIndex[0:ic*2], eleIntArray[0:ic*2])
   // sort elementary array

   // quickSortGetIndexArray(eleIntArrayReal, sortedIndex, 0, ic*2);
   radixsort(eleIntArray, sortedIndex, 0, ic * 2, ic * 2, MAX_INTERVAL_DIGIT_COUNT);

   #ifdef DG 
   cp3 = high_resolution_clock::now();
   #endif
   // for(int i=0; i<ic; ++i){
   // cout<<"++ "<<sortedIndex[i]<<" "<<eleIntArray[sortedIndex[i]]<<" "<<eleArray[sortedIndex[i]]<<endl;
   // }
   // cout<<"sort below"<<endl;

   REAL *elementaryArray2, *elementaryEdgeArray2;
   elementaryArray2 = (REAL *)malloc((icupper + 1) * sizeof(REAL));
   // elementaryEdgeArray=(Edge *)malloc((icupper)*sizeof(Edge));
   elementaryEdgeArray2 = (REAL *)malloc((icupper * 2) * sizeof(REAL));

   // remove duplicates and shrink elemnetary array
   currEle = eleArray[sortedIndex[0]];
   for (int i = 1; i < ic * 2; ++i)
   {
      prevCount = eleCount;
      if (currEle != eleArray[sortedIndex[i]])
      {
         elementaryArray2[eleCount] = currEle;
         eleCount += 1;
         currEle = eleArray[sortedIndex[i]];
         // cout<<i<<" ***"<<eleCount<<" "<<sortedIndex[i]<<" "<<endl;
      }
      else if (intervals[sortedIndex[i] / 2].start == intervals[sortedIndex[i] / 2].end && (lastVertical == -1 || intervals[lastVertical].start != intervals[sortedIndex[i] / 2].start))
      {
         elementaryArray2[eleCount] = intervals[sortedIndex[i] / 2].start;
         eleCount += 1;
         lastVertical = sortedIndex[i] / 2;
      }
      if (prevCount < eleCount && eleCount > 1)
      {
         // Edge e(elementaryArray2[eleCount-2], elementaryArray2[eleCount-1]);
         // elementaryEdgeArray[edgeCount]=e;
         elementaryEdgeArray2[edgeCount] = elementaryArray2[eleCount - 2];
         elementaryEdgeArray2[edgeCount + 1] = elementaryArray2[eleCount - 1];
         edgeCount += 2;
      }
   }
   #ifdef DG 
   cp4 = high_resolution_clock::now();
   #endif
   // adding last element
   elementaryArray2[eleCount++] = eleArray[sortedIndex[ic * 2 - 1]];
   // Edge e(elementaryArray2[eleCount-2], elementaryArray2[eleCount-1]);
   // elementaryEdgeArray[edgeCount]=e;
   elementaryEdgeArray2[edgeCount] = elementaryArray2[eleCount - 2];
   elementaryEdgeArray2[edgeCount + 1] = elementaryArray2[eleCount - 1];
   edgeCount += 2;

   // high_resolution_clock::time_point cp1, cp2;
   // if(DEBUG) cp1=high_resolution_clock::now();

   int upperBound = Util::getNPowOf2(eleCount);
   // normalize : Make the number of edges power of 2 for convenience
   int last = eleCount - 1;
   // Edge efinal(elementaryArray2[last], elementaryArray2[last]);

   // #pragma omp parallel for
   for (j = eleCount; j <= upperBound; j++)
   {
      elementaryArray2[j] = elementaryArray2[last];
      // elementaryEdgeArray[edgeCount]=efinal;
      elementaryEdgeArray2[edgeCount] = elementaryArray2[last];
      elementaryEdgeArray2[edgeCount + 1] = elementaryArray2[last];
      // #pragma omp atomic
      edgeCount += 2;
   }

   #ifdef DG 
   cp5 = high_resolution_clock::now();
   #endif
   eleCount = j;
   elementaryPointCount = eleCount;

   elementaryArray = (REAL *)malloc((eleCount) * sizeof(REAL));
   elementaryEdgeArray = (REAL *)malloc((edgeCount) * sizeof(REAL));
   for (int ii = 0; ii < eleCount; ++ii)
   {
      elementaryArray[ii] = elementaryArray2[ii];
   }
   // #pragma acc enter data create(this, elementaryEdgeArray[0:edgeCount])
   for (int ii = 0; ii < edgeCount; ++ii)
   {
      elementaryEdgeArray[ii] = elementaryEdgeArray2[ii];
   }
   // #pragma acc enter data copyin(elementaryArray[0:eleCount])
   // #pragma acc enter data copyin(elementaryEdgeArray[0:edgeCount])

   numLeaves = edgeCount / 2;
   // #pragma acc enter data copyin(elementaryArray[0:eleCount], elementaryEdgeArray[0:edgeCount]) copied in the constructTreeIntervalArray()

   #ifdef DG 
   cp6 = high_resolution_clock::now();
   #endif
   #ifdef DG
      // Calculating total time taken by the program.
      auto d1 = duration_cast<microseconds>(cp2 - cp1);
      auto d2 = duration_cast<microseconds>(cp3 - cp2);
      auto d3 = duration_cast<microseconds>(cp4 - cp3);
      auto d4 = duration_cast<microseconds>(cp5 - cp4);
      auto d5 = duration_cast<microseconds>(cp6 - cp5);

      cout << "CEIA: copy() " << fixed << d1.count() << setprecision(10) << endl;
      cout << "CEIA: radix sort() " << fixed << d2.count() << setprecision(10) << endl;
      cout << "CEIA: dup remove() " << fixed << d3.count() << setprecision(10) << endl;
      cout << "CEIA: last elements handling() " << fixed << d4.count() << setprecision(10) << endl;
      cout << "CEIA: to GPU copy() " << fixed << d5.count() << setprecision(10) << endl;
   #endif
}

/*
building tree interval array
*/
void SegTree::constructTreeIntervalArray()
{
   int ts = treeSize, ts2, ld = 2;
   int lcount = log2(ts), ipc;

   ipc = numLeaves;
   // // ipc=elementaryPointCount;
   // ts2 = ts * 2;

   // #ifdef DG 
   // high_resolution_clock::time_point cp1, cp2, cp3, cp4;
   // cp1 = high_resolution_clock::now();
   // #endif

   // #pragma acc enter data copyin(this)
   // #pragma acc enter data copyin(elementaryEdgeArray[0 : ipc * 2])

   // #ifdef DG 
   // cp4 = high_resolution_clock::now();
   // #endif

   ts2 = ts * 2;
   #pragma acc enter data copyin(this)
   #pragma acc enter data create(nodeIntervalArray[0:ts*2])
   #pragma acc data present(elementaryEdgeArray[0:ipc*2])
   #pragma acc parallel loop
   for (int i = 0; i < ipc * 2; ++i)
   {
      nodeIntervalArray[ts + i] = elementaryEdgeArray[i];
   }

   for (int lid = lcount - 1; lid >= 0; --lid)
   {
      int lstart = ts / (ld * 2);
      int lend = lstart * 2;
      // cout<<"\n"<<lid<<" **** "<<ld<<endl;

      #pragma acc parallel loop
      for (int nnid = 0; nnid < (lend-lstart); nnid++)
      {
         int nid = nnid+lstart;
         int nid1 = nid * 2;
         int nid2 = nid1 + 1;
         // cout<<"*** nid="<<nid<<" nid1="<<nid1<<" nid2="<<nid2;
         // parent start is start of first child. End is end of second child
         nodeIntervalArray[nid * 2] = nodeIntervalArray[nid1 * 2];   //start of first child
         nodeIntervalArray[nid * 2 + 1] = nodeIntervalArray[nid2 * 2 + 1]; //End is end of second child
      }
      // #pragma acc wait

      ld *= 2;
   }
   // #ifdef DG 
   // cp2 = high_resolution_clock::now();
   // #endif

   // REAL *nodeIntervalArray2;
   // nodeIntervalArray2=(REAL *)malloc((ts2)*sizeof(REAL));
   // for(int i=0; i<ts2; ++i){
   //    nodeIntervalArray2[i]=nodeIntervalArray[i];
   // }
   // cout<<"ts2 "<<ts2<<endl;
   int eleCount = elementaryPointCount;
   // #pragma acc update self(nodeIntervalArray[0 : ts2])

   // #ifdef DG 
   // cp3 = high_resolution_clock::now();
   // #endif
   // #ifdef DG
   //    // Calculating total time taken by the program.
   //    auto d1 = duration_cast<microseconds>(cp2 - cp1);
   //    // auto d2 = duration_cast<microseconds>(cp3 - cp2);
   //    // auto d3 = duration_cast<microseconds>(cp4 - cp1);

   //    cout << "constructIntervals: construct() " << fixed << d1.count() << setprecision(10) << endl;
   //    // cout << "constructIntervals: to GPU update() " << fixed << d2.count() << setprecision(10) << endl;
   //    // cout << "constructIntervals: to GPU copy() " << fixed << d3.count() << setprecision(10) << endl;
   // #endif
}

/*
building tree interval array. Multicore
*/
void SegTree::constructTreeIntervalArrayMulticore()
{
   int ts = treeSize, ts2, ld = 2;
   int lcount = log2(ts), ipc;

   ipc = numLeaves;
   // ipc=elementaryPointCount;
   ts2 = ts * 2;
   #ifdef DG 
   high_resolution_clock::time_point cp1, cp2, cp3, cp4;
   cp1 = high_resolution_clock::now();
   #endif

   // #pragma omp parallel for num_threads(20)
   for (int i = 0; i < ipc * 2; ++i)
   {
      nodeIntervalArray[ts + i] = elementaryEdgeArray[i];
   }

   for (int lid = lcount - 1; lid >= 0; --lid)
   {
      int lstart = ts / (ld * 2);
      int lend = lstart * 2;

      #pragma omp parallel for
      for (int nid = lstart; nid < lend; nid++)
      {
         int nid1 = nid * 2;
         int nid2 = nid1 + 1;
         // parent start is start of first child. End is end of second child
         nodeIntervalArray[nid * 2] = nodeIntervalArray[nid1 * 2];
         nodeIntervalArray[nid * 2 + 1] = nodeIntervalArray[nid2 * 2 + 1];
      }
      #pragma omp barrier

      ld *= 2;
   }
   #ifdef DG 
   cp2 = high_resolution_clock::now();
   #endif
   int eleCount = elementaryPointCount;
   #ifdef DG
      // Calculating total time taken by the program.
      auto d1 = duration_cast<microseconds>(cp2 - cp1);
      cout << "constructIntervals: construct() " << fixed << d1.count() << setprecision(10) << endl;
   #endif
}

vector<Edge> *SegTree ::getElementaryIntervals(Edge *intervals)
{
   vector<REAL> *pts = getPointsFromEdges2(intervals);
   m_elementaryPoints = pts;

   vector<Edge> *elementaryIntervals = new vector<Edge>();

   for (int i = 0; i < pts->size() - 1; i++)
   {
      REAL start = pts->at(i);
      REAL end = pts->at(i + 1);
      Edge e(start, end);
      elementaryIntervals->push_back(e);
   }

   // find upper bound on #intervals that is a 2 exponent
   int numIntervals = elementaryIntervals->size();
   int upperBound = Util::getNPowOf2(numIntervals);

   // normalize : Make the number of edges power of 2 for convenience
   int i = pts->size() - 1;
   for (int j = i; j < upperBound; j++)
   {
      REAL start = pts->at(i);
      REAL end = pts->at(i);
      Edge e(start, end);
      elementaryIntervals->push_back(e);

      // this part is useful for endlist for leaves. We want to normalize the points 2^integer
      m_elementaryPoints->push_back(start);
   }

   return elementaryIntervals;
}

/*
count end list sizes at each node
does the prefix sum as well
*/

void SegTree::countEndlistCoverListSizes()
{
   int ts = treeSize;
   int leafSize = 0;

   coverListSizes[0] = 0;
   coverListSizes[1] = 0;
   // cout<<"countEndlistCoverListSizes: start"<<endl;

   // if(treeNode->at(nodeIdx)==nullptr)

   for (int nid = 2; nid <= treeSize; ++nid)
   {
      // cout<<"nid="<<nid<<" "<<endl;
      // cout<<"nid="<<nid<<" "<<(treeNode->at(nid-1)==nullptr)<<endl;
      if (treeNode->at(nid-1)!=nullptr){
         // cout<<"nid="<<nid<<" "<<endl;
         // cout<<"nid="<<nid<<" "<<treeNode->at(nid-1)->coverList->size()<<endl;
         coverListSizes[nid]=coverListSizes[nid-1]+treeNode->at(nid-1)->coverList->size();
      }else{
         coverListSizes[nid]=coverListSizes[nid-1];
      }
   }

   sumCoverListSize = coverListSizes[treeSize];
   // cout<<"countEndlistCoverListSizes: counted"<<endl;
   #pragma acc enter data copyin(coverListSizes[0 : (ts + 1)])
}

/*
multicore
*/
void SegTree::countEndlistCoverListSizesMulticore()
{
   int ts = treeSize;
   int leafSize = 0;

   coverListSizes[0] = 0;
   coverListSizes[1] = 0;

   for (int nid = 2; nid <= treeSize; ++nid)
   {
      coverListSizes[nid] = coverListSizes[nid - 1] + treeNode->at(nid - 1)->coverList->size();
   }

   sumCoverListSize = coverListSizes[treeSize];
}

/*
Converts cover lists of the segment tree into a 2-D array which has nodes and their coverlists copy it to GPU
*/

void SegTree ::copyNodeCoverListToGPU(REAL *bPoly, REAL *cPoly, int bSize, int cPolysSize)
{
   coverArray = new REAL[sumCoverListSize * 4];
   // vector<Edge> *coverList;
   int scl = sumCoverListSize * 4;
   int sels = sumEndListSize;

   int bpolys = (bSize + 1) * 2, cpolys = (cPolysSize) * 2;

   int ts = treeSize;
   int bitArraySize = (ts + 8) / 8;
   // #pragma acc update self(coverListSizes2[0:bitArraySize])

   int nstart, clid2;
   for (int nid = 1, nid2 = 0; nid < treeSize; ++nid, nid2 += 2)
   {
      if(treeNode->at(nid)==nullptr)
         continue;
      // coverList=treeNode->at(nid)->coverList;
      auto cle = treeNode->at(nid)->coverList->begin();
      nstart = coverListSizes[nid] * 4;
      clid2 = nstart;

      for (int clid = 0, clid2 = nstart; clid < (coverListSizes[nid+1]-coverListSizes[nid]); ++clid, clid2 += 4)
      {
         // coverArray[clid2]=coverList->at(clid).id;
         // coverArray[clid2+1]=coverList->at(clid).type;
         // coverArray[clid2+2]=coverList->at(clid).start;
         // coverArray[clid2+3]=coverList->at(clid).end;
         coverArray[clid2] = cle->id;
         coverArray[clid2 + 1] = cle->type;
         coverArray[clid2 + 2] = cle->start;
         coverArray[clid2 + 3] = cle->end;
         cle++;
      }
   }
   // #pragma acc enter data copyin(this, coverArray[0:sumCoverListSize*4])
   // #pragma acc enter data copyin(coverArray[0 : scl])
   // #pragma acc update device(treeNodeEndListArray[0 : sels], endListSizes[0 : ts], endListCounts[0 : ts])
   // #pragma acc enter data copyin(coverArray[0 : scl], bPoly[0 : bpolys], cPoly[0 : cpolys])
   #pragma acc enter data copyin(coverArray[0 : scl])
   

}

/*
Converts cover lists of the segment tree into a 2-D array which has nodes and their coverlists copy it to GPU
Multicore
*/
void SegTree ::copyNodeCoverListMulticore()
{
   coverArray = new REAL[sumCoverListSize * 4];
   // vector<Edge> *coverList;
   int scl = sumCoverListSize * 4;

   int ts = treeSize;

   int nstart, clid2;
   // #pragma omp parallel for private(clid2, nstart) schedule(dynamic, 100)
   // for(int nid=1, nid2=0; nid<treeSize; ++nid, nid2+=2)
   for (int nid = 1; nid < treeSize; ++nid)
   {
      // coverList=iclid;
      // coverList=treeNode->at(nid)->coverList;
      nstart = coverListSizes[nid] * 4;
      clid2 = nstart;

      for (auto it = treeNode->at(nid)->coverList->begin(); it != treeNode->at(nid)->coverList->end(); ++it)
      // for (int clid=0, clid2=nstart; clid < coverList->size(); ++clid, clid2 += 4)
      {
         // #pragma omp critical
         {
            coverArray[clid2] = it->id;
            coverArray[clid2 + 1] = it->type;
            coverArray[clid2 + 2] = it->start;
            coverArray[clid2 + 3] = it->end;
            clid2 += 4;
         }
      }
   }
}

/*
Converts endEdgeSet list of the segment tree into a 2-D array which has nodes and their endlists copy it to GPU
*/
void SegTree ::copyNodeEndEdgeSetToGPU()
{
   // endArray=new REAL[sumEndListSize*4];
   cout << "sumEndListSize " << sumEndListSize << endl;
   endArray = (REAL *)malloc((sumEndListSize * 4) * sizeof(REAL));

   // set<Edge, edge_compare> *endEdgeSet;

   int elid2, nstart, end;
   // cout<<"****** end Array *******"<<endl;
   for (int nid = 1; nid < treeSize; ++nid)
   {
      // endEdgeSet=treeNode->at(nid)->endEdgeSet;
      nstart = endListSizes[nid] * 4;
      elid2 = nstart;

      end = (endListSizes[nid] + endListCounts[nid]) * 4;
      for (int i = endListSizes[nid]; i < (endListSizes[nid] + endListCounts[nid]); ++i)
      {
         *(endArray + (elid2)) = treeNodeEndListArray[i].id;
         *(endArray + (elid2 + 1)) = treeNodeEndListArray[i].type;
         *(endArray + (elid2 + 2)) = treeNodeEndListArray[i].start;
         *(endArray + (elid2 + 3)) = treeNodeEndListArray[i].end;
         // cout<<"nid="<<nid<<" e.id: "<<*(endArray + (nid * maxRowSize + elid2))<<endl;
         elid2 += 4;
      }
      if (elid2 > end)
      {
         cout << "ERROR: copyNodeEndEdgeSetToGPU(). End list overflow" << endl;
      }
   }
   // #pragma acc enter data copyin(this, endArray[0:sumEndListSize*4])
}

/*
Query using the 2-D array
*/
/*void SegTree ::recQueryArrayParallel(REAL qx)
{
   int ts = treeSize, scl = sumCoverListSize;
   int maxCoverListSize = sumCoverListSize / (ts * 4);
   int maxRowSize = maxCoverListSize * 4;

   int *resultSizes = new int[ts];
   int nid, clid, clid2;

   int *psCLSizes = new int[ts + 1]; // inclusive prefix sum

   #pragma acc enter data create(resultSizes [0:ts]) copyin(this)
   #pragma acc data present(coverListSizes [0:ts], nodeIntervalArray [0:ts * 2])

   #pragma acc parallel loop
   for (int nid = 0; nid < ts; ++nid)
   {
      resultSizes[nid] = 0;
      if (qx >= nodeIntervalArray[nid * 2] && qx <= nodeIntervalArray[nid * 2 + 1])
      {
         resultSizes[nid] = coverListSizes[nid];
      }
   }

   int nid2 = 0;
   // prefix sum of the resultsizes ************* need parallel version
   // ------------------------------------------------------
   #pragma acc update self(resultSizes [0:ts], nodeIntervalArray [0:ts * 2])
   psCLSizes[0] = 0;
   for (nid = 1, nid2 = 2; nid <= ts; ++nid, nid2 += 2)
   {
      psCLSizes[nid] = psCLSizes[nid - 1] + resultSizes[nid - 1];
      // cout<<nid<<">"<<resultSizes[nid]<<", >> ("<<nodeIntervalArray[nid2]<<" "<<nodeIntervalArray[nid2+1]<<") ";
   }
   // cout<<endl;
   // ------------------------------------------------------

   int rs = psCLSizes[ts];
   REAL *results = new REAL[rs];
   int pos;

   // #pragma acc enter data copyin(this, psCLSizes[0:ts+1])
   #pragma acc enter data create(results [0:rs]) copyin(psCLSizes [0:ts + 1])
   #pragma acc data present(coverListSizes [0:ts], nodeIntervalArray [0:ts * 2], coverArray [0:scl])

   #pragma acc parallel loop seq
   for (nid = 0; nid < ts; ++nid)
   {
      if (qx >= nodeIntervalArray[nid * 2] && qx <= nodeIntervalArray[nid * 2 + 1])
      {
         pos = psCLSizes[nid];
         for (clid = 0, clid2 = 0; clid < coverListSizes[nid]; ++clid, clid2 += 4)
         {
            results[pos + clid] = coverArray[nid * maxRowSize + clid2];
         }
      }
   }

   // testing the results
   #pragma acc update self(results [0:rs])
   cout << "Result CL size: " << rs << endl;
   for (int i = 0; i < rs; ++i)
   {
      cout << i << ">" << results[i] << ", ";
   }
   cout << endl;
}*/

/*
LSMF filter
*/
// #pragma acc routine seq
extern "C" __device__ bool LSMF(REAL l1p1x, REAL l1p1y, REAL l1p2x, REAL l1p2y,
                                REAL l2p1x, REAL l2p1y, REAL l2p2x, REAL l2p2y)
{
   REAL minPX = l1p1x, minPY = l1p1y;
   REAL maxPX = l1p2x, maxPY = l1p2y;
   REAL minQX = l2p1x, minQY = l2p1y;
   REAL maxQX = l2p2x, maxQY = l2p2y;
   // this toggle way optimizes this computation well compared to using 8 min max calls seperately
   if (minPX > l1p2x)
   {
      minPX = l1p2x;
      maxPX = l1p1x;
   }
   if (minPY > l1p2y)
   {
      minPY = l1p2y;
      maxPY = l1p1y;
   }
   if (minQX > l2p2x)
   {
      minQX = l2p2x;
      maxQX = l2p1x;
   }
   if (minQY > l2p2y)
   {
      minQY = l2p2y;
      maxQY = l2p1y;
   }
   // check intersection between MBRs
   if (minPX > maxQX || maxPX < minQX)
      return false;
   if (minPY > maxQY || maxPY < minQY)
      return false;
   return true;
}

extern "C" __device__ double A_(const point &P, const point &Q, const point &R)
// double A(const point &P, const point &Q, const point &R)
{
   return (Q.x - P.x) * (R.y - P.y) - (Q.y - P.y) * (R.x - P.x);
}
// multiply two 2D points
extern "C" __device__ double mul_(const point &a, const point &b)
// double mul(const point &a, const point &b)
{
   point r;
   r.x = a.x * b.x;
   r.y = a.y * b.y;
   return (r.x + r.y);
}

// extern "C" __device__ double A(const point &P, const point &Q, const point &R)
double A(const point &P, const point &Q, const point &R)
{
   return (Q.x - P.x) * (R.y - P.y) - (Q.y - P.y) * (R.x - P.x);
}
// multiply two 2D points
// extern "C" __device__ double mul(const point &a, const point &b)
double mul(const point &a, const point &b)
{
   point r;
   r.x = a.x * b.x;
   r.y = a.y * b.y;
   return (r.x + r.y);
}

// multiply scalar with 2D points
point mulScalar(const double c, const point &b)
{
   point r;
   r.x = c * b.x;
   r.y = c * b.y;
   return r;
}
// add two 2D points
point add(const point &a, const point &b)
{
   point r;
   r.x = a.x + b.x;
   r.y = a.y + b.y;
   return r;
}
// difference of two 2D points
extern "C" __device__ point sub_(const point &a, const point &b)
// point sub(const point &a, const point &b)
{
   point r;
   r.x = a.x - b.x;
   r.y = a.y - b.y;
   return r;
}

// double device_fabs(double val1)
extern "C" __device__ double device_fabs_(double val1)
{
   if (val1 >= 0)
      return val1;
   else
      return val1 * (-1);
}

#pragma acc routine seq
extern "C" double A_(const point &P, const point &Q, const point &R);
#pragma acc routine seq
extern "C" point sub_(const point &a, const point &b);
#pragma acc routine seq
extern "C" double mul_(const point &a, const point &b);
#pragma acc routine seq
extern "C" double device_fabs_(double val1);

point sub(const point &a, const point &b)
{
   point r;
   r.x = a.x - b.x;
   r.y = a.y - b.y;
   return r;
}

double device_fabs(double val1)
{
   if (val1 >= 0)
      return val1;
   else
      return val1 * (-1);
}

/*
intersect two given line segments
  NO_INTERSECTION, //0
  X_INTERSECTION,  //1
  T_INTERSECTION_Q, //2
  T_INTERSECTION_P, //3
  V_INTERSECTION, //4
  X_OVERLAP,      //5
  T_OVERLAP_Q,    //6
  T_OVERLAP_P,    //7
  V_OVERLAP       //8
*/
// #pragma acc routine seq
// extern "C" __device__ int lineSegmentIntersection
//       (REAL l1p1x, REAL l1p1y, REAL l1p2x, REAL l1p2y,
//       REAL l2p1x, REAL l2p1y, REAL l2p2x, REAL l2p2y,
//       double &alpha, double &beta){
// extern "C" __device__ int lineSegmentIntersection
int lineSegmentIntersection(REAL l1p1x, REAL l1p1y, REAL l1p2x, REAL l1p2y,
                            REAL l2p1x, REAL l2p1y, REAL l2p2x, REAL l2p2y,
                            double *alpha, double *beta)
{

   point l1p1, l1p2, l2p1, l2p2;
   l1p1.x = l1p1x;
   l1p1.y = l1p1y;
   l1p2.x = l1p2x;
   l1p2.y = l1p2y;
   l2p1.x = l2p1x;
   l2p1.y = l2p1y;
   l2p2.x = l2p2x;
   l2p2.y = l2p2y;
   double AP1 = A(l1p1, l2p1, l2p2);
   double AP2 = A(l1p2, l2p1, l2p2);

   if (fabs(AP1 - AP2) > EPSILON)
   {
      // from here: [P1,P2] and [Q1,Q2] are not parallel
      // analyse potential intersection
      double AQ1 = A(l2p1, l1p1, l1p2);
      double AQ2 = A(l2p2, l1p1, l1p2);
      // compute alpha and beta
      *alpha = AP1 / (AP1 - AP2);
      *beta = AQ1 / (AQ1 - AQ2);
      // classify alpha
      bool alpha_is_0 = false;
      bool alpha_in_0_1 = false;
      if ((*alpha > EPSILON) && (*alpha < 1.0 - EPSILON))
         alpha_in_0_1 = true;
      else if (fabs(*alpha) <= EPSILON)
         alpha_is_0 = true;
      // classify beta
      bool beta_is_0 = false;
      bool beta_in_0_1 = false;
      if ((*beta > EPSILON) && (*beta < 1.0 - EPSILON))
         beta_in_0_1 = true;
      else if (fabs(*beta) <= EPSILON)
         beta_is_0 = true;
      // distinguish intersection types
      if (alpha_in_0_1 && beta_in_0_1)
         return (1); // return (X_INTERSECTION);
      if (alpha_is_0 && beta_in_0_1)
         return (2); // return (T_INTERSECTION_Q);
      if (beta_is_0 && alpha_in_0_1)
         return (3); // return (T_INTERSECTION_P);
      if (alpha_is_0 && beta_is_0)
         return (4); // return (V_INTERSECTION);
   }
   else if (fabs(AP1) < EPSILON)
   {
      // from here: [P1,P2] and [Q1,Q2] are collinear
      // analyse potential overlap
      point dP = sub(l1p2, l1p1);
      point dQ = sub(l2p2, l2p1);
      point PQ = sub(l2p1, l1p1);
      *alpha = mul(PQ, dP) / mul(dP, dP);
      *beta = -mul(PQ, dQ) / mul(dQ, dQ);
      // classify alpha
      bool alpha_is_0 = false;
      bool alpha_in_0_1 = false;
      bool alpha_not_in_0_1 = false;
      if ((*alpha > EPSILON) && (*alpha < 1.0 - EPSILON))
         alpha_in_0_1 = true;
      else if (fabs(*alpha) <= EPSILON)
         alpha_is_0 = true;
      else
         alpha_not_in_0_1 = true;
      // classify beta
      bool beta_is_0 = false;
      bool beta_in_0_1 = false;
      bool beta_not_in_0_1 = false;
      if ((*beta > EPSILON) && (*beta < 1.0 - EPSILON))
         beta_in_0_1 = true;
      else if (fabs(*alpha) <= EPSILON)
         beta_is_0 = true;
      else
         beta_not_in_0_1 = true;

      // distinguish intersection types
      if (alpha_in_0_1 && beta_in_0_1)
         return (5); // return (X_OVERLAP);
      if (alpha_not_in_0_1 && beta_in_0_1)
         return (6); // return (T_OVERLAP_Q);
      if (beta_not_in_0_1 && alpha_in_0_1)
         return (7); // return (T_OVERLAP_P);
      if (alpha_is_0 && beta_is_0)
         return (8); // return (V_OVERLAP);
   }
   return (0); // return (NO_INTERSECTION);
}

// extern "C"  int lineSegmentIntersection2
extern "C" __device__ int lineSegmentIntersection2(REAL l1p1x, REAL l1p1y, REAL l1p2x, REAL l1p2y,
                                                   REAL l2p1x, REAL l2p1y, REAL l2p2x, REAL l2p2y, double epsilon)
{

   double alpha, beta;
   point l1p1, l1p2, l2p1, l2p2;
   l1p1.x = l1p1x;
   l1p1.y = l1p1y;
   l1p2.x = l1p2x;
   l1p2.y = l1p2y;
   l2p1.x = l2p1x;
   l2p1.y = l2p1y;
   l2p2.x = l2p2x;
   l2p2.y = l2p2y;
   double AP1 = A(l1p1, l2p1, l2p2);
   double AP2 = A(l1p2, l2p1, l2p2);

   if (device_fabs(AP1 - AP2) > epsilon)
   {
      // from here: [P1,P2] and [Q1,Q2] are not parallel
      // analyse potential intersection
      double AQ1 = A(l2p1, l1p1, l1p2);
      double AQ2 = A(l2p2, l1p1, l1p2);
      // compute alpha and beta
      alpha = AP1 / (AP1 - AP2);
      beta = AQ1 / (AQ1 - AQ2);
      // classify alpha
      bool alpha_is_0 = false;
      bool alpha_in_0_1 = false;
      if ((alpha > epsilon) && (alpha < 1.0 - epsilon))
         alpha_in_0_1 = true;
      else if (device_fabs(alpha) <= epsilon)
         alpha_is_0 = true;
      // classify beta
      bool beta_is_0 = false;
      bool beta_in_0_1 = false;
      if ((beta > epsilon) && (beta < 1.0 - epsilon))
         beta_in_0_1 = true;
      else if (device_fabs(beta) <= epsilon)
         beta_is_0 = true;
      // distinguish intersection types
      if (alpha_in_0_1 && beta_in_0_1)
         return (1); // return (X_INTERSECTION);
      if (alpha_is_0 && beta_in_0_1)
         return (2); // return (T_INTERSECTION_Q);
      if (beta_is_0 && alpha_in_0_1)
         return (3); // return (T_INTERSECTION_P);
      if (alpha_is_0 && beta_is_0)
         return (4); // return (V_INTERSECTION);
   }
   else if (device_fabs(AP1) < epsilon)
   {
      // from here: [P1,P2] and [Q1,Q2] are collinear
      // analyse potential overlap
      point dP = sub(l1p2, l1p1);
      point dQ = sub(l2p2, l2p1);
      point PQ = sub(l2p1, l1p1);
      alpha = mul(PQ, dP) / mul(dP, dP);
      beta = -mul(PQ, dQ) / mul(dQ, dQ);
      // classify alpha
      bool alpha_is_0 = false;
      bool alpha_in_0_1 = false;
      bool alpha_not_in_0_1 = false;
      if ((alpha > epsilon) && (alpha < 1.0 - epsilon))
         alpha_in_0_1 = true;
      else if (device_fabs(alpha) <= epsilon)
         alpha_is_0 = true;
      else
         alpha_not_in_0_1 = true;
      // classify beta
      bool beta_is_0 = false;
      bool beta_in_0_1 = false;
      bool beta_not_in_0_1 = false;
      if ((beta > epsilon) && (beta < 1.0 - epsilon))
         beta_in_0_1 = true;
      else if (device_fabs(alpha) <= epsilon)
         beta_is_0 = true;
      else
         beta_not_in_0_1 = true;

      // distinguish intersection types
      if (alpha_in_0_1 && beta_in_0_1)
         return (5); // return (X_OVERLAP);
      if (alpha_not_in_0_1 && beta_in_0_1)
         return (6); // return (T_OVERLAP_Q);
      if (beta_not_in_0_1 && alpha_in_0_1)
         return (7); // return (T_OVERLAP_P);
      if (alpha_is_0 && beta_is_0)
         return (8); // return (V_OVERLAP);
   }
   return (0); // return (NO_INTERSECTION);
}

extern "C" __device__ int lineSegmentIntersectionDev(REAL l1p1x, REAL l1p1y, REAL l1p2x, REAL l1p2y,
                                                     REAL l2p1x, REAL l2p1y, REAL l2p2x, REAL l2p2y,
                                                     double *alpha, double *beta)
{

   point l1p1, l1p2, l2p1, l2p2;
   l1p1.x = l1p1x;
   l1p1.y = l1p1y;
   l1p2.x = l1p2x;
   l1p2.y = l1p2y;
   l2p1.x = l2p1x;
   l2p1.y = l2p1y;
   l2p2.x = l2p2x;
   l2p2.y = l2p2y;
   double AP1 = A(l1p1, l2p1, l2p2);
   double AP2 = A(l1p2, l2p1, l2p2);

   if (fabs(AP1 - AP2) > EPSILON)
   {
      // from here: [P1,P2] and [Q1,Q2] are not parallel
      // analyse potential intersection
      double AQ1 = A(l2p1, l1p1, l1p2);
      double AQ2 = A(l2p2, l1p1, l1p2);
      // compute alpha and beta
      *alpha = AP1 / (AP1 - AP2);
      *beta = AQ1 / (AQ1 - AQ2);
      // classify alpha
      bool alpha_is_0 = false;
      bool alpha_in_0_1 = false;
      if ((*alpha > EPSILON) && (*alpha < 1.0 - EPSILON))
         alpha_in_0_1 = true;
      else if (fabs(*alpha) <= EPSILON)
         alpha_is_0 = true;
      // classify beta
      bool beta_is_0 = false;
      bool beta_in_0_1 = false;
      if ((*beta > EPSILON) && (*beta < 1.0 - EPSILON))
         beta_in_0_1 = true;
      else if (fabs(*beta) <= EPSILON)
         beta_is_0 = true;
      // distinguish intersection types
      if (alpha_in_0_1 && beta_in_0_1)
         return (1); // return (X_INTERSECTION);
      if (alpha_is_0 && beta_in_0_1)
         return (2); // return (T_INTERSECTION_Q);
      if (beta_is_0 && alpha_in_0_1)
         return (3); // return (T_INTERSECTION_P);
      if (alpha_is_0 && beta_is_0)
         return (4); // return (V_INTERSECTION);
   }
   else if (fabs(AP1) < EPSILON)
   {
      // from here: [P1,P2] and [Q1,Q2] are collinear
      // analyse potential overlap
      point dP = sub(l1p2, l1p1);
      point dQ = sub(l2p2, l2p1);
      point PQ = sub(l2p1, l1p1);
      *alpha = mul(PQ, dP) / mul(dP, dP);
      *beta = -mul(PQ, dQ) / mul(dQ, dQ);
      // classify alpha
      bool alpha_is_0 = false;
      bool alpha_in_0_1 = false;
      bool alpha_not_in_0_1 = false;
      if ((*alpha > EPSILON) && (*alpha < 1.0 - EPSILON))
         alpha_in_0_1 = true;
      else if (fabs(*alpha) <= EPSILON)
         alpha_is_0 = true;
      else
         alpha_not_in_0_1 = true;
      // classify beta
      bool beta_is_0 = false;
      bool beta_in_0_1 = false;
      bool beta_not_in_0_1 = false;
      if ((*beta > EPSILON) && (*beta < 1.0 - EPSILON))
         beta_in_0_1 = true;
      else if (fabs(*alpha) <= EPSILON)
         beta_is_0 = true;
      else
         beta_not_in_0_1 = true;

      // distinguish intersection types
      if (alpha_in_0_1 && beta_in_0_1)
         return (5); // return (X_OVERLAP);
      if (alpha_not_in_0_1 && beta_in_0_1)
         return (6); // return (T_OVERLAP_Q);
      if (beta_not_in_0_1 && alpha_in_0_1)
         return (7); // return (T_OVERLAP_P);
      if (alpha_is_0 && beta_is_0)
         return (8); // return (V_OVERLAP);
   }
   return (0); // return (NO_INTERSECTION);
}

/*
coverArray: 4*size
   [id, type, start, end]
*/

#pragma acc routine seq
extern "C" bool LSMF(REAL l1p1x, REAL l1p1y, REAL l1p2x, REAL l1p2y,
                     REAL l2p1x, REAL l2p1y, REAL l2p2x, REAL l2p2y);
#pragma acc routine seq
extern "C" int lineSegmentIntersection2(REAL l1p1x, REAL l1p1y, REAL l1p2x, REAL l1p2y,
                                        REAL l2p1x, REAL l2p1y, REAL l2p2x, REAL l2p2y, double epsilon);
// #pragma acc routine seq
// extern "C" int lineSegmentIntersection
//       (REAL l1p1x, REAL l1p1y, REAL l1p2x, REAL l1p2y,
//       REAL l2p1x, REAL l2p1y, REAL l2p2x, REAL l2p2y,
//       double *alpha, double *beta);
extern "C" int lineSegmentIntersectionDev(REAL l1p1x, REAL l1p1y, REAL l1p2x, REAL l1p2y,
                                          REAL l2p1x, REAL l2p1y, REAL l2p2x, REAL l2p2y,
                                          double *alpha, double *beta);

/*
count pairs at each node
*/
void SegTree ::countLineSegIntersectsAtNodes(int *iCount)
{

   int ts = treeSize, scl = sumCoverListSize * 4;
   int sels = sumEndListSize;

   for (int nid = 1; nid < ts; nid++)
   {
      int interCount = 0;
      int cs = coverListSizes[nid], ce = coverListSizes[nid + 1];
      for (int clid = cs; clid < ce; ++clid)
      {
         int clid2 = clid * 4;
         // b \in C_B(v) and o \in C_O(v)
         for (int clid3 = clid + 1; clid3 < ce; ++clid3)
         {
            int bId, cId;
            int clid4 = clid3 * 4;
            // check if line segments are from different polygons
            if (coverArray[clid2 + 1] != coverArray[clid4 + 1])
            {
               interCount += 1;
            }
         }
         //    // b \in C_B(v) and o \in E_O(v)
         int es = endListSizes[nid], ee = endListSizes[nid] + endListCounts[nid];
         for (int elid = es; elid < ee; ++elid)
         {
            int bId, cId;
            int elid2 = elid * 4;
            // check if line segments are from different polygons
            // if(coverArray[clid2+1]!=endArray[elid2+1]){
            if (coverArray[clid2 + 1] != treeNodeEndListArray[elid].type)
            {
               interCount += 1;
            }
         }
      }
      iCount[nid - 1] = interCount;
   }
}

/*
Count intersections
*/
void SegTree ::saveIntersectionsIds(REAL *bPoly, REAL *cPoly, int bSize, int cSize,
                                    int bType, int cType,
                                    vector<int> &bPolLineIds, vector<int> &cPolLineIds, vector<int> &intersectTypesDup)
{

   int ts = treeSize, scl = sumCoverListSize * 4;
   int sels = sumEndListSize;
   int bpolys = (bSize + 1) * 2, cpolys = (cSize + 1) * 2;
   double epsilon = EPSILON;
   bool bolpr;

   for (int nid = 1; nid < ts; nid++)
   {
      int cs = coverListSizes[nid], ce = coverListSizes[nid + 1];
      for (auto clid = treeNode->at(nid)->coverList->begin(); clid != treeNode->at(nid)->coverList->end(); ++clid)
      {
         // b \in C_B(v) and o \in C_O(v)
         // auto nclid=clid;
         // nclid++;
         for (auto clid3 = treeNode->at(nid)->coverList->begin(); clid3 != treeNode->at(nid)->coverList->end(); ++clid3)
         {
            // for(auto clid3=nclid; clid3!=treeNode->at(nid)->coverList->end(); ++clid3){
            int bId, cId;
            bolpr = false;
            // check if line segments are from different polygons
            // if (clid->type!=clid3->type){
            // If clid2 is base polygon. Else clid2 is subject polygon
            if (clid->type == bType && clid3->type == cType)
            {
               bId = clid->id;
               cId = clid3->id;
               bolpr = true;
            }
            else if (clid->type == cType && clid3->type == bType)
            {
               bId = clid3->id;
               cId = clid->id;
               bolpr = true;
            }
            if (bolpr)
            {
               if (LSMF(bPoly[bId], bPoly[bId + 1],
                        bPoly[bId + 2], bPoly[bId + 3],
                        cPoly[cId], cPoly[cId + 1],
                        cPoly[cId + 2], cPoly[cId + 3]))
               {
                  int intersecType = lineSegmentIntersection2(bPoly[bId], bPoly[bId + 1],
                                                              bPoly[bId + 2], bPoly[bId + 3],
                                                              cPoly[cId], cPoly[cId + 1],
                                                              cPoly[cId + 2], cPoly[cId + 3], epsilon);
                  if (intersecType)
                  {
                     bPolLineIds.push_back(bId);
                     cPolLineIds.push_back(cId);
                     intersectTypesDup.push_back(intersecType);
                  }
               }
            }
         }
         //    // b \in C_B(v) and o \in E_O(v)
         for (auto elid = treeNode->at(nid)->endEdgeSet->begin(); elid != treeNode->at(nid)->endEdgeSet->end(); ++elid)
         {
            int bId, cId;
            bolpr = false;
            // check if line segments are from different polygons
            // if(clid->type!=elid->type){
            // If clid2 is base polygon. Else clid2 is subject polygon
            if (clid->type == bType && elid->type == cType)
            {
               bId = clid->id;
               cId = elid->id;
               bolpr = true;
            }
            else if (clid->type == cType && elid->type == bType)
            {
               cId = clid->id;
               bId = elid->id;
               bolpr = true;
            }
            if (bolpr)
            {
               if (LSMF(bPoly[bId], bPoly[bId + 1],
                        bPoly[bId + 2], bPoly[bId + 3],
                        cPoly[cId], cPoly[cId + 1],
                        cPoly[cId + 2], cPoly[cId + 3]))
               {
                  int intersecType = lineSegmentIntersection2(bPoly[bId], bPoly[bId + 1],
                                                              bPoly[bId + 2], bPoly[bId + 3],
                                                              cPoly[cId], cPoly[cId + 1],
                                                              cPoly[cId + 2], cPoly[cId + 3], epsilon);
                  if (intersecType)
                  {
                     bPolLineIds.push_back(bId);
                     cPolLineIds.push_back(cId);
                     intersectTypesDup.push_back(intersecType);
                  }
               }
            }
         }
      }
   }
}

/*
save intersections multicore
*/
void SegTree ::saveIntersectionsIdsMulticore2(REAL *bPoly, REAL *cPoly, int bSize, int cSize,
                                              int bType, int cType,
                                              vector<int> &bPolLineIds, vector<int> &cPolLineIds, vector<int> &intersectTypesDup)
{

   int ts = treeSize, scl = sumCoverListSize * 4;
   // int nid;
   int sels = sumEndListSize;
   int bpolys = (bSize + 1) * 2, cpolys = (cSize + 1) * 2;
   double epsilon = EPSILON;

   // cout<<"Total copy items "<<bpolys+cpolys<<endl;

#pragma omp parallel for schedule(dynamic, 100)
   for (int nid = 1; nid < ts; nid++)
   {
      double alpha, beta;
      // int clid, clid3, elid;
      int interCount = 0;
      int cs = coverListSizes[nid], ce = coverListSizes[nid + 1];
      for (int clid = cs; clid < ce; ++clid)
      {
         int clid2 = clid * 4;
         // #pragma omp parallel for
         for (int clid3 = clid + 1; clid3 < ce; ++clid3)
         {
            int bId, cId;
            int clid4 = clid3 * 4;
            bool proceed = false;
            // check if line segments are from different polygons
            // if (coverArray[clid2+1]!=coverArray[clid4+1]){
            // If clid2 is base polygon. Else clid2 is subject polygon
            if (coverArray[clid2 + 1] == bType && coverArray[clid4 + 1] == cType)
            {
               bId = coverArray[clid2];
               cId = coverArray[clid4];
               proceed = true;
            }
            else if (coverArray[clid2 + 1] == cType && coverArray[clid4 + 1] == bType)
            {
               bId = coverArray[clid4];
               cId = coverArray[clid2];
               proceed = true;
            }
            if (proceed)
            {
               if (LSMF(bPoly[bId], bPoly[bId + 1],
                        bPoly[bId + 2], bPoly[bId + 3],
                        cPoly[cId], cPoly[cId + 1],
                        cPoly[cId + 2], cPoly[cId + 3]))
               {
                  int intersecType = lineSegmentIntersection2(bPoly[bId], bPoly[bId + 1],
                                                              bPoly[bId + 2], bPoly[bId + 3],
                                                              cPoly[cId], cPoly[cId + 1],
                                                              cPoly[cId + 2], cPoly[cId + 3], epsilon);
                  // cout<<"nid="<<nid<<" case1 bid: "<<bId/2<<" cid: "<<cId/2<<" type: "<<intersecType<<endl;
                  if (intersecType)
                  {
#pragma omp critical
                     // #pragma omp task
                     {
                        // try{
                        bPolLineIds.push_back(bId);
                        cPolLineIds.push_back(cId);
                        intersectTypesDup.push_back(intersecType);
                        // }catch(const std::bad_alloc &) {
                        //    // ec=1;
                        // }
                     }
                  }
               }
            }
         }

         int es = endListSizes[nid], ee = endListSizes[nid] + endListCounts[nid];
         // #pragma omp nowait
         // #pragma omp parallel for
         for (int elid = es; elid < ee; ++elid)
         {
            int bId, cId;
            int elid2 = elid * 4;
            bool proceed = false;
            // check if line segments are from different polygons
            if (coverArray[clid2 + 1] == bType && treeNodeEndListArray[elid].type == cType)
            {
               bId = coverArray[clid2];
               // cId=endArray[elid2];
               cId = treeNodeEndListArray[elid].id;
               proceed = true;
            }
            else if (coverArray[clid2 + 1] == cType && treeNodeEndListArray[elid].type == bType)
            {
               cId = coverArray[clid2];
               // bId=endArray[elid2];
               bId = treeNodeEndListArray[elid].id;
               proceed = true;
            }
            if (proceed)
            {
               // cout<<nid<<" case2 bid: "<<bId/2<<" cid: "<<cId/2<<endl;
               if (LSMF(bPoly[bId], bPoly[bId + 1],
                        bPoly[bId + 2], bPoly[bId + 3],
                        cPoly[cId], cPoly[cId + 1],
                        cPoly[cId + 2], cPoly[cId + 3]))
               {
                  int intersecType = lineSegmentIntersection2(bPoly[bId], bPoly[bId + 1],
                                                              bPoly[bId + 2], bPoly[bId + 3],
                                                              cPoly[cId], cPoly[cId + 1],
                                                              cPoly[cId + 2], cPoly[cId + 3], epsilon);
                  if (intersecType)
                  {
#pragma omp critical
                     // #pragma omp task
                     {
                        // try{
                        bPolLineIds.push_back(bId);
                        cPolLineIds.push_back(cId);
                        intersectTypesDup.push_back(intersecType);
                        //  }catch(const std::bad_alloc &) {
                        //    // ec=1;
                        // }
                     }
                  }
               }
            }
         }
      }
   }
}

void SegTree ::saveIntersectionsIdsMulticore(REAL *bPoly, REAL *cPoly, int bSize, int cSize,
                                             int bType, int cType,
                                             vector<int> &bPolLineIds, vector<int> &cPolLineIds, vector<int> &intersectTypesDup)
{
   // void SegTree ::saveIntersectionsIdsMulticore
   //       (REAL *bPoly, REAL *cPoly, int bSize, int cSize,
   //       int bType, int cType,
   //       list<int> &bPolLineIds, list<int> &cPolLineIds, list<int> &intersectTypesDup){

   // use pointers and check

   /*int ts=treeSize, scl=sumCoverListSize*4;
   int sels=sumEndListSize;
   int bpolys=(bSize+1)*2, cpolys=(cSize+1)*2;
   double epsilon=EPSILON;
   bool bolpr;

   // #pragma omp parallel for schedule(dynamic, 100)
   for(int nid=1; nid<ts; nid++){
      int cs=coverListSizes[nid], ce=coverListSizes[nid+1];
      for(auto clid=treeNode->at(nid)->coverList->begin(); clid!=treeNode->at(nid)->coverList->end(); ++clid){
         // b \in C_B(v) and o \in C_O(v)
         // auto nclid=clid;
         // nclid++;
         for(auto clid3=treeNode->at(nid)->coverList->begin(); clid3!=treeNode->at(nid)->coverList->end(); ++clid3){
         // for(auto clid3=nclid; clid3!=treeNode->at(nid)->coverList->end(); ++clid3){
            int bId, cId;
            bolpr=false;
            // check if line segments are from different polygons
            // if (clid->type!=clid3->type){
            // If clid2 is base polygon. Else clid2 is subject polygon
            if(clid->type==bType && clid3->type==cType){
               bId=clid->id;
               cId=clid3->id;
               bolpr=true;
            }else if(clid->type==cType && clid3->type==bType){
               bId=clid3->id;
               cId=clid->id;
               bolpr=true;
            }
            if(bolpr){
               if (LSMF(bPoly[bId], bPoly[bId+1],
                        bPoly[bId+2], bPoly[bId+3],
                        cPoly[cId], cPoly[cId+1],
                        cPoly[cId+2], cPoly[cId+3])){
                  int intersecType = lineSegmentIntersection2(bPoly[bId], bPoly[bId+1],
                                                               bPoly[bId+2], bPoly[bId+3],
                                                               cPoly[cId], cPoly[cId+1],
                                                               cPoly[cId+2], cPoly[cId+3], epsilon);
                  if (intersecType){
                     // #pragma omp critical
                     {
                     bPolLineIds.push_back(bId);
                     cPolLineIds.push_back(cId);
                     intersectTypesDup.push_back(intersecType);
                     }
                  }
               }
            }
         }
      //    // b \in C_B(v) and o \in E_O(v)
         for(auto elid=treeNode->at(nid)->endEdgeSet->begin(); elid!=treeNode->at(nid)->endEdgeSet->end(); ++elid){
            int bId, cId;
            bolpr=false;
            // check if line segments are from different polygons
            // if(clid->type!=elid->type){
            // If clid2 is base polygon. Else clid2 is subject polygon
            if(clid->type==bType && elid->type==cType){
               bId=clid->id;
               cId=elid->id;
               bolpr=true;
            }else if(clid->type==cType && elid->type==bType){
               cId=clid->id;
               bId=elid->id;
               bolpr=true;
            }
            if(bolpr){
               if (LSMF(bPoly[bId], bPoly[bId+1],
                        bPoly[bId+2], bPoly[bId+3],
                        cPoly[cId], cPoly[cId+1],
                        cPoly[cId+2], cPoly[cId+3])){
                  int intersecType = lineSegmentIntersection2(bPoly[bId], bPoly[bId+1],
                                                               bPoly[bId+2], bPoly[bId+3],
                                                               cPoly[cId], cPoly[cId+1],
                                                               cPoly[cId+2], cPoly[cId+3], epsilon);
                  if (intersecType){
                     // #pragma omp critical
                     {bPolLineIds.push_back(bId);
                     cPolLineIds.push_back(cId);
                     intersectTypesDup.push_back(intersecType);
                     }
                  }
               }
            }
         }
      }
   }*/

   int ts = treeSize, scl = sumCoverListSize * 4;
   // int nid;
   int sels = sumEndListSize;
   int bpolys = (bSize + 1) * 2, cpolys = (cSize + 1) * 2;
   double epsilon = EPSILON;
   // int candidateCount=0, candidateCount2=0, candidateCount3=0;
// cout<<"Total copy items "<<bpolys+cpolys<<endl;

// #pragma omp parallel for schedule(dynamic, 1000)
#pragma omp parallel for schedule(dynamic, 100)
   for (int nid = 1; nid < ts; nid++)
   {
      // if node is null, construct it before inserting
      if(treeNode->at(nid)==nullptr)
            continue;
      double alpha, beta;
      // int clid, clid3, elid;
      // int interCount=0;
      // int cs=coverListSizes[nid], ce=coverListSizes[nid+1];
      // for(int clid=cs; clid<ce; ++clid)
      for (auto it = treeNode->at(nid)->coverList->begin(); it != treeNode->at(nid)->coverList->end(); ++it)
      {
         // int clid2=clid*4;
         // for(int clid3=clid+1; clid3<ce; ++clid3)
         // #pragma omp parallel for
         for (auto itin = treeNode->at(nid)->coverList->begin(); itin != treeNode->at(nid)->coverList->end(); ++itin)
         {
            // candidateCount++;
            int bId, cId;
            // int clid4=clid3*4;
            bool proceed = false;
            // check if line segments are from different polygons
            // if (coverArray[clid2+1]!=coverArray[clid4+1]){
            // If clid2 is base polygon. Else clid2 is subject polygon
            if (it->type == bType && itin->type == cType)
            {
               bId = it->id;
               cId = itin->id;
               proceed = true;
            }
            else if (it->type == cType && itin->type == bType)
            {
               bId = itin->id;
               cId = it->id;
               proceed = true;
            }
            if (proceed)
            {
               // candidateCount2++;
               if (LSMF(bPoly[bId], bPoly[bId + 1],
                        bPoly[bId + 2], bPoly[bId + 3],
                        cPoly[cId], cPoly[cId + 1],
                        cPoly[cId + 2], cPoly[cId + 3]))
               {
                  // candidateCount3++;
                  int intersecType = lineSegmentIntersection2(bPoly[bId], bPoly[bId + 1],
                                                              bPoly[bId + 2], bPoly[bId + 3],
                                                              cPoly[cId], cPoly[cId + 1],
                                                              cPoly[cId + 2], cPoly[cId + 3], epsilon);
                  // cout<<"nid="<<nid<<" case1 bid: "<<bId/2<<" cid: "<<cId/2<<" type: "<<intersecType<<endl;
                  if (intersecType)
                  {
                     #pragma omp critical
                     // #pragma omp task
                     {
                        // try{
                        bPolLineIds.push_back(bId);
                        cPolLineIds.push_back(cId);
                        intersectTypesDup.push_back(intersecType);
                        // }catch(const std::bad_alloc &) {
                        //    // ec=1;
                        // }
                     }
                  }
               }
            }
         }

         int es = endListSizes[nid], ee = endListSizes[nid] + endListCounts[nid];
         // #pragma omp parallel for
         for (int elid = es; elid < ee; ++elid)
         {
            // candidateCount++;
            int bId, cId;
            int elid2 = elid * 4;
            bool proceed = false;
            // check if line segments are from different polygons
            if (it->type == bType && treeNodeEndListArray[elid].type == cType)
            {
               bId = it->id;
               // cId=endArray[elid2];
               cId = treeNodeEndListArray[elid].id;
               proceed = true;
            }
            else if (it->type == cType && treeNodeEndListArray[elid].type == bType)
            {
               cId = it->id;
               // bId=endArray[elid2];
               bId = treeNodeEndListArray[elid].id;
               proceed = true;
            }
            if (proceed)
            {
               // candidateCount2++;
               // cout<<nid<<" case2 bid: "<<bId/2<<" cid: "<<cId/2<<endl;
               if (LSMF(bPoly[bId], bPoly[bId + 1],
                        bPoly[bId + 2], bPoly[bId + 3],
                        cPoly[cId], cPoly[cId + 1],
                        cPoly[cId + 2], cPoly[cId + 3]))
               {
                  // candidateCount3++;
                  int intersecType = lineSegmentIntersection2(bPoly[bId], bPoly[bId + 1],
                                                              bPoly[bId + 2], bPoly[bId + 3],
                                                              cPoly[cId], cPoly[cId + 1],
                                                              cPoly[cId + 2], cPoly[cId + 3], epsilon);
                  // cout<<"nid="<<nid<<" case2 bid: "<<bId/2<<" cid: "<<cId/2<<" type: "<<intersecType<<endl;
                  if (intersecType)
                  {
                     #pragma omp critical
                     // #pragma omp task
                     {
                        // try{
                        bPolLineIds.push_back(bId);
                        cPolLineIds.push_back(cId);
                        intersectTypesDup.push_back(intersecType);
                        // }catch(const std::bad_alloc &) {
                        //    // ec=1;
                        // }
                     }
                  }
               }
            }
         }
      }
   }
   // cout<<"***********candidateCount= "<<candidateCount<<" candidateCount2="<<candidateCount2<<" candidateCount3="<<candidateCount3<<endl;
}

/*
Count intersections parallel
*/

void SegTree ::countIntersectionsParallel(REAL *bPoly, REAL *cPoly, int bSize, int cSize, int cPolyCount,
                                          int bType, int cType,
                                          int *iCount)
{

   int ts = treeSize, scl = sumCoverListSize * 4;
   // int nid;
   int sels = sumEndListSize;
   int bpolys = (bSize + 1) * 2, cpolys = (cSize + cPolyCount) * 2;
   double epsilon = EPSILON;
   // cout<<"---- "<<*(*bPoly+4)<<" "<<*(*bPoly+5)<<endl;
   // cout<<"---- "<<bPoly[0][4]<<" "<<bPoly[0][5]<<endl;


   // cout<<"Total copy items "<<bpolys+cpolys<<endl;
   // #ifdef DG2
   // high_resolution_clock::time_point cp1, cp2;
   // cp1 = high_resolution_clock::now();
   // #endif
   // #pragma acc update device(endListCounts[0:ts])
   // #pragma acc update device(treeNodeEndListArray[0 : sels], endListSizes[0 : ts], endListCounts[0 : ts])

   // for device
   #pragma acc enter data create(iCount[0 : ts]) /*copyin(bPoly[0:(bSize+1)*2], cPoly[0:(cSize+1)*2])*/
   // #pragma acc enter data copyin(bPoly[0 : bpolys], cPoly[0 : cpolys], epsilon)

   // #pragma acc data present(coverArray [0:scl], treeNodeEndListArray[0:sels], coverListSizes [0:ts+1], endListSizes[0:ts+1], endListCounts[0:ts])
   #pragma acc data present(coverArray[0 : scl], treeNodeEndListArray[0 : sels], coverListSizes[0 : ts + 1], endListSizes[0 : ts + 1], endListCounts[0:ts], nodeIntervalArray[0 : ts * 2], bPoly[0 : bpolys], cPoly[0 : cpolys])

   // #ifdef DG2
   // cp2 = high_resolution_clock::now();
   // #endif

   #pragma acc parallel loop
   for (int nid = 1; nid < ts; nid++)
   {
      double alpha, beta;
      // int clid, clid3, elid;
      int interCount = 0;
      int cs = coverListSizes[nid], ce = coverListSizes[nid + 1];

      // for(int i=cs; i<ce; ++i){
      //    interCount+=1;
      // }
      
      // #pragma acc for async
      // for(clid=coverListSizes[nid], clid2=clid*4; clid<coverListSizes[nid+1]; ++clid, clid2+=4)
      for (int clid = cs; clid < ce; ++clid)
      {
         int clid2 = clid * 4;
         // b \in C_B(v) and o \in C_O(v)
         // for(clid3=clid+1, clid4=clid3*4; clid3<coverListSizes[nid+1]; ++clid3, clid4+=4)
         for (int clid3 = clid + 1; clid3 < ce; ++clid3)
         {
            int bId, cId;
            int clid4 = clid3 * 4;
            bool proceed = false;
            // check if line segments are from different polygons
            // if (coverArray[clid2+1]!=coverArray[clid4+1]){
            // If clid2 is base polygon. Else clid2 is subject polygon
            if (coverArray[clid2 + 1] == bType && coverArray[clid4 + 1] == cType)
            {
               bId = coverArray[clid2];
               cId = coverArray[clid4];
               proceed = true;
            }
            else if (coverArray[clid2 + 1] == cType && coverArray[clid4 + 1] == bType)
            {
               bId = coverArray[clid4];
               cId = coverArray[clid2];
               proceed = true;
            }
            if (proceed)
            {
               // cout<<nid<<" case1 bid: "<<bId/2<<" cid: "<<cId/2<<endl;

               // if (LSMF(*(*bPoly + (bId)), *(*bPoly + (bId + 1)),
               //          *(*bPoly + (bId + 2)), *(*bPoly + (bId + 3)),
               //          *(*cPoly + (cId)), *(*cPoly + (cId + 1)),
               //          *(*cPoly + (cId + 2)), *(*cPoly + (cId + 3)))){
               //    // cout<<"nid="<<nid<<" case1 bid: "<<bId/2<<" cid: "<<cId/2<<endl;
               //    int intersecType = lineSegmentIntersection
               //                         (*(*bPoly + (bId)), *(*bPoly + (bId + 1)),
               //                         *(*bPoly + (bId + 2)), *(*bPoly + (bId + 3)),
               //                         *(*cPoly + (cId)), *(*cPoly + (cId + 1)),
               //                         *(*cPoly + (cId + 2)), *(*cPoly + (cId + 3)), &alpha, &beta);

               if (LSMF(bPoly[bId], bPoly[bId + 1],
                        bPoly[bId + 2], bPoly[bId + 3],
                        cPoly[cId], cPoly[cId + 1],
                        cPoly[cId + 2], cPoly[cId + 3]))
               {
                  int intersecType = lineSegmentIntersection2(bPoly[bId], bPoly[bId + 1],
                                                              bPoly[bId + 2], bPoly[bId + 3],
                                                              cPoly[cId], cPoly[cId + 1],
                                                              cPoly[cId + 2], cPoly[cId + 3], epsilon);
                  // cout<<"nid="<<nid<<" case1 bid: "<<bId/2<<" cid: "<<cId/2<<" type: "<<intersecType<<endl;
                  if (intersecType)
                  {
                     interCount += 1;
                  }
                  // /*REAL intersectX, intersectY;
                  // point b1, b2, c1, I;
                  // // int bId=bPolyLineIds[id];
                  // b1.x=bPoly[bId];
                  // b1.y=bPoly[bId+1];
                  // b2.x=bPoly[bId+2];
                  // b2.y=bPoly[bId+3];

                  // // int cId=cPolyLineIds[id];
                  // c1.x=cPoly[cId];
                  // c1.y=cPoly[cId+1];

                  // if(intersecType){
                  //    // if((intersecType==1 || intersecType==3 || intersecType==5 || intersecType==7)){
                  //    //    nonDegenCount++;
                  //    //    count2=nonDegenCount;
                  //    // }
                  //    // else if((intersecType==2 || intersecType==4 || intersecType==6 || intersecType==8)){
                  //    //    count2=0;
                  //    // }

                  //    switch(intersecType){
                  //       // case X_INTERSECTION:
                  //       // I and I
                  //       case 1:
                  //          I = add(mulScalar((1.0-alpha), b1), mulScalar(alpha, b2));
                  //          intersectX=I.x;       //consider edge for the intersection array
                  //          intersectY=I.y;
                  //          break;
                  //       case 3:
                  //       case 4:
                  //       case 5:
                  //       case 7:
                  //       case 8:
                  //          intersectX=c1.x;
                  //          intersectY=c1.y;
                  //          break;
                  //       case 2:
                  //       case 6:
                  //          intersectX=b1.x;
                  //          intersectY=b1.y;
                  //       break;

                  //    }
                  //    // if(isContain(intersectX, intersectY, nodeIntervalArray[nid*2], nodeIntervalArray[nid*2+1])){
                  //    if(isContain(nodeIntervalArray[nid*2], nodeIntervalArray[nid*2+1], intersectX, intersectY)){
                  //       interCount+=1;
                  //    }
                  // } */
               }
            }
         }
         //    // b \in C_B(v) and o \in E_O(v)
         //    // for(elid=endListSizes[nid], elid2=elid*4; elid<(endListSizes[nid]+endListCounts[nid]); ++elid, elid2+=4)

         int es = endListSizes[nid], ee = /*endListSizes[nid+1]*/ endListSizes[nid] + endListCounts[nid];
         // #pragma acc for async
         for (int elid = es; elid < ee; ++elid)
         {
            int bId, cId;
            int elid2 = elid * 4;
            bool proceed = false;
            // check if line segments are from different polygons
            // if(coverArray[clid2+1]!=endArray[elid2+1]){
            // if(coverArray[clid2+1]!=treeNodeEndListArray[elid].type){
            // If clid2 is base polygon. Else clid2 is subject polygon
            if (coverArray[clid2 + 1] == bType && treeNodeEndListArray[elid].type == cType)
            {
               bId = coverArray[clid2];
               // cId=endArray[elid2];
               cId = treeNodeEndListArray[elid].id;
               proceed = true;
            }
            else if (coverArray[clid2 + 1] == cType && treeNodeEndListArray[elid].type == bType)
            {
               cId = coverArray[clid2];
               // bId=endArray[elid2];
               bId = treeNodeEndListArray[elid].id;
               proceed = true;
            }
            if (proceed)
            {
               // cout<<nid<<" case2 bid: "<<bId/2<<" cid: "<<cId/2<<endl;
               if (LSMF(bPoly[bId], bPoly[bId + 1],
                        bPoly[bId + 2], bPoly[bId + 3],
                        cPoly[cId], cPoly[cId + 1],
                        cPoly[cId + 2], cPoly[cId + 3]))
               {
                  int intersecType = lineSegmentIntersection2(bPoly[bId], bPoly[bId + 1],
                                             bPoly[bId + 2], bPoly[bId + 3],
                                             cPoly[cId], cPoly[cId + 1],
                                             cPoly[cId + 2], cPoly[cId + 3], epsilon);

                  // if (LSMF(bPoly[0][bId], bPoly[0][bId+1],
                  //          bPoly[0][bId+2], bPoly[0][bId+3],
                  //          cPoly[0][cId], cPoly[0][cId+1],
                  //          cPoly[0][cId+2], cPoly[0][cId+3])){
                  //    int intersecType = lineSegmentIntersection(bPoly[0][bId], bPoly[0][bId+1],
                  //                                                 bPoly[0][bId+2], bPoly[0][bId+3],
                  //                                                 cPoly[0][cId], cPoly[0][cId+1],
                  //                                                 cPoly[0][cId+2], cPoly[0][cId+3], &alpha, &beta);
                  // if (LSMF(*(*bPoly + (bId)), *(*bPoly + (bId + 1)),
                  //          *(*bPoly + (bId + 2)), *(*bPoly + (bId + 3)),
                  //          *(*cPoly + (cId)), *(*cPoly + (cId + 1)),
                  //          *(*cPoly + (cId + 2)), *(*cPoly + (cId + 3)))){
                  // int intersecType = lineSegmentIntersection(*(*bPoly + (bId)), *(*bPoly + (bId + 1)),
                  //                                              *(*bPoly + (bId + 2)), *(*bPoly + (bId + 3)),
                  //                                              *(*cPoly + (cId)), *(*cPoly + (cId + 1)),
                  //                                              *(*cPoly + (cId + 2)), *(*cPoly + (cId + 3)), &alpha, &beta);
                  if (intersecType)
                  {
                     interCount += 1;
                  }

                  // REAL intersectX, intersectY;
                  // point b1, b2, c1, I;
                  // // int bId=bPolyLineIds[id];
                  // b1.x=bPoly[bId];
                  // b1.y=bPoly[bId+1];
                  // b2.x=bPoly[bId+2];
                  // b2.y=bPoly[bId+3];

                  // // int cId=cPolyLineIds[id];
                  // c1.x=cPoly[cId];
                  // c1.y=cPoly[cId+1];

                  // if(intersecType){
                  //    // if((intersecType==1 || intersecType==3 || intersecType==5 || intersecType==7)){
                  //    //    nonDegenCount++;
                  //    //    count2=nonDegenCount;
                  //    // }
                  //    // else if((intersecType==2 || intersecType==4 || intersecType==6 || intersecType==8)){
                  //    //    count2=0;
                  //    // }

                  //    switch(intersecType){
                  //       // case X_INTERSECTION:
                  //       // I and I
                  //       case 1:
                  //          I = add(mulScalar((1.0-alpha), b1), mulScalar(alpha, b2));
                  //          intersectX=I.x;       //consider edge for the intersection array
                  //          intersectY=I.y;
                  //          break;
                  //       case 3:
                  //       case 4:
                  //       case 5:
                  //       case 7:
                  //       case 8:
                  //          intersectX=c1.x;
                  //          intersectY=c1.y;
                  //          break;
                  //       case 2:
                  //       case 6:
                  //          intersectX=b1.x;
                  //          intersectY=b1.y;
                  //       break;

                  //    }
                  //    // if(isContain(intersectX, intersectY, nodeIntervalArray[nid*2], nodeIntervalArray[nid*2+1])){
                  //    if(isContain(nodeIntervalArray[nid*2], nodeIntervalArray[nid*2+1], intersectX, intersectY)){
                  //       interCount+=1;
                  //    }
                  // }
               }
            }
         }
      }
      iCount[nid - 1] = interCount;
   }

   // #ifdef DG2
   //    // Calculating total time taken by the program.
   //    auto d1 = duration_cast<microseconds>(cp2 - cp1);
   //    cout << "countIntersect: to GPU copy() " << fixed << d1.count() << setprecision(10) << endl;
   // #endif
}

/*
GPU parallel
*/

void SegTree ::polyClipArrayGPUParallel
      (REAL *bPoly, REAL *cPoly, int bSize, int cSize,
      int bType, int cType,
      int *iCount, int *bPolLineIds, int *cPolLineIds, int *intersectTypes){

   int ts = treeSize, scl=sumCoverListSize*4;
   int sels=sumEndListSize;
   int intc=iCount[ts-1];
   int bpolys=(bSize+1)*2, cpolys=(cSize+1)*2;

   int nid, clid;
   double epsilon=EPSILON;
   // int interCount=0;
   // int cs=coverListSizes[nid], ce=coverListSizes[nid+1];
   // int es=endListSizes[nid], ee=endListSizes[nid]+endListCounts[nid];

   #ifdef TIME_L2
   high_resolution_clock::time_point cp1, cp2;
   #endif

   // #pragma acc data present(coverArray [0:scl], treeNodeEndListArray[0:sels], coverListSizes [0:ts+1], endListSizes[0:ts+1], endListCounts[0:ts])

   #pragma acc enter data create(bPolLineIds[0:intc], cPolLineIds[0:intc], intersectTypes[0:intc])
   #pragma acc data present(coverArray [0:scl], treeNodeEndListArray[0:sels], coverListSizes [0:ts+1], endListSizes[0:ts+1], endListCounts[0:ts], bPoly[0:bpolys], cPoly[0:cpolys], iCount[0 : ts])
   #pragma acc parallel loop
   for(nid=1; nid<ts; nid++)
   {
      int interCount=0;
      int einterCount=iCount[nid]-iCount[nid-1]-1;
      for(clid=coverListSizes[nid]; clid<coverListSizes[nid+1]; ++clid)
      {
         int  clid2=clid*4;
         // b \in C_B(v) and o \in C_O(v)
         // #pragma acc for async
         // #pragma omp for nowait
         // #pragma omp parallel for
         for(int clid3=clid+1; clid3<coverListSizes[nid+1]; ++clid3)
         {
            int bId, cId;
            int clid4=clid3*4;
            bool proceed=false;
            // check if line segments are from different polygons
            // if(coverArray[clid2+1]!=coverArray[clid4+1]){
            // If clid2 is base polygon. Else clid2 is subject polygon
            if(coverArray[clid2+1]==bType && coverArray[clid4+1]==cType){
               bId=coverArray[clid2];
               cId=coverArray[clid4];
               proceed=true;
            }else if (coverArray[clid2+1]==cType && coverArray[clid4+1]==bType){
               bId=coverArray[clid4];
               cId=coverArray[clid2];
               proceed=true;
            }
            if(proceed){
               // cout<<"case1 bid: "<<bId/2<<" cid: "<<cId/2<<endl;
               // double alpha, beta;
               if (LSMF(bPoly[bId], bPoly[bId+1],
                        bPoly[bId+2], bPoly[bId+3],
                        cPoly[cId], cPoly[cId+1],
                        cPoly[cId+2], cPoly[cId+3]))
               {
                  int intersecType = lineSegmentIntersection2(bPoly[bId], bPoly[bId+1],
                                                               bPoly[bId+2], bPoly[bId+3],
                                                               cPoly[cId], cPoly[cId+1],
                                                               cPoly[cId+2], cPoly[cId+3], epsilon);
                  // int intersecType=1;
                  if (intersecType){
                     // #pragma acc atomic
                     {
                        bPolLineIds[iCount[nid-1]+interCount]=bId; 
                        cPolLineIds[iCount[nid-1]+interCount]=cId; 
                        intersectTypes[iCount[nid-1]+interCount]=intersecType;
                        interCount+=1;
                     }
                  }
               }
            }
         }
         // b \in C_B(v) and o \in E_O(v)
         // #pragma acc for async
         // #pragma omp for nowait
         // #pragma omp parallel for
         for(int elid=endListSizes[nid]; elid<endListSizes[nid]+endListCounts[nid]; ++elid)
         // for(int elid=endListSizes[nid]; elid<endListSizes[nid+1]; ++elid)
         {
            int bId, cId;
            int elid2=elid*4;
            bool proceed=false;
            // check if line segments are from different polygons
            // if(coverArray[clid2+1]!=treeNodeEndListArray[elid].type){
            // If clid2 is base polygon. Else clid2 is subject polygon
            if(coverArray[clid2+1]==bType && treeNodeEndListArray[elid].type==cType){
               bId=coverArray[clid2];
               cId=treeNodeEndListArray[elid].id;
               proceed=true;
            }else if(coverArray[clid2+1]==cType && treeNodeEndListArray[elid].type==bType){
               cId=coverArray[clid2];
               bId=treeNodeEndListArray[elid].id;
               proceed=true;
            }
            if(proceed){
               // double alpha, beta;
               // cout<<nid<<" case2 bid: "<<bId/2<<" cid: "<<cId/2<<endl;
               if (LSMF(bPoly[bId], bPoly[bId+1],
                        bPoly[bId+2], bPoly[bId+3],
                        cPoly[cId], cPoly[cId+1],
                        cPoly[cId+2], cPoly[cId+3]))
               {
                  int intersecType = lineSegmentIntersection2(bPoly[bId], bPoly[bId+1],
                                                               bPoly[bId+2], bPoly[bId+3],
                                                               cPoly[cId], cPoly[cId+1],
                                                               cPoly[cId+2], cPoly[cId+3], epsilon);
                  
                  // int intersecType=1;
                  if (intersecType){
                     // #pragma acc atomic
                     {
                        bPolLineIds[iCount[nid-1]+einterCount]=bId; 
                        cPolLineIds[iCount[nid-1]+einterCount]=cId; 
                        intersectTypes[iCount[nid-1]+einterCount]=intersecType;
                        einterCount-=1;
                     }
                     // {
                     //    bPolLineIds[iCount[nid-1]+interCount]=bId;
                     //    cPolLineIds[iCount[nid-1]+interCount]=cId;
                     //    intersectTypes[iCount[nid-1]+interCount]=intersecType;
                     //    interCount+=1;
                     // }
                     // cout<<"*********-----------**** nid="<<nid<<" bid="<<bId<<" cId="<<cId<<endl;
                  }
               }
            }
         }
      }
   }

   #ifdef TIME_L2
      cp1 = high_resolution_clock::now();
   #endif
   #pragma acc update self(bPolLineIds[0:intc], cPolLineIds[0:intc], intersectTypes[0:intc])
   #ifdef TIME_L2
      cp2 = high_resolution_clock::now();
      auto d1 = duration_cast<microseconds>(cp2 - cp1);
      #if DT_INFO
         cout<< " |-Time: interscting edge IDs and types copy to CPU: " << fixed << d1.count() << setprecision(10) << endl;
      #else
         cout<< fixed << d1.count() << setprecision(10) << " ";
      #endif
   #endif
}

/*
Count intersections parallel: Multicore
*/
void SegTree ::countIntersectionsMulticore(REAL *bPoly, REAL *cPoly, int bSize, int cSize,
                                           int bType, int cType,
                                           int *iCount)
{

   int ts = treeSize, scl = sumCoverListSize * 4;
   // int nid;
   int sels = sumEndListSize;
   int bpolys = (bSize + 1) * 2, cpolys = (cSize + 1) * 2;
   double epsilon = EPSILON;


   // cout<<"Total copy items "<<bpolys+cpolys<<endl;
   #ifdef DG2
   high_resolution_clock::time_point cp1, cp2;
   cp1 = high_resolution_clock::now();
   #endif

   #ifdef DG2
   cp2 = high_resolution_clock::now();
   #endif

   #pragma omp parallel for schedule(dynamic, 100)
   for (int nid = 1; nid < ts; nid++)
   {
      double alpha, beta;
      // int clid, clid3, elid;
      int interCount = 0;
      int cs = coverListSizes[nid], ce = coverListSizes[nid + 1];

      // for(int i=cs; i<ce; ++i){
      //    interCount+=1;
      // }

      // for(clid=coverListSizes[nid], clid2=clid*4; clid<coverListSizes[nid+1]; ++clid, clid2+=4)
      for (int clid = cs; clid < ce; ++clid)
      {
         int clid2 = clid * 4;
         // b \in C_B(v) and o \in C_O(v)
         // for(clid3=clid+1, clid4=clid3*4; clid3<coverListSizes[nid+1]; ++clid3, clid4+=4)
         for (int clid3 = clid + 1; clid3 < ce; ++clid3)
         {
            int bId, cId;
            int clid4 = clid3 * 4;
            bool proceed = false;
            // check if line segments are from different polygons
            // if (coverArray[clid2+1]!=coverArray[clid4+1]){
            // If clid2 is base polygon. Else clid2 is subject polygon
            if (coverArray[clid2 + 1] == bType && coverArray[clid4 + 1] == cType)
            {
               bId = coverArray[clid2];
               cId = coverArray[clid4];
               proceed = true;
            }
            else if (coverArray[clid2 + 1] == cType && coverArray[clid4 + 1] == bType)
            {
               bId = coverArray[clid4];
               cId = coverArray[clid2];
               proceed = true;
            }
            if (proceed)
            {
               if (LSMF(bPoly[bId], bPoly[bId + 1],
                        bPoly[bId + 2], bPoly[bId + 3],
                        cPoly[cId], cPoly[cId + 1],
                        cPoly[cId + 2], cPoly[cId + 3]))
               {
                  int intersecType = lineSegmentIntersection2(bPoly[bId], bPoly[bId + 1],
                                                              bPoly[bId + 2], bPoly[bId + 3],
                                                              cPoly[cId], cPoly[cId + 1],
                                                              cPoly[cId + 2], cPoly[cId + 3], epsilon);
                  // cout<<"nid="<<nid<<" case1 bid: "<<bId/2<<" cid: "<<cId/2<<" type: "<<intersecType<<endl;
                  if (intersecType)
                  {
                     interCount += 1;
                  }
               }
            }
         }

         int es = endListSizes[nid], ee = endListSizes[nid] + endListCounts[nid];
         for (int elid = es; elid < ee; ++elid)
         {
            int bId, cId;
            int elid2 = elid * 4;
            bool proceed = false;
            // check if line segments are from different polygons
            if (coverArray[clid2 + 1] == bType && treeNodeEndListArray[elid].type == cType)
            {
               bId = coverArray[clid2];
               // cId=endArray[elid2];
               cId = treeNodeEndListArray[elid].id;
               proceed = true;
            }
            else if (coverArray[clid2 + 1] == cType && treeNodeEndListArray[elid].type == bType)
            {
               cId = coverArray[clid2];
               // bId=endArray[elid2];
               bId = treeNodeEndListArray[elid].id;
               proceed = true;
            }
            if (proceed)
            {
               // cout<<nid<<" case2 bid: "<<bId/2<<" cid: "<<cId/2<<endl;
               if (LSMF(bPoly[bId], bPoly[bId + 1],
                        bPoly[bId + 2], bPoly[bId + 3],
                        cPoly[cId], cPoly[cId + 1],
                        cPoly[cId + 2], cPoly[cId + 3]))
               {
                  int intersecType = lineSegmentIntersection2(bPoly[bId], bPoly[bId + 1],
                                                              bPoly[bId + 2], bPoly[bId + 3],
                                                              cPoly[cId], cPoly[cId + 1],
                                                              cPoly[cId + 2], cPoly[cId + 3], epsilon);
                  if (intersecType)
                  {
                     interCount += 1;
                  }
               }
            }
         }
      }
      iCount[nid - 1] = interCount;
   }

   #ifdef DG2
      // Calculating total time taken by the program.
      auto d1 = duration_cast<microseconds>(cp2 - cp1);

      cout << "countIntersect: to GPU copy() " << fixed << d1.count() << setprecision(10) << endl;
   #endif
}

/*
intersection betwwen edges in each node
*/
// -------------------------
// void SegTree ::polyClipArrayParallel
//       (REAL *bPoly, REAL *cPoly, int bSize, int cSize,
//       int *iCount, int *bPolLineIds, int *cPolLineIds, int *intersectTypes){
// // void SegTree ::polyClipArrayParallel
// //       (REAL *bPoly, REAL *cPoly, int bSize, int cSize,
// //       int *iCount, dsintersection *dsint){

//    int ts = treeSize, scl=sumCoverListSize*4;
//    int sels=sumEndListSize;
//    int intc=iCount[ts-1];
//    int bpolys=(bSize+1)*2, cpolys=(cSize+1)*2;

//    /*int *interCountArray;
//    interCountArray=(int *)malloc((ts)*sizeof(int));
//    for(int i=0; i<ts; ++i){
//       interCountArray[i]=0;
//    }*/

//    // for host
//    // #pragma acc update self(coverArray [0:sumEndListSize*4], endArray[0:sumEndListSize*4])

//    // for device
//    // #pragma acc enter data create(intsections[0:intc]) /*copyin(interCountArray[0:ts])*/

//    // #pragma acc enter data create(bPolLineIds[0:intc], cPolLineIds[0:intc], intersectTypes[0:intc])
//    // #pragma acc parallel loop present(coverArray [0:scl], treeNodeEndListArray[0:sels], coverListSizes[0:ts+1], endListSizes[0:ts+1], endListCounts[0:ts], iCount[0:ts], bPoly[0:bpolys], cPoly[0:cpolys])

//    int nid, clid;
//    // omp_set_dynamic(0);
//    // omp_set_num_threads(20);
//    #pragma omp parallel for private(clid)
//    // #pragma omp parallel for shared(coverArray, treeNodeEndListArray, coverListSizes, endListSizes, endListCounts, iCount, bPoly, cPoly, bPolLineIds, cPolLineIds, intersectTypes)
//    for(nid=1; nid<ts; nid++)
//    {
//       int interCount=0;
//       int cs=coverListSizes[nid], ce=coverListSizes[nid+1];

//       // interCountArray[nid-1]=0;

//       // #pragma acc loop
//       // #pragma omp parallel for
//       for(clid=cs; clid<ce; ++clid)
//       {
//          int  clid2=clid*4;
//          // b \in C_B(v) and o \in C_O(v)
//          // #pragma omp for nowait
//          // #pragma omp parallel for
//          for(int clid3=clid+1; clid3<ce; ++clid3)
//          {
//             int bId, cId;
//             int clid4=clid3*4;
//             // check if line segments are from different polygons
//             if(coverArray[clid2+1]!=coverArray[clid4+1]){
//                // If clid2 is base polygon. Else clid2 is subject polygon
//                if(coverArray[clid2+1]==0){
//                   bId=coverArray[clid2];
//                   cId=coverArray[clid4];
//                }else{
//                   bId=coverArray[clid4];
//                   cId=coverArray[clid2];
//                }
//                // cout<<"case1 bid: "<<bId/2<<" cid: "<<cId/2<<endl;
//                double alpha, beta;
//                if (LSMF(bPoly[bId], bPoly[bId+1],
//                         bPoly[bId+2], bPoly[bId+3],
//                         cPoly[cId], cPoly[cId+1],
//                         cPoly[cId+2], cPoly[cId+3]))
//                {
//                   int intersecType = lineSegmentIntersection(bPoly[bId], bPoly[bId+1],
//                                                                bPoly[bId+2], bPoly[bId+3],
//                                                                cPoly[cId], cPoly[cId+1],
//                                                                cPoly[cId+2], cPoly[cId+3], &alpha, &beta);
//                   // int intersecType=1;
//                   if (intersecType){
//                      // dsintersection temp;
//                      // temp.type=intersecType;
//                      // temp.bPolId=bId;
//                      // temp.cPolId=cId;
//                      // #pragma acc atomic update
//                      // dsint[iCount[nid-1]+(interCountArray[nid]+=1)]=temp;
//                      // dsint[iCount[nid-1]+interCount]=temp;
//                      // dsint[iCount[nid-1]]=temp;
//                      #pragma omp critical
//                      {
//                         bPolLineIds[iCount[nid-1]+interCount]=bId;
//                         cPolLineIds[iCount[nid-1]+interCount]=cId;
//                         intersectTypes[iCount[nid-1]+interCount]=intersecType;
//                         interCount+=1;
//                      }
//                      // cout<<"*********-----------**** nid="<<nid<<" bid="<<bId<<" cId="<<cId<<endl;
//                   }

//                   /*REAL intersectX, intersectY;
//                   point b1, b2, c1, I;
//                   // int bId=bPolyLineIds[id];
//                   b1.x=bPoly[bId];
//                   b1.y=bPoly[bId+1];
//                   b2.x=bPoly[bId+2];
//                   b2.y=bPoly[bId+3];

//                   // int cId=cPolyLineIds[id];
//                   c1.x=cPoly[cId];
//                   c1.y=cPoly[cId+1];

//                   if(intersecType){
//                      // if((intersecType==1 || intersecType==3 || intersecType==5 || intersecType==7)){
//                      //    nonDegenCount++;
//                      //    count2=nonDegenCount;
//                      // }
//                      // else if((intersecType==2 || intersecType==4 || intersecType==6 || intersecType==8)){
//                      //    count2=0;
//                      // }

//                      switch(intersecType){
//                         // case X_INTERSECTION:
//                         // I and I
//                         case 1:
//                            I = add(mulScalar((1.0-alpha), b1), mulScalar(alpha, b2));
//                            intersectX=I.x;       //consider edge for the intersection array
//                            intersectY=I.y;
//                            break;
//                         case 3:
//                         case 4:
//                         case 5:
//                         case 7:
//                         case 8:
//                            intersectX=c1.x;
//                            intersectY=c1.y;
//                            break;
//                         case 2:
//                         case 6:
//                            intersectX=b1.x;
//                            intersectY=b1.y;
//                         break;

//                      }
//                      // if(isContain(intersectX, intersectY, nodeIntervalArray[nid*2], nodeIntervalArray[nid*2+1])){
//                      if(isContain(nodeIntervalArray[nid*2], nodeIntervalArray[nid*2+1], intersectX, intersectY)){
//                         #pragma omp critical
//                         {
//                            bPolLineIds[iCount[nid-1]+interCount]=bId;
//                            cPolLineIds[iCount[nid-1]+interCount]=cId;
//                            intersectTypes[iCount[nid-1]+interCount]=intersecType;
//                            interCount+=1;
//                         }
//                      }
//                   }*/

//                }
//             }
//          }
//          // b \in C_B(v) and o \in E_O(v)
//          int es=endListSizes[nid], ee=endListSizes[nid]+endListCounts[nid];
//          // #pragma omp for nowait
//          // #pragma omp parallel for
//          for(int elid=es; elid<ee; ++elid)
//          {
//             int bId, cId;
//             int elid2=elid*4;
//             // check if line segments are from different polygons
//             // if(coverArray[clid2+1]!=endArray[elid2 + 1]){
//             if(coverArray[clid2+1]!=treeNodeEndListArray[elid].type){
//                // If clid2 is base polygon. Else clid2 is subject polygon
//                if(coverArray[clid2+1]==0){
//                   bId=coverArray[clid2];
//                   // cId=endArray[elid2];
//                   cId=treeNodeEndListArray[elid].id;

//                }else{
//                   cId=coverArray[clid2];
//                   // bId=endArray[elid2];
//                   bId=treeNodeEndListArray[elid].id;
//                }
//                double alpha, beta;
//                // cout<<nid<<" case2 bid: "<<bId/2<<" cid: "<<cId/2<<endl;
//                if (LSMF(bPoly[bId], bPoly[bId+1],
//                         bPoly[bId+2], bPoly[bId+3],
//                         cPoly[cId], cPoly[cId+1],
//                         cPoly[cId+2], cPoly[cId+3]))
//                {
//                   int intersecType = lineSegmentIntersection(bPoly[bId], bPoly[bId+1],
//                                                                bPoly[bId+2], bPoly[bId+3],
//                                                                cPoly[cId], cPoly[cId+1],
//                                                                cPoly[cId+2], cPoly[cId+3], &alpha, &beta);
//                   // int intersecType=1;
//                   if (intersecType){
//                      // dsintersection temp;
//                      // temp.type=intersecType;
//                      // temp.bPolId=bId;
//                      // temp.cPolId=cId;
//                      // #pragma acc atomic capture
//                      // {
//                      //    dsint[iCount[nid-1]+interCount]=temp;
//                      //    interCount += 1;
//                      // }
//                      // #pragma acc atomic capture
//                      // {
//                      #pragma omp critical
//                      {
//                         bPolLineIds[iCount[nid-1]+interCount]=bId;
//                         cPolLineIds[iCount[nid-1]+interCount]=cId;
//                         intersectTypes[iCount[nid-1]+interCount]=intersecType;
//                         interCount+=1;
//                      }
//                      // cout<<"*********-----------**** nid="<<nid<<" bid="<<bId<<" cId="<<cId<<endl;
//                   }

//                   /*REAL intersectX, intersectY;
//                   point b1, b2, c1, I;
//                   // int bId=bPolyLineIds[id];
//                   b1.x=bPoly[bId];
//                   b1.y=bPoly[bId+1];
//                   b2.x=bPoly[bId+2];
//                   b2.y=bPoly[bId+3];

//                   // int cId=cPolyLineIds[id];
//                   c1.x=cPoly[cId];
//                   c1.y=cPoly[cId+1];

//                   if(intersecType){
//                      // if((intersecType==1 || intersecType==3 || intersecType==5 || intersecType==7)){
//                      //    nonDegenCount++;
//                      //    count2=nonDegenCount;
//                      // }
//                      // else if((intersecType==2 || intersecType==4 || intersecType==6 || intersecType==8)){
//                      //    count2=0;
//                      // }

//                      switch(intersecType){
//                         // case X_INTERSECTION:
//                         // I and I
//                         case 1:
//                            I = add(mulScalar((1.0-alpha), b1), mulScalar(alpha, b2));
//                            intersectX=I.x;       //consider edge for the intersection array
//                            intersectY=I.y;
//                            break;
//                         case 3:
//                         case 4:
//                         case 5:
//                         case 7:
//                         case 8:
//                            intersectX=c1.x;
//                            intersectY=c1.y;
//                            break;
//                         case 2:
//                         case 6:
//                            intersectX=b1.x;
//                            intersectY=b1.y;
//                         break;

//                      }
//                      // if(isContain(intersectX, intersectY, nodeIntervalArray[nid*2], nodeIntervalArray[nid*2+1])){
//                      if(isContain(nodeIntervalArray[nid*2], nodeIntervalArray[nid*2+1], intersectX, intersectY)){
//                         #pragma omp critical
//                         {
//                            bPolLineIds[iCount[nid-1]+interCount]=bId;
//                            cPolLineIds[iCount[nid-1]+interCount]=cId;
//                            intersectTypes[iCount[nid-1]+interCount]=intersecType;
//                            interCount+=1;
//                         }
//                      }
//                   }*/

//                }
//             }
//          }
//       }
//    }
//    // high_resolution_clock::time_point cp1, cp2;
//    // if(DEBUG_TIME) cp1=high_resolution_clock::now();

//    // #pragma acc update self(bPolLineIds[0:intc], cPolLineIds[0:intc], intersectTypes[0:intc])

//    // if(DEBUG_TIME){
//    //    cp2=high_resolution_clock::now();
//    //    auto d1 = duration_cast<microseconds>(cp2-cp1);
//    //    cout<<" copy saved ids time: "<<d1.count()<<endl;
//    // }
// }
//-----------------------------

// cleaned version for OpenMP
void SegTree ::polyClipArrayParallel(REAL *bPoly, REAL *cPoly, int bSize, int cSize,
                                     int bType, int cType,
                                     int *iCount, int *bPolLineIds, int *cPolLineIds, int *intersectTypes)
{

   int ts = treeSize, scl = sumCoverListSize * 4;
   int sels = sumEndListSize;
   int intc = iCount[ts - 1];
   int bpolys = (bSize + 1) * 2, cpolys = (cSize + 1) * 2;

   int nid, clid;
   // int interCount=0;
   // int cs=coverListSizes[nid], ce=coverListSizes[nid+1];
   // int es=endListSizes[nid], ee=endListSizes[nid]+endListCounts[nid];

   // omp_set_dynamic(0);
   // omp_set_num_threads(20);

   // check thread times

   // #pragma omp parallel for private(clid)
   #pragma omp parallel for private(clid) schedule(dynamic, 100)
   // #pragma omp parallel for private(clid) schedule(auto)
   // #pragma omp parallel for private(clid) schedule(guided)
   for (nid = 1; nid < ts; nid++)
   {
      int interCount = 0;
      int einterCount = iCount[nid] - iCount[nid - 1] - 1;
      for (clid = coverListSizes[nid]; clid < coverListSizes[nid + 1]; ++clid)
      {
         // #pragma omp parallel
         // #pragma omp parallel sections num_threads(2)
         {
            int clid2 = clid * 4;
            // b \in C_B(v) and o \in C_O(v)
            // #pragma omp section
            // #pragma omp for nowait
            for (int clid3 = clid + 1; clid3 < coverListSizes[nid + 1]; ++clid3)
            {
               int bId, cId;
               int clid4 = clid3 * 4;
               bool proceed = false;
               // check if line segments are from different polygons
               // if(coverArray[clid2+1]!=coverArray[clid4+1]){
               // If clid2 is base polygon. Else clid2 is subject polygon
               if (coverArray[clid2 + 1] == bType && coverArray[clid4 + 1] == cType)
               {
                  bId = coverArray[clid2];
                  cId = coverArray[clid4];
                  proceed = true;
               }
               else if (coverArray[clid2 + 1] == cType && coverArray[clid4 + 1] == bType)
               {
                  bId = coverArray[clid4];
                  cId = coverArray[clid2];
                  proceed = true;
               }
               if (proceed)
               {
                  // cout<<"case1 bid: "<<bId/2<<" cid: "<<cId/2<<endl;
                  double alpha, beta;
                  if (LSMF(bPoly[bId], bPoly[bId + 1],
                           bPoly[bId + 2], bPoly[bId + 3],
                           cPoly[cId], cPoly[cId + 1],
                           cPoly[cId + 2], cPoly[cId + 3]))
                  {
                     int intersecType = lineSegmentIntersection(bPoly[bId], bPoly[bId + 1],
                                                                bPoly[bId + 2], bPoly[bId + 3],
                                                                cPoly[cId], cPoly[cId + 1],
                                                                cPoly[cId + 2], cPoly[cId + 3], &alpha, &beta);
                     // int intersecType=1;
                     if (intersecType)
                     {
                        // #pragma omp critical
                        {
                           bPolLineIds[iCount[nid - 1] + interCount] = bId;
                           cPolLineIds[iCount[nid - 1] + interCount] = cId;
                           intersectTypes[iCount[nid - 1] + interCount] = intersecType;
                           interCount += 1;
                        }
                     }
                  }
               }
            }
            // b \in C_B(v) and o \in E_O(v)
            // #pragma omp section
            // #pragma omp for nowait
            for (int elid = endListSizes[nid]; elid < endListSizes[nid] + endListCounts[nid]; ++elid)
            {
               int bId, cId;
               int elid2 = elid * 4;
               bool proceed = false;
               // check if line segments are from different polygons
               // if(coverArray[clid2+1]!=treeNodeEndListArray[elid].type){
               // If clid2 is base polygon. Else clid2 is subject polygon
               if (coverArray[clid2 + 1] == bType && treeNodeEndListArray[elid].type == cType)
               {
                  bId = coverArray[clid2];
                  cId = treeNodeEndListArray[elid].id;
                  proceed = true;
               }
               else if (coverArray[clid2 + 1] == cType && treeNodeEndListArray[elid].type == bType)
               {
                  cId = coverArray[clid2];
                  bId = treeNodeEndListArray[elid].id;
                  proceed = true;
               }
               if (proceed)
               {
                  double alpha, beta;
                  // cout<<nid<<" case2 bid: "<<bId/2<<" cid: "<<cId/2<<endl;
                  if (LSMF(bPoly[bId], bPoly[bId + 1],
                           bPoly[bId + 2], bPoly[bId + 3],
                           cPoly[cId], cPoly[cId + 1],
                           cPoly[cId + 2], cPoly[cId + 3]))
                  {
                     int intersecType = lineSegmentIntersection(bPoly[bId], bPoly[bId + 1],
                                                                bPoly[bId + 2], bPoly[bId + 3],
                                                                cPoly[cId], cPoly[cId + 1],
                                                                cPoly[cId + 2], cPoly[cId + 3], &alpha, &beta);
                     // int intersecType=1;
                     if (intersecType)
                     {
                        // #pragma omp critical
                        {
                           bPolLineIds[iCount[nid - 1] + einterCount] = bId;
                           cPolLineIds[iCount[nid - 1] + einterCount] = cId;
                           intersectTypes[iCount[nid - 1] + einterCount] = intersecType;
                           einterCount -= 1;
                        }
                        // cout<<"*********-----------**** nid="<<nid<<" bid="<<bId<<" cId="<<cId<<endl;
                     }
                  }
               }
            }
         }
      }
   }
}

/*
return alpha value after converting to int
*/
int alphaToInt(double alpha)
{
   if (alpha < EPSILON)
      return 0;
   return (int)pow(10, EPSILON_POSITIONS) * alpha;
}

/*
save only intersections in new arrays
*/
void SegTree ::saveIntersectionsParallel(
    REAL *bPoly, REAL *cPoly,
    int *bPolyLineIds, int *cPolyLineIds, int polySize,
    REAL *intersections, int *intersectTypes,
    int *intersectAlphaValuesB, int *intersectAlphaValuesC)
{
   int id;
   double alpha, beta;

   #pragma omp parallel for private(id, alpha, beta)
   for (id = 0; id < polySize; id++)
   {
      point b1, b2, c1, I;
      // int count2=0, nonDegenCount=0;

      int bId = bPolyLineIds[id];
      b1.x = bPoly[bId];
      b1.y = bPoly[bId + 1];
      b2.x = bPoly[bId + 2];
      b2.y = bPoly[bId + 3];

      int cId = cPolyLineIds[id];
      c1.x = cPoly[cId];
      c1.y = cPoly[cId + 1];
      // c2.x=*(*cPoly + (cId+2));
      // c2.y=*(*cPoly + (cId+3));

      int intersecType = lineSegmentIntersection(bPoly[bId], bPoly[bId + 1],
                                                 bPoly[bId + 2], bPoly[bId + 3],
                                                 cPoly[cId], cPoly[cId + 1],
                                                 cPoly[cId + 2], cPoly[cId + 3], &alpha, &beta);

      // cout<<"intersection:\nB(("<<bPoly[bId]<<", "<<bPoly[bId+1]<<"), ("<<bPoly[bId+2]<<", "<<bPoly[bId+3]<<")) X C(("<<cPoly[cId]<<", "<<cPoly[cId+1]<<"), ("<<cPoly[cId+2]<<", "<<cPoly[cId+3]<<"))\nalpha="<<alpha<<" beta="<<beta<<endl;
      // cout<<"testC E(("<<*(*cPoly + (cId))<<", "<<*(*cPoly+ (cId + 1))<<"), ("
      //    <<*(*cPoly + (cId + 2))<<", "<<*(*cPoly + (cId + 3))<<")) ";

      // cout<<intersecType<<" types ** "<<intersectTypes[bid]<<endl;
      if (intersecType)
      {
         // if((intersecType==1 || intersecType==3 || intersecType==5 || intersecType==7)){
         //    nonDegenCount++;
         //    count2=nonDegenCount;
         // }
         // else if((intersecType==2 || intersecType==4 || intersecType==6 || intersecType==8)){
         //    count2=0;
         // }

         switch (intersecType)
         {
         // case X_INTERSECTION:
         // I and I
         case 1:
            I = add(mulScalar((1.0 - alpha), b1), mulScalar(alpha, b2));
            intersections[id * 2] = I.x; // consider edge for the intersection array
            intersections[(id) * 2 + 1] = I.y;
            // intersectionsB[id*2]=I.x;
            // intersectionsB[id*2+1]=I.y;
            // intersectionsC[id*2]=I.x;
            // intersectionsC[id*2+1]=I.y;
            break;
         case 3:
         case 4:
         case 5:
         case 7:
         case 8:
            intersections[id * 2] = c1.x;
            intersections[(id) * 2 + 1] = c1.y;
            // intersectionsB[id*2]=c1.x;
            // intersectionsB[id*2+1]=c1.y;
            break;
         case 2:
         case 6:
            intersections[id * 2] = b1.x;
            intersections[(id) * 2 + 1] = b1.y;
            // intersectionsC[id*2]=b1.x;
            // intersectionsC[id*2+1]=b1.y;
            break;
         }
         intersectTypes[id] = intersecType;
         intersectAlphaValuesB[id] = alphaToInt(alpha);
         intersectAlphaValuesC[id] = alphaToInt(beta);
         // intersectAlphaValuesC[id]=(int)pow(10, EPSILON_POSITIONS)*beta;
      }
   }
}

void SegTree ::saveIntersections(
    REAL *bPoly, REAL *cPoly,
    int *bPolyLineIds, int *cPolyLineIds, int polySize,
    REAL *intersections, int *intersectTypes,
    int *intersectAlphaValuesB, int *intersectAlphaValuesC)
{
   int id;
   double alpha, beta;

   for (id = 0; id < polySize; id++)
   {
      point b1, b2, c1, I;
      // int count2=0, nonDegenCount=0;

      int bId = bPolyLineIds[id];
      b1.x = bPoly[bId];
      b1.y = bPoly[bId + 1];
      b2.x = bPoly[bId + 2];
      b2.y = bPoly[bId + 3];

      int cId = cPolyLineIds[id];
      c1.x = cPoly[cId];
      c1.y = cPoly[cId + 1];
      // c2.x=*(*cPoly + (cId+2));
      // c2.y=*(*cPoly + (cId+3));

      int intersecType = lineSegmentIntersection(bPoly[bId], bPoly[bId + 1],
                                                 bPoly[bId + 2], bPoly[bId + 3],
                                                 cPoly[cId], cPoly[cId + 1],
                                                 cPoly[cId + 2], cPoly[cId + 3], &alpha, &beta);

      // cout<<"testB E(("<<*(*bPoly + (bId))<<", "<<*(*bPoly+ (bId + 1))<<"), ("
      //    <<*(*bPoly + (bId + 2))<<", "<<*(*bPoly + (bId + 3))<<")) ";
      // cout<<"testC E(("<<*(*cPoly + (cId))<<", "<<*(*cPoly+ (cId + 1))<<"), ("
      //    <<*(*cPoly + (cId + 2))<<", "<<*(*cPoly + (cId + 3))<<")) ";

      // cout<<intersecType<<" types ** "<<intersectTypes[bid]<<endl;
      if (intersecType)
      {
         // if((intersecType==1 || intersecType==3 || intersecType==5 || intersecType==7)){
         //    nonDegenCount++;
         //    count2=nonDegenCount;
         // }
         // else if((intersecType==2 || intersecType==4 || intersecType==6 || intersecType==8)){
         //    count2=0;
         // }

         switch (intersecType)
         {
         // case X_INTERSECTION:
         // I and I
         case 1:
            I = add(mulScalar((1.0 - alpha), b1), mulScalar(alpha, b2));
            intersections[id * 2] = I.x; // consider edge for the intersection array
            intersections[(id) * 2 + 1] = I.y;
            // intersectionsB[id*2]=I.x;
            // intersectionsB[id*2+1]=I.y;
            // intersectionsC[id*2]=I.x;
            // intersectionsC[id*2+1]=I.y;
            break;
         case 3:
         case 4:
         case 5:
         case 7:
         case 8:
            intersections[id * 2] = c1.x;
            intersections[(id) * 2 + 1] = c1.y;
            // intersectionsB[id*2]=c1.x;
            // intersectionsB[id*2+1]=c1.y;
            break;
         case 2:
         case 6:
            intersections[id * 2] = b1.x;
            intersections[(id) * 2 + 1] = b1.y;
            // intersectionsC[id*2]=b1.x;
            // intersectionsC[id*2+1]=b1.y;
            break;
         }
         intersectTypes[id] = intersecType;
         intersectAlphaValuesB[id] = alphaToInt(alpha);
         intersectAlphaValuesC[id] = alphaToInt(beta);
         // intersectAlphaValuesC[id]=(int)pow(10, EPSILON_POSITIONS)*beta;
      }
   }
}

/*
sort intersections at each edge using alpha values
edgeCountPS: prefix sum of counts of edges
uniqueIntersectCount: Uniqueue edges in the intersections
sortedIndicies: sorted index order by the polygon Ids
sortedAlphaIndicies: Sortted index order of sortedIndicies by alpha values
*/
void SegTree::sortIntersectionsAtEdgesParallel(int *intersectAlphaValues, int *sortedIndicies, int cStart,
                                               int *allIntersectCounts, int *nonDegenIntersectCounts, int *polLineUniqueIds, int uniqueIntersectCount, int intersectCount,
                                               int *sortedAlphaIndicies)
{

   int eid, id;
   LONG *intersectAlphaValuesCp;
   intersectAlphaValuesCp = (LONG *)malloc((intersectCount) * sizeof(LONG));

#pragma omp parallel for
   for (int i = 0; i < intersectCount; ++i)
   {
      sortedAlphaIndicies[i] = i;
      intersectAlphaValuesCp[i] = intersectAlphaValues[sortedIndicies[i]];
      // cout<<"i="<<i<<" "<<intersectAlphaValuesCp[i]<<", ";
   }
   // for(int i=0; i<intersectCount; ++i){
   //    intersectAlphaValues[i]=intersectAlphaValuesCp[i];
   //    // cout<<"i="<<i<<" "<<intersectAlphaValuesCp[i]<<", ";
   // }
   // cout<<"copied "<<endl;
   // #pragma omp parallel for private(eid, id)
   for (id = 0; id < uniqueIntersectCount; ++id)
   {
      eid=(polLineUniqueIds[id]/2)-cStart;
      int s = allIntersectCounts[eid];
      int e = allIntersectCounts[eid + 1];
      int ns = nonDegenIntersectCounts[eid];
      int ne = nonDegenIntersectCounts[eid + 1];
      // cout<<id<<" "<<eid<<" sort2 start="<<s<<" end="<<e<<" ns="<<ns<<" ne="<<ne<<endl;
      if ((ne - ns) > 1)
         // if((e-s)>1)
         radixsort(intersectAlphaValuesCp, sortedAlphaIndicies, s, e, intersectCount, EPSILON_POSITIONS);
   }
   // cout<<"Done"<<endl;
   free(intersectAlphaValuesCp);
   // cout<<"Done"<<endl;
}

void SegTree::sortIntersectionsAtEdges(int *intersectAlphaValues, int *sortedIndicies, int cStart,
                                       int *allIntersectCounts, int *nonDegenIntersectCounts, int *polLineUniqueIds, int uniqueIntersectCount, int intersectCount,
                                       int *sortedAlphaIndicies)
{

   int eid, id;
   LONG *intersectAlphaValuesCp;
   intersectAlphaValuesCp = (LONG *)malloc((intersectCount) * sizeof(LONG));

   for (int i = 0; i < intersectCount; ++i)
   {
      sortedAlphaIndicies[i] = i;
      intersectAlphaValuesCp[i] = intersectAlphaValues[sortedIndicies[i]];
   }

   for (id = 0; id < uniqueIntersectCount; ++id)
   {
      eid = (polLineUniqueIds[id]/2)-cStart;
      int s = allIntersectCounts[eid];
      int e = allIntersectCounts[eid + 1];
      int ns = nonDegenIntersectCounts[eid];
      int ne = nonDegenIntersectCounts[eid + 1];
      if ((ne - ns) > 1)
         radixsort(intersectAlphaValuesCp, sortedAlphaIndicies, s, e, intersectCount, EPSILON_POSITIONS);
   }
   free(intersectAlphaValuesCp);
}

/*
copy intersections into source polygons at correct positions
Save neighbor IDs into arrays
can not to handle duplicate occurrances of intersections.
*/
void SegTree ::copyIntersectionsParallel(
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
    int *alphaValuesB, int *alphaValuesC)
{
   // int size=bSize*2;
   int *bNeighborIntersect, *cNeighborIntersect;
   bNeighborIntersect = (int *)malloc((polySize) * sizeof(int));
   cNeighborIntersect = (int *)malloc((polySize) * sizeof(int));

   int id;

// copy source vertices
// assume base is the largest polygon
#pragma omp parallel for private(id)
   for (id = 0; id < bSize; id++)
   {
      int id2 = id * 2;
      int pId = id + bNonDegenIntersectCountPS[id];
      int pId2 = pId * 2;

      intersectedB[pId2] = bPoly[id2];
      intersectedB[pId2 + 1] = bPoly[id2 + 1];
      alphaValuesB[pId] = -100; // default value to indentify non intersecting edges
      neighborsB[pId] = -100;   // default value to indentify non intersecting edges
      if (id < cSize)
      {
         int cid2 = id2 + cStart * 2;
         // int cid2 = id2;
         pId = id + cNonDegenIntersectCountPS[id];
         pId2 = 2 * pId;
         intersectedC[pId2] = cPoly[cid2];
         intersectedC[pId2 + 1] = cPoly[cid2 + 1];
         alphaValuesC[pId] = -100; // default value to indentify non intersecting edges
         neighborsC[pId] = -100;   // default value to indentify non intersecting edges
      }
   }

   // copy intersections into intersected array at the correct position
   // run time O(k') time, where k= # unique edge ids of intersections
   for (int id = 0; id < polySize; ++id)
   {
      int i, bid, cid, bId, nBId, nBId2, cId, nCId, nCId2, blid, blid2, prevblid, clid, clid2, prevclid;
      // save only non-degenerate cases: base polygon
      blid = bPolLineIdsSortedIndex[sortedIndiciesB[id]];
      clid = cPolLineIdsSortedIndex[sortedIndiciesC[id]];
      bid = id;
      cid = id;

      // if previous edge id and current edge id is same. It has been taken care of before
      // id=0 is saved without checkeing previous
      if (id > 0)
      {
         prevblid = bPolLineIdsSortedIndex[sortedIndiciesB[id - 1]];
         prevclid = cPolLineIdsSortedIndex[sortedIndiciesC[id - 1]];
      }
      if (id <= 0 || bPolyLineIds[prevblid] != bPolyLineIds[blid])
      {

         // cout<<"****dup "<<bPolyLineIds[blid]<<" "<<cPolyLineIds[blid]<<endl;
         // handles degenerate case at starting vertex of the edge
         if (intersectTypes[blid] == 2 || intersectTypes[blid] == 4 || intersectTypes[blid] == 6 || intersectTypes[blid] == 8)
         {
            // cout<<"**source bId="<<bPolyLineIds[blid]<<" cId="<<cPolyLineIds[blid]<<endl;
            bId = bPolyLineIds[blid] / 2;
            bId += bNonDegenIntersectCountPS[bId];
            alphaValuesB[bId] = intersectAlphaValuesB[blid];
            bNeighborIntersect[blid] = bId;
            if (bid < polySize - 1)
            {
               // update clid
               blid = bPolLineIdsSortedIndex[sortedIndiciesB[++bid]];
            }
         }
         blid2 = 2 * blid;
         bId = bPolyLineIds[blid] / 2;
         int start = bNonDegenIntersectCountPS[bId];
         int end = bNonDegenIntersectCountPS[bId + 1];
         // handles non-degenerate cases at a given edge
         for (i = 0, nBId = start + bId + 1, nBId2 = 2 * nBId; i < (end - start); ++i, ++nBId, nBId2 += 2)
         {
            // cout<<"***blid2="<<blid2<<" id="<<id<<" bId="<<bPolyLineIds[blid]<<" cId="<<cPolyLineIds[blid]<<" nBId="<<nBId<<" nBId2="<<nBId2<<" P1("<<intersections[blid2]<<", "<<intersections[blid2+1]<<")"<<endl;

            intersectedB[nBId2] = intersections[blid2];
            intersectedB[nBId2 + 1] = intersections[blid2 + 1];
            alphaValuesB[nBId] = intersectAlphaValuesB[blid];
            bNeighborIntersect[blid] = nBId;

            if (bid < polySize - 1)
            {
               // update blid
               blid = bPolLineIdsSortedIndex[sortedIndiciesB[++bid]];
               blid2 = 2 * blid;
            }
         }
      }

      if (id <= 0 || cPolyLineIds[prevclid] != cPolyLineIds[clid])
      {
         // save only non-degenerate cases: clipped polygon
         // handles degenerate case. First vertex
         if (intersectTypes[clid] == 3 || intersectTypes[clid] == 4 || intersectTypes[clid] == 5 || intersectTypes[clid] == 7 || intersectTypes[clid] == 8)
         {
            // cout<<"==source bId="<<bPolyLineIds[clid]<<" cId="<<cPolyLineIds[clid]<<" cId="<<cId<<" clid="<<clid<<endl;
            cId=cPolyLineIds[clid]/2-cStart;
            // cId=cPolyLineIds[clid]/2;
            cId += cNonDegenIntersectCountPS[cId];
            alphaValuesC[cId] = intersectAlphaValuesC[clid];
            cNeighborIntersect[clid] = cId;
            if (cid < polySize - 1)
            {
               // update clid
               clid = cPolLineIdsSortedIndex[sortedIndiciesC[++cid]];
            }
         }
         // handles non-degenerate cases
         clid2 = 2 * clid;
         cId=cPolyLineIds[clid]/2-cStart;
         // cId=cPolyLineIds[clid]/2;
         int start = cNonDegenIntersectCountPS[cId];
         int end = cNonDegenIntersectCountPS[cId + 1];

         // cout<<"cId="<<cId<<" Before start="<<start<<" end="<<end<<endl;

         for(i=0, nCId=start+cId+1, nCId2=2*nCId; i<(end-start); i++, ++nCId, nCId2+=2)
         {
            // cout<<"==clid2="<<clid2<<" id="<<id<<" bId="<<bPolyLineIds[clid]<<" cId="<<cPolyLineIds[clid]<<" nCId2="<<nCId2<<" nCId="<<nCId<<" P1("<<intersections[clid2]<<", "<<intersections[clid2+1]<<") "<<clid<<endl;
            intersectedC[nCId2] = intersections[clid2];
            intersectedC[nCId2 + 1] = intersections[clid2 + 1];
            alphaValuesC[nCId] = intersectAlphaValuesC[clid];
            cNeighborIntersect[clid] = nCId;

            if (cid < polySize - 1)
            {
               // update clid
               clid = cPolLineIdsSortedIndex[sortedIndiciesC[++cid]];
               clid2 = 2 * clid;
            }
         }
      }
   }
   // copy neighbor ids
   for (int i = 0; i < polySize; ++i)
   {
      // cout<<"i="<<i<<" bNeighbor="<<bNeighborIntersect[i]<<" cneighbor="<<cNeighborIntersect[i]<<endl;
      neighborsB[bNeighborIntersect[i]]=cNeighborIntersect[i]+1;
      neighborsC[cNeighborIntersect[i]]=bNeighborIntersect[i]+1;
   }
}

void SegTree ::copyIntersections(
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
    int *alphaValuesB, int *alphaValuesC)
{
   // int size=bSize*2;
   int *bNeighborIntersect, *cNeighborIntersect;
   bNeighborIntersect = (int *)malloc((polySize) * sizeof(int));
   cNeighborIntersect = (int *)malloc((polySize) * sizeof(int));

   int id;
   // copy source vertices
   // assume base is the largest polygon
   for (id = 0; id < bSize; id++)
   {
      int id2 = id * 2;
      int pId = id + bNonDegenIntersectCountPS[id];
      int pId2 = pId * 2;

      intersectedB[pId2] = bPoly[id2];
      intersectedB[pId2 + 1] = bPoly[id2 + 1];
      alphaValuesB[pId] = -100; // default value to indentify non intersecting edges
      neighborsB[pId] = -100;   // default value to indentify non intersecting edges
      if (id < cSize)
      {
          int cid2 = id2 + cStart * 2;
         // int cid2 = id2;
         pId = id + cNonDegenIntersectCountPS[id];
         pId2 = 2 * pId;
         intersectedC[pId2] = cPoly[cid2];
         intersectedC[pId2 + 1] = cPoly[cid2 + 1];
         alphaValuesC[pId] = -100; // default value to indentify non intersecting edges
         neighborsC[pId] = -100;   // default value to indentify non intersecting edges
      }
   }

   // copy intersections into intersected array at the correct position
   // run time O(k') time, where k= # unique edge ids of intersections
   for (int id = 0; id < polySize; ++id)
   {
      int i, bid, cid, bId, nBId, nBId2, cId, nCId, nCId2, blid, blid2, prevblid, clid, clid2, prevclid;
      // save only non-degenerate cases: base polygon
      blid = bPolLineIdsSortedIndex[sortedIndiciesB[id]];
      clid = cPolLineIdsSortedIndex[sortedIndiciesC[id]];
      bid = id;
      cid = id;

      // if previous edge id and current edge id is same. It has been taken care of before
      // id=0 is saved without checkeing previous
      if (id > 0)
      {
         prevblid = bPolLineIdsSortedIndex[sortedIndiciesB[id - 1]];
         prevclid = cPolLineIdsSortedIndex[sortedIndiciesC[id - 1]];
      }
      if (id <= 0 || bPolyLineIds[prevblid] != bPolyLineIds[blid])
      {

         // cout<<"****dup "<<bPolyLineIds[blid]<<" "<<cPolyLineIds[blid]<<endl;
         // handles degenerate case at starting vertex of the edge
         if (intersectTypes[blid] == 2 || intersectTypes[blid] == 4 || intersectTypes[blid] == 6 || intersectTypes[blid] == 8)
         {
            // cout<<"**source bId="<<bPolyLineIds[blid]<<" cId="<<cPolyLineIds[blid]<<endl;
            bId = bPolyLineIds[blid] / 2;
            bId += bNonDegenIntersectCountPS[bId];
            alphaValuesB[bId] = intersectAlphaValuesB[blid];
            bNeighborIntersect[blid] = bId;
            if (bid < polySize - 1)
            {
               // update clid
               blid = bPolLineIdsSortedIndex[sortedIndiciesB[++bid]];
            }
         }
         blid2 = 2 * blid;
         bId = bPolyLineIds[blid] / 2;
         int start = bNonDegenIntersectCountPS[bId];
         int end = bNonDegenIntersectCountPS[bId + 1];
         // handles non-degenerate cases at a given edge
         for (i = 0, nBId = start + bId + 1, nBId2 = 2 * nBId; i < (end - start); ++i, ++nBId, nBId2 += 2)
         {
            // cout<<"***blid2="<<blid2<<" id="<<id<<" bId="<<bPolyLineIds[blid]<<" cId="<<cPolyLineIds[blid]<<" nBId="<<nBId<<" nBId2="<<nBId2<<" P1("<<intersections[blid2]<<", "<<intersections[blid2+1]<<")"<<endl;

            intersectedB[nBId2] = intersections[blid2];
            intersectedB[nBId2 + 1] = intersections[blid2 + 1];
            alphaValuesB[nBId] = intersectAlphaValuesB[blid];
            bNeighborIntersect[blid] = nBId;

            if (bid < polySize - 1)
            {
               // update blid
               blid = bPolLineIdsSortedIndex[sortedIndiciesB[++bid]];
               blid2 = 2 * blid;
            }
         }
      }

      if (id <= 0 || cPolyLineIds[prevclid] != cPolyLineIds[clid])
      {
         // save only non-degenerate cases: clipped polygon
         // handles degenerate case. First vertex
         if (intersectTypes[clid] == 3 || intersectTypes[clid] == 4 || intersectTypes[clid] == 5 || intersectTypes[clid] == 7 || intersectTypes[clid] == 8)
         {
            // cout<<"==source bId="<<bPolyLineIds[clid]<<" cId="<<cPolyLineIds[clid]<<endl;
            cId=cPolyLineIds[clid]/2-cStart;
            cId += cNonDegenIntersectCountPS[cId];
            alphaValuesC[cId] = intersectAlphaValuesC[clid];
            cNeighborIntersect[clid] = cId;
            if (cid < polySize - 1)
            {
               // update clid
               clid = cPolLineIdsSortedIndex[sortedIndiciesC[++cid]];
            }
         }
         // handles non-degenerate cases
         clid2 = 2 * clid;
         cId=cPolyLineIds[clid]/2-cStart;
         int start = cNonDegenIntersectCountPS[cId];
         int end = cNonDegenIntersectCountPS[cId + 1];

         for (i = 0, nCId = start + cId + 1, nCId2 = 2 * nCId; i < (end - start); i++, ++nCId, nCId2 += 2)
         {
            // cout<<"==clid2="<<clid2<<" id="<<id<<" bId="<<bPolyLineIds[clid]<<" cId="<<cPolyLineIds[clid]<<" nCId2="<<nCId2<<" P1("<<intersections[clid2]<<", "<<intersections[clid2+1]<<")"<<endl;
            intersectedC[nCId2] = intersections[clid2];
            intersectedC[nCId2 + 1] = intersections[clid2 + 1];
            alphaValuesC[nCId] = intersectAlphaValuesC[clid];
            cNeighborIntersect[clid] = nCId;

            if (cid < polySize - 1)
            {
               // update clid
               clid = cPolLineIdsSortedIndex[sortedIndiciesC[++cid]];
               clid2 = 2 * clid;
            }
         }
      }
   }
   // cout<<"here"<<endl;
   // copy neighbor ids
   for (int i = 0; i < polySize; ++i)
   {
      neighborsB[bNeighborIntersect[i]] = cNeighborIntersect[i] + 1;
      neighborsC[cNeighborIntersect[i]] = bNeighborIntersect[i] + 1;
   }
}

// ---------------------------------------------
/*
-----------------------------------------------------------------
Function to get circular id of a given id
Runs in GPU
Called from Device
-------------------------------------------------------------------
*/
int getCircularId(int id, int maxCount)
{
   if (maxCount == id)
      return 0;
   else if (id == -1)
      return maxCount - 1;
   else
      return id;
}

/*
-----------------------------------------------------------------
Function to get relative position type
Runs in GPU
Called from Device
0 -> LEFT,
1 -> RIGHT,
2 -> IS_P_m,
3 -> IS_P_p
-------------------------------------------------------------------
*/
int oracle(int pMNId, int pPNId, int qId, const point &Q, const point &P1, const point &P2, const point &P3)
{
   // is Q linked to P1 ?
   if (pMNId != -100 && pMNId == qId)
      return 2;
   // is Q linked to P2 ?
   else if (pPNId != -100 && pPNId == qId)
      return 3;
   // check relative position of Q with respect to chain (P1,P2,P3)
   double s1 = A(Q, P1, P2);
   double s2 = A(Q, P2, P3);
   double s3 = A(P1, P2, P3);
   if (s3 > 0)
   {
      // chain makes a left turn
      if (s1 > 0 && s2 > 0)
         return 0;
      else
         return 1;
   }
   else
   {
      // chain makes a right turn (or is straight)
      if (s1 < 0 && s2 < 0)
         return 1;
      else
         return 0;
   }
}
// ---------------------------------------------

/*
-----------------------------------------------------------------
Function to get initial classification label
Runs in GPU
Called from Device
Intersection Labels
0  NONE,
1  CROSSING,
2  BOUNCING,
3  LEFT_ON,
4  RIGHT_ON,
5  ON_ON,
6  ON_LEFT,
7  ON_RIGHT,
8  DELAYED_CROSSING,
9  DELAYED_BOUNCING
-------------------------------------------------------------------
*/
int getInitialLabel(int qMType, int qPType)
{
   // check non-overlapping cases
   if ((qMType == 0 && qPType == 1) || (qMType == 1 && qPType == 0))
   {
      return 1;
   }
   if ((qMType == 0 && qPType == 0) || (qMType == 1 && qPType == 1))
   {
      return 2;
   }
   // check overlapping cases
   if (((qPType == 3) && (qMType == 1)) || ((qMType == 3) && (qPType == 1)))
      return 3;
   if (((qPType == 3) && (qMType == 0)) || ((qMType == 3) && (qPType == 0)))
      return 4;
   if (((qPType == 3) && (qMType == 2)) || ((qMType == 3) && (qPType == 2)))
      return 5;
   if (((qMType == 2) && (qPType == 1)) || ((qPType == 2) && (qMType == 1)))
      return 6;
   if (((qMType == 2) && (qPType == 0)) || ((qPType == 2) && (qMType == 0)))
      return 7;
   else
      return -102;
}

/*
Initial labeling
*/
void SegTree::calculateInitLabelParallel(int bSize, int *bNonDegenIntersectCounts,
                                         REAL *intersectedB, REAL *intersectedC, int *alphaValuesB,
                                         int *neighborsB,
                                         int bIntersectedCount, int cIntersectedCount, int *bInitLabels, int *cInitLabels)
{

   int bid;

#pragma omp parallel for private(bid)
   for (bid = 0; bid < bSize; ++bid)
   {
      int start = bid + bNonDegenIntersectCounts[bid], end = bid + bNonDegenIntersectCounts[bid + 1] + 1;
      // cout<<"=== "<<start<<" "<<end<<endl;
      // int start=psP2[id], end=psP2[id+1];
      int tmpId, nId, pMNId, pPNId;
      point pM, pP, qM, qP, current;
      int qMType, qPType, tmpIniLabel;
      int i;
      for (i = start; i < end; i++)
      {
         bInitLabels[i] = 0;
         if (alphaValuesB[i] != -100)
         { // consider intersections only
            current.x = intersectedB[i * 2];
            current.y = intersectedB[i * 2 + 1];
            tmpId = getCircularId(i - 1, bIntersectedCount);
            // determine local configuration at this intersection vertex
            pM.x = intersectedB[tmpId * 2];     // P-, predecessor of I on P
            pM.y = intersectedB[tmpId * 2 + 1]; // P-, predecessor of I on P
            // if(intersectedB[tmpId*2+2]!=-100)
            if (alphaValuesB[tmpId] != -100)
               pMNId = neighborsB[tmpId] - 1; // get neighbor id of P_m vertex
            else
               pMNId = -100;

            tmpId = getCircularId(i + 1, bIntersectedCount);
            pP.x = intersectedB[tmpId * 2];     // P+, successor of I on P
            pP.y = intersectedB[tmpId * 2 + 1]; // P+, successor of I on P
            if (alphaValuesB[tmpId] != -100)
               pPNId = neighborsB[tmpId] - 1; // get neighbor id of P_p vertex
            else
               pPNId = -100;

            // nId=getNeighborIndex(i, neighborMapP, neighborQ);
            nId = neighborsB[i] - 1;
            // cout<<"/-/ "<<nId<<endl;
            tmpId = getCircularId(nId - 1, cIntersectedCount);
            qM.x = intersectedC[tmpId * 2];     // Q-, predecessor of I on Q
            qM.y = intersectedC[tmpId * 2 + 1]; // Q-, predecessor of I on Q
            qMType = oracle(pMNId, pPNId, tmpId, qM, pM, current, pP);

            tmpId = getCircularId(nId + 1, cIntersectedCount);
            qP.x = intersectedC[tmpId * 2];     // Q+, successor of I on P
            qP.y = intersectedC[tmpId * 2 + 1]; // Q+, successor of I on P
            qPType = oracle(pMNId, pPNId, tmpId, qP, pM, current, pP);

            tmpIniLabel = getInitialLabel(qMType, qPType);
            bInitLabels[i] = tmpIniLabel;
            cInitLabels[nId] = tmpIniLabel;
            // cout<<"// "<<i<<endl;
         }
      }
   }
}

void SegTree::calculateInitLabel(int bSize, int *bNonDegenIntersectCounts,
                                 REAL *intersectedB, REAL *intersectedC, int *alphaValuesB,
                                 int *neighborsB,
                                 int bIntersectedCount, int cIntersectedCount, int *bInitLabels, int *cInitLabels)
{

   int bid;

   for (bid = 0; bid < bSize; ++bid)
   {
      int start = bid + bNonDegenIntersectCounts[bid], end = bid + bNonDegenIntersectCounts[bid + 1] + 1;
      // cout<<"=== "<<start<<" "<<end<<endl;
      // int start=psP2[id], end=psP2[id+1];
      int tmpId, nId, pMNId, pPNId;
      point pM, pP, qM, qP, current;
      int qMType, qPType, tmpIniLabel;
      int i;
      for (i = start; i < end; i++)
      {
         bInitLabels[i] = -100;
         if (alphaValuesB[i] != -100)
         { // consider intersections only
            current.x = intersectedB[i * 2];
            current.y = intersectedB[i * 2 + 1];
            tmpId = getCircularId(i - 1, bIntersectedCount);
            // determine local configuration at this intersection vertex
            pM.x = intersectedB[tmpId * 2];     // P-, predecessor of I on P
            pM.y = intersectedB[tmpId * 2 + 1]; // P-, predecessor of I on P
            // if(intersectedB[tmpId*2+2]!=-100)
            if (alphaValuesB[tmpId] != -100)
               pMNId = neighborsB[tmpId] - 1; // get neighbor id of P_m vertex
            else
               pMNId = -100;

            tmpId = getCircularId(i + 1, bIntersectedCount);
            pP.x = intersectedB[tmpId * 2];     // P+, successor of I on P
            pP.y = intersectedB[tmpId * 2 + 1]; // P+, successor of I on P
            if (alphaValuesB[tmpId] != -100)
               pPNId = neighborsB[tmpId] - 1; // get neighbor id of P_p vertex
            else
               pPNId = -100;

            // nId=getNeighborIndex(i, neighborMapP, neighborQ);
            nId = neighborsB[i] - 1;
            // cout<<"/-/ "<<nId<<endl;
            tmpId = getCircularId(nId - 1, cIntersectedCount);
            qM.x = intersectedC[tmpId * 2];     // Q-, predecessor of I on Q
            qM.y = intersectedC[tmpId * 2 + 1]; // Q-, predecessor of I on Q
            qMType = oracle(pMNId, pPNId, tmpId, qM, pM, current, pP);

            tmpId = getCircularId(nId + 1, cIntersectedCount);
            qP.x = intersectedC[tmpId * 2];     // Q+, successor of I on P
            qP.y = intersectedC[tmpId * 2 + 1]; // Q+, successor of I on P
            qPType = oracle(pMNId, pPNId, tmpId, qP, pM, current, pP);

            tmpIniLabel = getInitialLabel(qMType, qPType);
            bInitLabels[i] = tmpIniLabel;
            cInitLabels[nId] = tmpIniLabel;
            // cout<<"// "<<i<<endl;
         }
      }
   }
}

/*
Converts end lists of the segment tree into a 2-D array which has nodes and their endlists copy it to GPU
*/
/*void SegTree ::copyNodeEndListToGPU()
{
   int maxEndListSize = treeNode->at(0)->endList->size();

   for (int nid = 1; nid < treeSize; ++nid)
   {
      endListSizes[nid] = treeNode->at(nid)->endList->size();
      if (maxEndListSize < treeNode->at(nid)->endList->size())
      {
         maxEndListSize = treeNode->at(nid)->endList->size();
      }
   }

   sumEndListSize=treeSize*maxEndListSize;
   sumEndListSize*=4;
   int sel = sumEndListSize;
   int maxRowSize = maxEndListSize * 4;
   endArray = new REAL[sumEndListSize];
   vector<Edge> *endList;

   for (int nid = 0; nid < treeSize; ++nid)
   {
      endList = treeNode->at(nid)->endList;
      for (int elid = 0, elid2 = 0; elid < endList->size(); ++elid, elid2 += 4)
      {
         *(endArray + (nid * maxRowSize + elid2)) = endList->at(elid).id;
         *(endArray + (nid * maxRowSize + elid2 + 1)) = endList->at(elid).type;
         *(endArray + (nid * maxRowSize + elid2 + 2)) = endList->at(elid).start;
         *(endArray + (nid * maxRowSize + elid2 + 3)) = endList->at(elid).end;
      }
   }
// #pragma acc enter data copyin(this, endListSizes [0:treeSize], endArray [0:sel])
}*/

/*
Converts endEdgeSet list of the segment tree into a 2-D array which has nodes and their endlists copy it to GPU
this method is not space efficient and does not wor for larger data sets
*/
/*void SegTree ::copyNodeEndEdgeSetToGPU()
{
   int maxEndListSize = treeNode->at(0)->endEdgeSet->size();

   // cout<<"size "<<endListSizes[0]<<endl;
   for (int nid = 1; nid < treeSize; ++nid)
   {
      endListSizes[nid] = treeNode->at(nid)->endEdgeSet->size();
      if (maxEndListSize < treeNode->at(nid)->endEdgeSet->size())
      {
         maxEndListSize = treeNode->at(nid)->endEdgeSet->size();
         // cout<<"max size "<<maxEndListSize<<" @"<<nid<<endl;
      }
   }

   sumEndListSize=treeSize*maxEndListSize;
   sumEndListSize*=4;
   int sel = sumEndListSize;
   int maxRowSize = maxEndListSize * 4;

   // cout<<treeSize<<" "<<maxEndListSize<<" sumEndListSize "<<sumEndListSize<<endl;
   endArray = new REAL[sumEndListSize];
   set<Edge, edge_compare> *endEdgeSet;

   int elid2 = 0;
   // cout<<"****** end Array *******"<<endl;
   for (int nid = 0; nid < treeSize; ++nid)
   {
      endEdgeSet = treeNode->at(nid)->endEdgeSet;
      elid2 = 0;
      // cout<<"nid="<<nid<<" size "<<treeNode->at(nid)->endEdgeSet->size()<<endl;
      // for (Edge e : *treeNode->at(nid)->endEdgeSet)
      for(auto s=treeNode->at(nid)->endEdgeSet->begin(); s!=treeNode->at(nid)->endEdgeSet->end(); s++)
      {
         // *(endArray + (nid * maxRowSize + elid2)) =  e.id;
         // *(endArray + (nid * maxRowSize + elid2 + 1)) = e.type;
         // *(endArray + (nid * maxRowSize + elid2 + 2)) = e.start;
         // *(endArray + (nid * maxRowSize + elid2 + 3)) = e.end;
         // if(nid>34945) cout<<"t "<<s->id<<endl;
         *(endArray + (nid * maxRowSize + elid2)) =  s->id;
         *(endArray + (nid * maxRowSize + elid2 + 1)) = s->type;
         *(endArray + (nid * maxRowSize + elid2 + 2)) = s->start;
         *(endArray + (nid * maxRowSize + elid2 + 3)) = s->end;
         // cout<<"nid="<<nid<<" e.id: "<<*(endArray + (nid * maxRowSize + elid2))<<endl;
         elid2 += 4;
      }
   }
#pragma acc enter data copyin(this, endListSizes [0:treeSize], endArray [0:sel])
}*/

/*
prefix of an array
*/
/*void SegTree :: prefixSumSizeArray(){
   int *pSum = new int[treeSize];
   int i, j, e, ts=treeSize;
   int numWorkers=10;
   int buf=ts/numWorkers;
   int *localSum = new int[numWorkers];

   pSum[0]=coverListSizes[0];

   #pragma acc data create(pSum[0:treeSize], localsum[0:numWorkers])

   #pragma acc parallel loop
   for(i=0; i<numWorkers; ++i){
      e=i*buf+buf;
      #pragma acc parallel loop seq
      for(j=e-buf; j<e; ++j){
         pSum[j] = pSum[j-1]+coverListSizes[j];
      }
      localSum[i]=pSum[e-1];
      #pragma acc parallel loop
      for(j=e-buf; j<e; ++j){
         pSum[j] = pSum[j]+localSum[j-1];
      }
   }


   for(i=1; i<ts; ++i){
      cout<<coverListSizes[i]<<", ";
   }
   cout<<endl;
   for(i=1; i<ts; ++i){
      cout<<pSum[i]<<", ";
   }
   cout<<endl;
}*/

/*
copy intersections into source polygons at correct positions
Save neighbor IDs into arrays
able to handle duplicate occurrances of intersections. But will not work with init labeling. DO NOT USE!
*/
/*
void SegTree ::copyIntersectionsParallel2(
      REAL **bPoly, REAL **cPoly, int bSize, int cSize,
      int *bPolyLineIds, int *cPolyLineIds, int polySize,
      int *bPolLineIdsSortedIndex, int *cPolLineIdsSortedIndex,
      int *sortedIndiciesB, int *sortedIndiciesC,
      REAL *intersections, int *intersectTypes,
      int *bAllIntersectCounts, int *cAllIntersectCounts,
      int *bNonDegenIntersectCountPS, int *cNonDegenIntersectCountPS,
      int *intersectAlphaValuesB, int *intersectAlphaValuesC,
      REAL *intersectedB, REAL *intersectedC,
      int *neighborsB, int *neighborsC,
      int *alphaValuesB, int *alphaValuesC){
   // int id, id2;
   int size=bSize*2;
   int bNeighborIntersect[polySize], cNeighborIntersect[polySize];

   // copy source vertices
   // assume base is the largest polygon
   for(int id=0, id2=0; id<bSize; id2+=2, id++){
      int pId=id+bNonDegenIntersectCountPS[id];
      int pId2=pId*2;
      // cout<<"pId "<<pId<<endl;
      intersectedB[pId2]=*(*bPoly+id2);
      intersectedB[pId2+1]=*(*bPoly+id2+1);
      alphaValuesB[pId]=-100;    //default value to indentify non intersecting edges
      neighborsB[pId]=-100;    //default value to indentify non intersecting edges
      if(id<cSize){
         pId=id+cNonDegenIntersectCountPS[id];
         pId2=2*pId;
         intersectedC[pId2]=*(*cPoly+id2);
         intersectedC[pId2+1]=*(*cPoly+id2+1);
         // cout<<"c "<<pId2<<" "<<id2<<" "<<*(*cPoly+id2)<<" "<<*(*cPoly+id2+1);
         // cout<<" "<<intersectedC[pId]<<" "<<intersectedC[pId+1]<<endl;
         alphaValuesC[pId]=-100;    //default value to indentify non intersecting edges
         neighborsC[pId]=-100;    //default value to indentify non intersecting edges
      }
   }

   // copy intersections into intersected array at the correct position
   // run time O(k') time, where k= # unique edge ids of intersections
   for(int id=0; id<polySize; ++id){
      int i, bId, nBId, nBId2, cId, nCId, nCId2, blid, blid2, prevblid, clid, clid2, prevclid;
      // save only non-degenerate cases: base polygon
      blid=bPolLineIdsSortedIndex[sortedIndiciesB[id]];
      clid=cPolLineIdsSortedIndex[sortedIndiciesC[id]];

      // if previous edge id and current edge id is same. It has been taken care of before
      // id=0 is saved without checkeing previous
      if(id>0){
         prevblid=bPolLineIdsSortedIndex[sortedIndiciesB[id-1]];
         prevclid=cPolLineIdsSortedIndex[sortedIndiciesC[id-1]];
      }
      if(id<=0 || bPolyLineIds[prevblid]!=bPolyLineIds[blid]){

         // cout<<"****dup "<<bPolyLineIds[blid]<<" "<<cPolyLineIds[blid]<<endl;
         // handles degenerate cases at starting vertex of the edge
         if(intersectTypes[blid]==2 || intersectTypes[blid]==4 || intersectTypes[blid]==6 || intersectTypes[blid]==8){
            cout<<"**source bId="<<bPolyLineIds[blid]<<" cId="<<cPolyLineIds[blid]<<endl;
            bId=bPolyLineIds[blid]/2;
            bId+=bNonDegenIntersectCountPS[bId];
            alphaValuesB[bId]=intersectAlphaValuesB[clid];
            bNeighborIntersect[blid]=bId;
            // cout<<"bId "<<bId<<endl;
         }
         blid2=2*blid;
         bId=bPolyLineIds[blid]/2;
         int start=bAllIntersectCounts[bId];
         int end=bAllIntersectCounts[bId+1];
         int nonStart=bNonDegenIntersectCountPS[bId];
         // handles non-degenerate cases at a given edge
         nBId=nonStart+bId+1;
         nBId2=2*nBId;
         for(i=id; i<id+(end-start); ++i){
            // nBId=bId+id+1;
            // nBId2=2*nBId;
            if(intersectTypes[blid]==1 || intersectTypes[blid]==3 || intersectTypes[blid]==5 || intersectTypes[blid]==7){
               if(i-id>0){
                  // cout<<"id="<<id<<" i="<<i<<" blid="<<blid<<endl;
                  prevblid=bPolLineIdsSortedIndex[sortedIndiciesB[i-1]];
                  // prevclid=cPolLineIdsSortedIndex[sortedIndiciesC[i-1]];
                  // cout<<"blid="<<blid<<" prevblid="<<prevblid<<endl;
               }

               if(i-id<=0 || (bPolyLineIds[prevblid]!=bPolyLineIds[blid] || cPolyLineIds[prevblid]!=cPolyLineIds[blid])){
                  cout<<"***blid2="<<blid2<<" id="<<id<<" bId="<<bPolyLineIds[blid]<<" cId="<<cPolyLineIds[blid]<<" nBId="<<nBId<<" nBId2="<<nBId2<<" P1("<<intersections[blid2]<<", "<<intersections[blid2+1]<<")"<<endl;

                  intersectedB[nBId2]=intersections[blid2];
                  intersectedB[nBId2+1]=intersections[blid2+1];
                  alphaValuesB[nBId]=intersectAlphaValuesB[blid];
                  bNeighborIntersect[blid]=nBId;
                  // cout<<"saved"<<endl;
               }else{
                  // cout<<"--- bdup "<<bPolyLineIds[blid]<<" "<<cPolyLineIds[blid]<<endl;
                  // bNeighborIntersect[blid]=bNeighborIntersect[prevblid];
                  // cNeighborIntersect[clid]=cNeighborIntersect[prevclid];
                  intersectAlphaValuesB[blid]=-404;    //to identify a duplicate and to ignore
                  bNeighborIntersect[blid]=-404;
                  alphaValuesB[nBId]=-404;    //to identify a duplicate and to ignore
                  neighborsB[nBId]=-404;
               }
               ++nBId;
               nBId2+=2;
            }else{
               // first vertex is already added as a non-duplicate
               if(i>id){
                  bNeighborIntersect[blid]=-404;
                  intersectAlphaValuesB[blid]=-404;    //to identify a duplicate and to ignore
                  alphaValuesB[nBId]=-404;    //to identify a duplicate and to ignore
                  neighborsB[nBId]=-404;
               }
            }
            if(i<polySize-1){
               // update blid
               blid=bPolLineIdsSortedIndex[sortedIndiciesB[i+1]];
               blid2=2*blid;
               // clid=cPolLineIdsSortedIndex[sortedIndiciesC[i+1]];
            }
         }
      }
      // cout<<">>original bId="<<bPolyLineIds[clid]<<" cId="<<cPolyLineIds[clid]<<endl;
      if(id<=0 || cPolyLineIds[prevclid]!=cPolyLineIds[clid]){
         // save only non-degenerate cases: clipped polygon

         // cout<<"==original bId="<<bPolyLineIds[clid]<<" cId="<<cPolyLineIds[clid]<<endl;
         // handles degenerate cases
         if(intersectTypes[clid]==3 || intersectTypes[clid]==4 || intersectTypes[clid]==5 || intersectTypes[clid]==7 || intersectTypes[clid]==8){
            cout<<"==source bId="<<bPolyLineIds[clid]<<" cId="<<cPolyLineIds[clid]<<endl;
            cId=cPolyLineIds[clid]/2;
            cId+=cNonDegenIntersectCountPS[cId];
            alphaValuesC[cId]=intersectAlphaValuesC[clid];
            cNeighborIntersect[clid]=cId;
            // cout<<clid<<" cId "<<cId<<endl;
         }
         // handles non-degenerate cases
         clid2=2*clid;
         cId=cPolyLineIds[clid]/2;
         int start=cAllIntersectCounts[cId];
         int end=cAllIntersectCounts[cId+1];
         int nonStart=cNonDegenIntersectCountPS[cId];
         // nCId=cId+id+1;
         // nCId2=2*nCId;
         // cout<<cId<<" c st "<<start<<" "<<end<<endl;

         nCId=nonStart+cId+1;
         nCId2=2*nCId;
         for(i=id; i<id+(end-start); i++){
            cout<<id<<" cId "<<cId<<" nCId "<<nCId<<" "<<nCId2<<endl;
            if(intersectTypes[clid]==1 || intersectTypes[clid]==2 || intersectTypes[clid]==6){
               if(i-id>0){
                  // cout<<"id="<<id<<" i="<<i<<" clid="<<clid<<endl;
                  // prevblid=bPolLineIdsSortedIndex[sortedIndiciesB[i-1]];
                  prevclid=cPolLineIdsSortedIndex[sortedIndiciesC[i-1]];
                  // cout<<"clid="<<clid<<" prevclid="<<prevclid<<endl;
               }
               // cout<<"original cId="<<cPolyLineIds[clid]<<" bId="<<bPolyLineIds[clid]<<endl;
               if(i-id<=0 || (bPolyLineIds[prevclid]!=bPolyLineIds[clid] || cPolyLineIds[prevclid]!=cPolyLineIds[clid])){
                  cout<<"==clid2="<<clid2<<" id="<<id<<" bId="<<bPolyLineIds[clid]<<" cId="<<cPolyLineIds[clid]<<" nCId2="<<nCId2<<" P1("<<intersections[clid2]<<", "<<intersections[clid2+1]<<")"<<endl;
                  intersectedC[nCId2]=intersections[clid2];
                  intersectedC[nCId2+1]=intersections[clid2+1];
                  alphaValuesC[nCId]=intersectAlphaValuesC[clid];
                  cNeighborIntersect[clid]=nCId;
               }else{
                  // cout<<"--- c dup "<<bPolyLineIds[clid]<<" "<<cPolyLineIds[clid]<<endl;
                  // bNeighborIntersect[blid]=bNeighborIntersect[prevblid];
                  // cNeighborIntersect[clid]=cNeighborIntersect[prevclid];
                  intersectAlphaValuesC[clid]=-404;    //to identify a duplicate and to ignore
                  cNeighborIntersect[clid]=-404;
                  alphaValuesC[nCId]=-404;    //to identify a duplicate and to ignore
                  neighborsC[nCId]=-404;
               }
               ++nCId;
               nCId2+=2;
            }else{
               // first vertex is already added as a non-duplicate
               if(i>id){
                  cNeighborIntersect[clid]=-404;
                  intersectAlphaValuesC[clid]=-404;    //to identify a duplicate and to ignore
                  alphaValuesC[nCId]=-404;
                  neighborsC[nCId]=-404;
               }
            }
            if(i<polySize-1){
               // update clid
               // blid=bPolLineIdsSortedIndex[sortedIndiciesB[i+1]];
               // cout<<i<<" here "<<sortedIndiciesC[i+1]<<endl;
               clid=cPolLineIdsSortedIndex[sortedIndiciesC[i+1]];
               clid2=2*clid;
            }
         }
      }
   }
   cout<<"here"<<endl;
   int dupCount=0;
   // copy neighbor ids
   for(int i=0; i<polySize; ++i){
      // cout<<i<<" neighbors b="<<bNeighborIntersect[i]<<" alphaB="<<intersectAlphaValuesB[i]<<" c="<<cNeighborIntersect[i]<<" alphaC="<<intersectAlphaValuesC[i]<<endl;
      // +1 acts as an offset
      if(intersectAlphaValuesB[i]!=-404){
         neighborsB[bNeighborIntersect[i]]=cNeighborIntersect[i]+1;
         neighborsC[cNeighborIntersect[i]]=bNeighborIntersect[i]+1;
      }else{
         dupCount++;
      }
      if(DEBUG && intersectAlphaValuesB[i]==-100 && intersectAlphaValuesC[i]!=-100){
         cout<<"ERROR: Base and clipping have different combinations from duplicates. Need to fix this!!!"<<endl;
      }

   }
   cout<<"Duplicate count: "<<dupCount<<endl;
}*/

/*
original. Errors
map<REAL, pair<Node *, Node *>> *SegTree ::preparePointsToNodeMap()
{
   if (m_elementaryPoints == NULL && false == m_elementaryPoints->empty())
   {
      cout << "Error: m_elementaryPoints in insertToLeafLevelEndList" << endl;
      return NULL;
   }

   // special handling for 1st point
   REAL pt0 = m_elementaryPoints->at(0);
   // cout<<"pt0 is "<<pt0<<endl;

   map<REAL, pair<Node *, Node *>> *ptToAdjNodeMap = new map<REAL, pair<Node *, Node *>>();

   int firstLeafIndex = treeSize / 2;
   Node *firstLeaf = treeNode->at(firstLeafIndex);
   pair<Node *, Node *> nodePair1 = make_pair(firstLeaf, firstLeaf);

   pair<REAL, pair<Node *, Node *>> kv1 = make_pair(pt0, nodePair1);
   ptToAdjNodeMap->insert(kv1);

   int lastLeafIndex = treeSize - 1;
   // cout << endl
      //   << " firstLeafIndex " << firstLeafIndex << " lastLeafIndex " << lastLeafIndex << " m_elementaryPoints->size() " << m_elementaryPoints->size() << endl;

   int leafIndex = firstLeafIndex;
   int i;
   for (i = 1; i < m_elementaryPoints->size() - 1; i++, leafIndex++)
   {
      // cout<<"m_elemetry points "<<m_elementaryPoints->at(i)<<endl;
      Node *leftLeaf = treeNode->at(leafIndex);
      Node *rightLeaf = treeNode->at(leafIndex + 1);

      pair<Node *, Node *> ithNodePair = make_pair(leftLeaf, rightLeaf);

      REAL ithPt = m_elementaryPoints->at(i);
      pair<REAL, pair<Node *, Node *>> kv = make_pair(ithPt, ithNodePair);
      ptToAdjNodeMap->insert(kv);
   }
   // special handling for the last point
   // cout<<endl<<"i = "<<i<<endl;

   Node *lastLeaf = treeNode->at(lastLeafIndex);
   pair<Node *, Node *> lastNodePair = make_pair(lastLeaf, lastLeaf);

   REAL lastPt = m_elementaryPoints->at(m_elementaryPoints->size() - 1);
   pair<REAL, pair<Node *, Node *>> lastKV = make_pair(lastPt, lastNodePair);
   ptToAdjNodeMap->insert(lastKV);

   return ptToAdjNodeMap;
}*/

// ---------------------------
/*
radix sort - V1. Sorts one array and other array
*/
/*void countSort(REAL *otherData,  REAL *outputOther, int *input, int start, int n, int exp){
   int output[n]; // output array
   int lsVals[n];  //save least significant digits of the input
   int i, count[10] = {0};

   // Store count of occurrences in count[]
   for(i=start; i<(start+n); i++){
      lsVals[i]=input[i]%exp;
      count[lsVals[i]]++;
   }

   // Change count[i] so that count[i] now contains actual
   //  position of this digit in output[]
   for(i=1; i<10; i++)
      count[i]+=count[i-1];

   // Build the output array
   for (i=(start+n-1); i>=start; i--) {
      output[count[lsVals[i]]-1]=input[i];
      outputOther[count[lsVals[i]]-1]=otherData[i];
      count[lsVals[i]]--;
   }

   // Copy the output array to arr[], so that arr[] now
   // contains sorted numbers according to current digit
   for(i=start; i<(start+n); i++){
      input[i]=output[i];
      otherData[i]=outputOther[i];
      // cout<<"i "<<i<<" "<<output[i]<<" "<<input[i]<<endl;
   }
}

// The main function to that sorts arr[] of size n using
// Radix Sort
void radixsort(REAL *otherData, int *data, int start, int n, int m){
   // Do counting sort for every digit. Note that instead
   // of passing digit number, exp is passed. exp is 10^i
   // where i is current digit number
   for (int exp=1, i=0; m>i; exp*=10){
      countSort(otherData, data, start, n, exp*10);
      // cout<<"round "<<exp<<"\n\n"<<endl;
   }
}*/
// ---------------------------

/*
check edge intersections between two given edges
returns intersection type
*/
/*int edgeIntersection
      (REAL l1p1x, REAL l1p1y, REAL l1p2x, REAL l1p2y,
      REAL l2p1x, REAL l2p1y, REAL l2p2x, REAL l2p2y,
      double &alpha, double &beta,
      int bPolyId, int cPolyId, int &bPolyIdIntersect, int &cPolyIdIntersect){
   point lines[4];
   int type;

   lines[0].x = l1p1x;
   lines[0].y = l1p1y;
   lines[1].x = l1p2x;
   lines[1].y = l1p2y;
   lines[2].x = l2p1x;
   lines[2].y = l2p1y;
   lines[3].x = l2p2x;
   lines[3].y = l2p2y;
   type=lineSegmentIntersection
               (lines[0].x, lines[0].y,  lines[1].x, lines[1].y,
               lines[2].x, lines[2].y, lines[3].x, lines[3].y,
               alpha, beta);
   if(type){
      bPolyIdIntersect=bPolyId;
      cPolyIdIntersect=cPolyId;
      return type;
   }
   // type=lineSegmentIntersection
   //             (lines[0].x, lines[0].y,  lines[1].x, lines[1].y,
   //             lines[3].x, lines[3].y, lines[2].x, lines[2].y,
   //             alpha, beta);
   // if(type){
   //    bPolyIdIntersect=bPolyId;
   //    cPolyIdIntersect=cPolyId+2;
   //    return type;
   // }
   // type=lineSegmentIntersection
   //             (lines[1].x, lines[1].y,  lines[0].x, lines[0].y,
   //             lines[2].x, lines[2].y, lines[3].x, lines[3].y,
   //             alpha, beta);
   // if(type){
   //    bPolyIdIntersect=bPolyId+2;
   //    cPolyIdIntersect=cPolyId;
   //    return type;
   // }type=lineSegmentIntersection
   //             (lines[1].x, lines[1].y,  lines[0].x, lines[0].y,
   //             lines[3].x, lines[3].y, lines[2].x, lines[2].y,
   //             alpha, beta);
   // if(type){
   //    bPolyIdIntersect=bPolyId+2;
   //    cPolyIdIntersect=cPolyId+2;
   //    return type;
   // }

   return type; // default exit when no intersections are reported

}*/

// -----------------------------------------------
// Query segment tree methods. Partially completed

/*vector<Edge> *SegTree ::querySegTree(REAL qx)
{
   // Input: the root of a (subtree of a) segment tree and query point qx.
   // Output: All intervals in the tree containing qx.
   vector<Edge> *edges = new vector<Edge>();
   // if query is out of root's range, it does not include in the tree
   if (qx < treeNode->at(1)->interval.start || qx > treeNode->at(1)->interval.end)
   {
      return edges;
   }

   recQueryTree(qx, treeNode->at(1), 1, edges);
   recQueryTreeIterative(qx, treeNode->at(1), 1, edges);

   return edges;
}

void SegTree ::recQueryTree(REAL qx, Node *root, int nodeIdx, vector<Edge> *edges)
{
   // edges->push_back(root.coverList);
   string edge = root->interval.toString();
   // cout<<"cover list size "<<root->coverList->size()<<" "<<nodeIdx<<" "<<edge<<endl;
   if (root->coverList != NULL && root->coverList->size() > 0)
   {
      // cout<<"cover list size "<<root.coverList->size()<<endl;
      edges->insert(edges->end(), root->coverList->begin(), root->coverList->end());
   }
   // if root is not a leaf
   int leftIdx = 2 * nodeIdx;

   if (leftIdx < treeSize)
   {
      Node *leftChild = treeNode->at(leftIdx);
      if (leftChild->interval.contains(qx))
      {
         recQueryTree(qx, leftChild, leftIdx, edges);
      }
   }

   int rightIdx = 2 * nodeIdx + 1;
   if (rightIdx < treeSize)
   {
      Node *rightChild = treeNode->at(rightIdx);
      if (rightChild->interval.contains(qx))
      {
         recQueryTree(qx, rightChild, rightIdx, edges);
      }
   }
}

void SegTree ::recQueryTreeIterative(REAL qx, Node *root, int nodeIdx, vector<Edge> *edges)
{
   int maxCoverListSize = 0;
   int coverListSizes[treeSize];
   // #pragma acc parallel loop copyin(treeNode[:treeSize])
   for (int nid = 1; nid < treeSize; ++nid)
   {
      coverListSizes[nid] = 0;
      if (treeNode->at(nid)->interval.contains(qx))
      {
         if (treeNode->at(nid)->coverList != NULL && treeNode->at(nid)->coverList->size() > 0)
         {
            // edges->insert(edges->end(), treeNode->at(nid)->coverList->begin(), treeNode->at(nid)->coverList->end());
            coverListSizes[nid] = treeNode->at(nid)->coverList->size();
            if (maxCoverListSize < treeNode->at(nid)->coverList->size())
            {
               maxCoverListSize = treeNode->at(nid)->coverList->size();
            }
         }
      }
   }

   int resultSize = treeSize * maxCoverListSize;
   Edge **qresults = new Edge *[resultSize];
   vector<Edge> *coverList;

   // #pragma acc parallel loop copyin(treeNode[0:treeSize]) present(coverListSizes[:treeSize]) copyout(qresults[:resultSize])
   for (int nid = 1; nid < treeSize; ++nid)
   {
      // cout<<"nn "<<nid<<"-> "<<treeNode->at(nid)->interval.start<<" "<<treeNode->at(nid)->interval.end<<endl;
      if (treeNode->at(nid)->interval.contains(qx))
      {
         if (treeNode->at(nid)->coverList != NULL && treeNode->at(nid)->coverList->size() > 0)
         {
            coverList = treeNode->at(nid)->coverList;
            for (int clid = 0; clid < coverList->size(); ++clid)
            {
               *(qresults + (nid * maxCoverListSize + clid)) = &coverList->at(clid);
            }
         }
      }
   }

   if (DEBUG_MODE)
      cout << "Max cover list size: " << maxCoverListSize << endl;
   // for(int i=1; i<treeSize; ++i){
   //    cout<<"Node="<<i<<" ["<<treeNode->at(i)->interval.start<<"-"<<treeNode->at(i)->interval.end <<"] cl_size="<<coverListSizes[i]<<endl;
   //    for(int j=0; j<coverListSizes[i]; ++j){
   //       cout<<"\tline_id> "<<(*(qresults+(i*maxCoverListSize+j)))->id<<endl;
   //    }
   //    cout<<endl;
   // }
   // cout<<"\ntotal "<<val<<endl;
}*/
// -----------------------------------------------

// original
/*
void SegTree ::insertToLeafLevelEdgeList(vector<Edge> *edges)
{
   map<REAL, pair<Node *, Node *>> *ptToNodeMap = preparePointsToNodeMap();

   if (ptToNodeMap == NULL)
   {
      cout << "BUG: map empty" << endl;
      exit(1);
   }
   else
   {
      cout << "Map Size " << ptToNodeMap->size() << endl;
   }

   std::map<REAL, pair<Node *, Node *>>::iterator it;
   // cout<<"before for loop all edges "<<endl;
   for (Edge e : *edges)
   {
      // start point
      // if(e.type ==1) cout<<"Start point "<<e.start<<" "<<e.end<<" type: "<<e.type<<" id: "<<e.id<<endl;
      it = ptToNodeMap->find(e.start);
      if (it != ptToNodeMap->end())
      // if (e.start!=e.end && it != ptToNodeMap->end())
      {
         pair<Node *, Node *> nodePair = ptToNodeMap->at(e.start);

         Node *after = nodePair.second;
         after->addToEndList(e);
         after->addToEndEdgeSet(e);

         // if(after->interval.start==e.start || after->interval.start==e.end || after->interval.end==e.start || after->interval.end==e.end)
         // {
         // // if(after->interval.start==e.first || after->interval.end==e.first){
         //    Node *before = nodePair.first;
         //    before->addToEndList(e);
         //    before->addToEndEdgeSet(e);
         // }
         // if(e.type ==1) cout<<"Added point "<<e.start<<" "<<e.end<<" type: "<<e.type<<" id: "<<e.id<<endl;
      }

      // end point
      it = ptToNodeMap->find(e.end);
      if (it != ptToNodeMap->end())
      // if (e.start!=e.end && it != ptToNodeMap->end())
      {
         pair<Node *, Node *> nodePair = ptToNodeMap->at(e.end);
         Node *before = nodePair.first;

         before->addToEndList(e);
         before->addToEndEdgeSet(e);

         // if(before->interval.start==e.start || before->interval.start==e.end || before->interval.end==e.start || before->interval.end==e.end)
         // {
         //    Node *after = nodePair.second;
         //    after->addToEndList(e);
         //    after->addToEndEdgeSet(e);
         // }
      }
   }

   // convert set ds to vector ds
   // int firstLeafIndex = treeSize/2;
   // int lastLeafIndex = treeSize-1;

   // for(int nodeIdx = firstLeafIndex; nodeIdx <= lastLeafIndex; nodeIdx++)
   // {
   //    Node *leaf = treeNode->at(nodeIdx);
   //    //leaf->endEdgeSet
   //    set<Edge,edge_compare> *endSet = leaf->endEdgeSet;
   //    leaf->endList = new vector<Edge>( endSet->begin(), endSet->end());
   // }

}*/

// Buddhi
/*
void SegTree ::insertAllToEdgeList2(vector<Edge> *edges)
{
   for (Edge e : *edges){
      insertToEdgeSet(e);
   }
}

void SegTree ::insertEdgeToEdegSet(Edge x, Node *root, int nodeIdx)
{
   if (root->interval.contains(x.start) || root->interval.contains(x.end)){
      root->addToEndEdgeSet(x);
   }
   else{
      int leftIdx = 2 * nodeIdx;
      if (leftIdx < treeSize){
         Node *leftChild = treeNode->at(leftIdx);
         if (x.intersects(leftChild->interval))
         {
            insertEdgeToEdegSet(x, leftChild, leftIdx);
         }
      }
      int rightIdx = 2 * nodeIdx + 1;
      if (rightIdx < treeSize){
         Node *rightChild = treeNode->at(rightIdx);
         if (x.intersects(rightChild->interval))
         {
            insertEdgeToEdegSet(x, rightChild, rightIdx);
         }
      }
   }
}

void SegTree ::insertToEdgeSet(Edge x)
{
   insertEdgeToEdegSet(x, treeNode->at(1), 1);
}*/