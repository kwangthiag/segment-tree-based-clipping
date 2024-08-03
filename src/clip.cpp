/*
clipping two given polygons
Step1: Convert polygons into segment tree DS
Step2: for each scan beam, find lines segments intersect
    Stab query to find line segments pass through an imaginary line
    Step3: Label line segments as L, R
    Step4: Label line segments MX, MN
Step5: Merge scan beams pair-wise
    Step6: Use sorted order of x coordinates and MX, MN, L, R labels
*/

// #include <stdio.h>
#include <iomanip>
#include <iostream>
#include <chrono>

#include <thrust/scan.h>
#include <thrust/execution_policy.h>
// #include <thrust/device_malloc.h>
// #include <thrust/device_ptr.h>
// #include <thrust/copy.h>

#ifdef _OPENACC
#include <openacc.h>
#endif

#include "../p_seg_tree/segment_tree/segment_tree.h"
#include "../p_seg_tree/util/thrust_func.h"

// #include "lib/optimizedFostersAlgorithm/polyclip_time.h"
string outputFile = string("../results/Fig-def-R.txt");
string outputFile1 = string("../results/Fig-def-R");
#include "../foster/polyclip.h"

#define MALLOC(p, b, s, t)                              \
  {                                                     \
    if ((b) > 0)                                        \
    {                                                   \
      p = (t *)malloc(b);                               \
      if (!(p))                                         \
      {                                                 \
        fprintf(stderr, "gpc malloc failure: %s\n", s); \
        exit(0);                                        \
      }                                                 \
    }                                                   \
    else                                                \
      p = NULL;                                         \
  }

using namespace std::chrono;

void gpc_read_polygon(FILE *bfp, FILE *cfp, REAL **bPoly, REAL **cPoly, int &bSize, int &cSize, Edge **earr, int &intervalCount, REAL &minInterval, REAL &maxInterval)
{
  int v, v2, id;
  int tmpSize;
  FILE *tmpFp;
  Edge curr;
  point2D vtex;
  polygon bPolyList, cPolyList;

  fscanf(bfp, "%d", &bSize);
  fscanf(cfp, "%d", &cSize);
  // if subject polygon is larger, we assign it as the base polygon
  if (bSize < cSize)
  {
    tmpSize = bSize;
    tmpFp = bfp;

    bSize = cSize;
    bfp = cfp;

    cSize = tmpSize;
    cfp = tmpFp;
    #ifdef DG
    cout << "Base and clipping swapped" << endl;
    #endif
  }
  MALLOC(*bPoly, (bSize + 1) * 2 * sizeof(REAL), "Base polygon allocating", REAL);
  MALLOC(*cPoly, (cSize + 1) * 2 * sizeof(REAL), "Clipping polygon allocating", REAL);

  *earr = (Edge *)malloc(((bSize) + (cSize) + 2) * sizeof(Edge));
  // read first row befor loop. Required for x coordinate copying
  // use v2 instead v. Later one multiplication can be reduced using this way
  fscanf(bfp, "%lf,%lf", (*bPoly), (*bPoly + 1));
  fscanf(cfp, "%lf,%lf", (*cPoly), (*cPoly + 1));
  vtex = point2D(*(*bPoly), *(*bPoly + 1));
  bPolyList.newVertex(vtex, true);
  vtex = point2D(*(*cPoly), *(*cPoly + 1));
  cPolyList.newVertex(vtex, true);

  for (v = 1, v2 = 2, id = 0; v < bSize; v++, v2 += 2)
  {
    fscanf(bfp, "%lf,%lf", (*bPoly + v2), (*bPoly + v2 + 1));
    curr.getData(v2 - 2, 0, *(*bPoly + v2 - 2), *(*bPoly + v2));
    *(*earr + id) = curr;
    id++;

    if (v == 1)
    {
      minInterval = curr.start;
      maxInterval = curr.end;
    }

    if (minInterval > curr.start)
    {
      minInterval = curr.start;
    }
    if (maxInterval < curr.end)
    {
      maxInterval = curr.end;
    }

    vtex = point2D(*(*bPoly + v2), *(*bPoly + v2 + 1));
    bPolyList.newVertex(vtex, true);
    if (v < cSize)
    {
      fscanf(cfp, "%lf,%lf", (*cPoly + v2), (*cPoly + v2 + 1));
      curr.getData(v2 - 2, 1, *(*cPoly + v2 - 2), *(*cPoly + v2));
      *(*earr + id) = curr;
      id++;

      if (minInterval > curr.start)
      {
        minInterval = curr.start;
      }
      if (maxInterval < curr.end)
      {
        maxInterval = curr.end;
      }

      vtex = point2D(*(*cPoly + v2), *(*cPoly + v2 + 1));
      cPolyList.newVertex(vtex, true);
    }
  }
  fclose(bfp);
  fclose(cfp);

  // push line segments that closes the polygons
  *(*bPoly + bSize * 2) = *(*bPoly);
  *(*bPoly + bSize * 2 + 1) = *(*bPoly + 1);
  curr.getData(bSize * 2 - 2, 0, *(*bPoly + bSize * 2 - 2), *(*bPoly + bSize * 2));
  *(*earr + id) = curr;
  id++;

  *(*cPoly + cSize * 2) = *(*cPoly);
  *(*cPoly + cSize * 2 + 1) = *(*cPoly + 1);
  curr.getData(cSize * 2 - 2, 1, *(*cPoly + cSize * 2 - 2), *(*cPoly + cSize * 2));
  *(*earr + id) = curr;
  id++;

  bPolysList.push_back(bPolyList);
  cPolysList.push_back(cPolyList);
  intervalCount = id;

  if (intervalCount >= (bSize) + (cSize) + 2)
  {
    cout << "ERROR: interval overflow. Check gpc_read_polygon()" << endl;
  }

  #ifdef DG
    printf("Base polygon with %d vertices reading completed!\n", bSize++);
    printf("Clipping polygon with %d vertices reading completed!\n", cSize++);
    printf("%d intervals saved\n", intervalCount);
    printf("Min val: %f. Max val: %f.\n", minInterval, maxInterval);
  #endif
  
}

/*
muti polygon reading
*/
void gpc_read_multi_polygon(FILE *bfp, FILE **cfp, int cPolyCount, REAL **bPoly, REAL **cPoly, int &bSize, int &cSize, int *cSizes, Edge **earr, int &intervalCount, REAL &minInterval, REAL &maxInterval)
{
  int v, v2, id, c;
  int tmpSize;
  FILE *tmpFp;
  Edge curr;
  point2D vtex;
  polygon bPolyList;

  fscanf(bfp, "%d", &bSize);

  intervalCount = 0;
  cSize = 0;
  cSizes[0] = 0;
  for (int i = 0; i < cPolyCount; ++i)
  {
    fscanf(*(cfp + i), "%d", &c);
    cSizes[i + 1] = cSizes[i] + c + 1;
    cSize += c;

    #ifdef DG
    cout<<"clipping poly id="<<i<<" size="<<c<<" prefix sum at "<<i+1<<"="<<cSizes[i+1]<<endl;
    #endif
    
  }
  MALLOC(*bPoly, (bSize + 1) * 2 * sizeof(REAL), "Base polygon allocating", REAL);
  MALLOC(*cPoly, (cSize + cPolyCount) * 2 * sizeof(REAL), "Clipping polygon allocating", REAL);

  *earr = (Edge *)malloc(((bSize) + (cSize) + 1 + cPolyCount) * sizeof(Edge));
  // read first row befor loop. Required for x coordinate copying
  // use v2 instead v. Later one multiplication can be reduced using this way
  fscanf(bfp, "%lf,%lf", (*bPoly), (*bPoly + 1));
  vtex = point2D(*(*bPoly), *(*bPoly + 1));
  bPolyList.newVertex(vtex, true);

  for (v = 1, v2 = 2, id = 0; v < bSize; v++, v2 += 2)
  {
    fscanf(bfp, "%lf,%lf", (*bPoly + v2), (*bPoly + v2 + 1));
    curr.getData(v2 - 2, 0, *(*bPoly + v2 - 2), *(*bPoly + v2));
    *(*earr + id) = curr;
    id++;

    if (v == 1)
    {
      minInterval = curr.start;
      maxInterval = curr.end;
    }
    if (minInterval > curr.start)
    {
      minInterval = curr.start;
    }
    if (maxInterval < curr.end)
    {
      maxInterval = curr.end;
    }

    vtex = point2D(*(*bPoly + v2), *(*bPoly + v2 + 1));
    bPolyList.newVertex(vtex, true);
  }
  fclose(bfp);
  // push line segments that closes the polygons
  *(*bPoly + bSize * 2) = *(*bPoly);
  *(*bPoly + bSize * 2 + 1) = *(*bPoly + 1);
  curr.getData(bSize * 2 - 2, 0, *(*bPoly + bSize * 2 - 2), *(*bPoly + bSize * 2));
  *(*earr + id) = curr;
  id++;

  bPolysList.push_back(bPolyList);

  v2 = 0;
  REAL fx, fy;
  int lindex;
  for (int i = 0; i < cPolyCount; ++i)
  {
    polygon cPolyList;

    fscanf(*(cfp+i), "%lf,%lf", (*cPoly+v2), (*cPoly+v2+1));
    // cout<<v2<<" ==== "<<*(*cPoly+v2)<<" "<<*(*cPoly+v2+1)<<endl;

    fx = *(*cPoly+v2);
    fy = *(*cPoly+v2+1);
    vtex = point2D(*(*cPoly+v2), *(*cPoly+v2+1));
    cPolyList.newVertex(vtex, true);
    v2 += 2;

    // cout<<cSizes[i+1]-cSizes[i]<<" first "<<fx<<" "<<fy<<endl;

    for(v=1; v<(cSizes[i+1]-cSizes[i]-1); v++, v2+=2)
    {
      fscanf(*(cfp+i), "%lf,%lf", (*cPoly+v2), (*cPoly+v2+1));
      // cout<<v2<<" ==== "<<*(*cPoly+v2)<<" "<<*(*cPoly+v2+1)<<endl;
      curr.getData(v2-2, i+1, *(*cPoly+v2-2), *(*cPoly+v2));
      *(*earr + id) = curr;
      // cout<<(*earr+id)->id<<" -- "<<(*earr+id)->type<<" *** "<<(*earr+id)->start<<" "<<(*earr+id)->end<<endl;
      id++;
      if (minInterval > curr.start)
      {
        minInterval = curr.start;
      }
      if (maxInterval < curr.end)
      {
        maxInterval = curr.end;
      }
      vtex = point2D(*(*cPoly + v2), *(*cPoly + v2 + 1));
      cPolyList.newVertex(vtex, true);
    }
    // cout<<"b done "<<id<<endl;
    // push line segments that closes the polygons
    // *(*cPoly+v2)=fx;
    // *(*cPoly+v2+1)=fy;
    // curr.getData(v2-2, i+1, *(*cPoly+v2-2), *(*cPoly+v2));
    // lindex=cSizes[i+1]*2-2;
    // cout<<"== "<<lindex<<" "<<v2<<endl;
    lindex=v2;
    *(*cPoly+lindex)=fx;
    *(*cPoly+lindex+1)=fy;
    curr.getData(lindex-2, i+1, *(*cPoly+lindex-2), *(*cPoly+lindex));
    *(*earr+id)=curr;
    // cout<<(*earr+id)->id<<" -- "<<(*earr+id)->type<<" *** "<<(*earr+id)->start<<" "<<(*earr+id)->end<<endl;
    // cout<<v2<<" ==== "<<*(*cPoly+v2)<<" "<<*(*cPoly+v2+1)<<endl;
    
    id++;
    v2 += 2;
    // cout<<"b done"<<endl;
    cPolysList.push_back(cPolyList);
    fclose(*(cfp + i));
  }
  intervalCount = id;
  if (intervalCount >= (bSize) + (cSize) + 1 + cPolyCount)
  {
    cout << "ERROR: interval overflow. Check gpc_read_multi_polygon()" << endl;
  } 

  #ifdef DG
    printf("Base polygon size=%d vertices reading completed!\n", bSize);
    printf("%d clipping polygon(s) with %d vertices reading completed!\n", cPolyCount, cSize);
    printf("%d intervals saved\n", intervalCount);
    printf("Min val: %f. Max val: %f.\n", minInterval, maxInterval);
  #endif
}

/*
duplicate the base layer polygon allowing each clipping layer polygon to calculate intersection independently
*/
void duplicateBasePoly(int cPolyCount){
  for(int cid=0; cid<cPolyCount-1; ++cid){
    polygon bPolyList;
    vertex *curr;
    curr=bPolysList[0].root;
    for(int i=0; i<bPolysList[0].size; ++i){
      bPolyList.newVertex(curr->p, true);
      curr=curr->next;
    }
    bPolysList.push_back(bPolyList);
  }
}

/*
helper function to print an array
*/
void printDoubleArray(double *arr, int size, string msg)
{
  cout << "\n*****" << msg << ": " << size << endl;
  for (int i = 0; i < size; ++i)
  {
    cout << arr[i] << " ";
  }
  cout << endl;
}
void printIntArray(int *arr, int size, string msg)
{
  cout << "\n*****" << msg << ": " << size << endl;
  for (int i = 0; i < size; ++i)
  {
    cout << arr[i] << " ";
  }
  cout << endl;
}

/*
print using sorted index array
*/
void printIntIndexArray(int *arr, int *indexArr, int size, string msg)
{
  cout << "\n*******" << msg << ": " << size << endl;
  for (int i = 0; i < size; ++i)
  {
    cout << arr[indexArr[i]] << " ";
  }
  cout << endl;
}

/*
helper function to sort candidate polyid arrays
*/
// ------------------------------------------
void countSort(int *input, int *indexArr, int n)
{
  int *output;      // output array
  int *lsVals;      // save least significant digits of the input
  int *outputIndex; // current sorted index order
  output = (int *)malloc((n) * sizeof(int));
  lsVals = (int *)malloc((n) * sizeof(int));
  outputIndex = (int *)malloc((n) * sizeof(int));

  int i, count[10] = {0};

  // Store count of occurrences in count[]
  for (i = 0; i < n; i++)
  {
    lsVals[i] = input[i] % 10;
    count[lsVals[i]]++;
    // currOutputIndex[i]=outputIndex[i];
  }

  // Change count[i] so that count[i] now contains actual
  //  position of this digit in output[]
  for (i = 1; i < 10; i++)
    count[i] += count[i - 1];

  // Build the output array
  for (i = n - 1; i >= 0; i--)
  {
    output[count[lsVals[i]] - 1] = input[i];
    outputIndex[count[lsVals[i]] - 1] = indexArr[i];
    count[lsVals[i]]--;
  }

  // Copy the output array to arr[], so that arr[] now
  // contains sorted numbers according to current digit
  for (i = 0; i < n; i++)
  {
    input[i] = output[i] / 10;
    indexArr[i] = outputIndex[i];
    // cout<<"i "<<i<<" "<<output[i]<<" "<<input[i]<<endl;
  }
}

// The main function to that sorts arr[] of size n using
// Radix Sort
void radixsort(int *data, int *outputIndex, int n, int m)
{
  // copy input data to use as the output array
  int *input;
  input = (int *)malloc((n) * sizeof(int));
  // if(DEBUG) cout<<"allocated for sort array"<<endl;

  for (int i = 0; i < n; ++i)
  {
    input[i] = data[i];
    outputIndex[i] = i;
  }
  // if(DEBUG) cout<<"initialize sort array"<<endl;
  // Do counting sort for every digit. Note that instead
  // of passing digit number, exp is passed. exp is 10^i
  // where i is current digit number
  for (int exp = 1; m / exp > 0; exp *= 10)
  {
    // cout<<n<<" round "<<exp<<"\n\n"<<endl;
    countSort(input, outputIndex, n);
    // cout<<"round "<<exp<<"\n\n"<<endl;
  }
  free(input);
}
// ------------------------------------------

/*
remove duplicate intersections at a given edge
*/
void removeDuplicatesAtEdge(int *bPolLineIdsDup, int *cPolLineIdsDup, int *intersectTypesDup, int *bPolLineIdsSortedIndexDup, int intersectCountDup,
                            int *bPolLineIds, int *cPolLineIds, int *intersectTypes, int *bPolLineIdsSortedIndex, int &intersectCount)
{

  int *dupIntersections;
  dupIntersections = (int *)calloc(intersectCountDup, sizeof(int));
  int j = 0, sid, sid2, exitid;

  // dupIntersections[bPolLineIdsSortedIndexDup[0]]=0;
  for (int i = 0; i < intersectCountDup; ++i)
    dupIntersections[i] = 0;
  // cout<<bPolLineIdsSortedIndexDup[0]<<" dup "<<dupIntersections[bPolLineIdsSortedIndexDup[0]]<<endl;

  for (int id = 0; id < intersectCountDup;)
  {
    sid = bPolLineIdsSortedIndexDup[id];

    // cout<<bPolLineIdsDup[sid]<<" //// "<<cPolLineIdsDup[sid]<<endl;
    if (dupIntersections[sid] != 100)
    {
      dupIntersections[sid] = 0;
      // cout<<sid<<" dup "<<dupIntersections[sid]<<endl;
      exitid = id + 1;
      for (int id2 = id + 1; id2 < intersectCountDup; ++id2)
      {
        sid2 = bPolLineIdsSortedIndexDup[id2];
        if (bPolLineIdsDup[sid] == bPolLineIdsDup[sid2] && cPolLineIdsDup[sid] == cPolLineIdsDup[sid2])
        {
          dupIntersections[sid2] = 100;
          // cout<<sid2<<" dup "<<dupIntersections[sid2]<<endl;
        }
        // exit search when bPolId is different
        if (bPolLineIdsDup[sid] != bPolLineIdsDup[sid2])
        {
          // exitid=id2;
          id2 = intersectCountDup;
        }
      }
      id = exitid;
      bPolLineIds[j] = bPolLineIdsDup[sid];
      cPolLineIds[j] = cPolLineIdsDup[sid];
      bPolLineIdsSortedIndex[j] = j;
      intersectTypes[j] = intersectTypesDup[sid];
      j++;
    }
    else
    {
      ++id;
    }
  }

  intersectCount = j;
}

/*
helper function find max
*/
int findMax(int val1, int val2)
{
  if (val1 > val2)
    return val1;
  return val2;
}

/*
Check degen type
*/
bool isNonDegenerateCase(char polyName, int type)
{
  // cout<<"is "<<polyName<<" "<<type<<endl;
  if (polyName == 'b' && (type == 1 || type == 3 || type == 5 || type == 7))
    return true;
  else if (polyName == 'c' && (type == 1 || type == 2 || type == 6))
    return true;

  return false;
}

/*
Function to count intersections at each edge
*/
void countEdgeIntersections(int *intersectEdges, int *edgeIdsSortedIndex, int cStart, int *intersectTypes, int intersectSize,
                            char polyName, int *allIntersectCounts, int *nonDegenIntersectCounts, int *edgeIds, int polysize,
                            int *uniqueIntersectEdgeCount)
{
  int i;
  int j = 0, edgeCount = 1, nonDegeneEdgeCount = 0;
  int *edgeCounts;
  int *nonDegenEdgeCounts;
  edgeCounts = (int *)malloc((intersectSize) * sizeof(int));
  nonDegenEdgeCounts = (int *)malloc((intersectSize) * sizeof(int));

  // int edgeIds[intersectSize];
  int current = intersectEdges[edgeIdsSortedIndex[0]];
  edgeCounts[0] = edgeCount;
  edgeIds[0] = current;

  if (isNonDegenerateCase(polyName, intersectTypes[edgeIdsSortedIndex[0]]))
  {
    nonDegeneEdgeCount = 1;
  }
  nonDegenEdgeCounts[0] = nonDegeneEdgeCount;

  // assign a default value to check if the last edge count is assigned or not
  // edgeIds[intersectSize-1]=-1;

  // count for only for intersections
  for (i = 1; i < intersectSize; ++i)
  {
    if (intersectEdges[edgeIdsSortedIndex[i]] == current)
    {
      edgeCount++;
      // cout<<i<<" "<<j<<" =="<<edgeCount<<endl;
      if (isNonDegenerateCase(polyName, intersectTypes[edgeIdsSortedIndex[i]]))
        nonDegeneEdgeCount++;
    }
    else
    {
      edgeCounts[j] = edgeCount;
      nonDegenEdgeCounts[j] = nonDegeneEdgeCount;
      edgeIds[j + 1] = intersectEdges[edgeIdsSortedIndex[i]];

      // cout<<j<<" cc "<<edgeIds[j]<<" "<<edgeCounts[j]<<endl;
      current = intersectEdges[edgeIdsSortedIndex[i]];
      j++;
      edgeCount = 1;
      nonDegeneEdgeCount = 0;
      if (isNonDegenerateCase(polyName, intersectTypes[edgeIdsSortedIndex[i]]))
      {
        nonDegeneEdgeCount = 1;
      }
      // cout<<polyName<<" +++ "<<intersectTypes[edgeIdsSortedIndex[i]]<<" "<<nonDegeneEdgeCount<<endl;
    }
    if (i == intersectSize - 1)
    {
      edgeCounts[j] = edgeCount;
      nonDegenEdgeCounts[j] = nonDegeneEdgeCount;
    }
  }

  j = 0;
  // copy counts to the all edge array
  for (i = 0; i < polysize; ++i)
  {
    if ((i+cStart) == edgeIds[j] / 2)
    {
      allIntersectCounts[i] = edgeCounts[j];
      nonDegenIntersectCounts[i] = nonDegenEdgeCounts[j];
      // cout<<"---" <<i<<" "<<allIntersectCounts[i]<<" "<<nonDegenIntersectCounts[i]<<endl;
      j++;
    }
    else
    {
      allIntersectCounts[i] = 0;
      nonDegenIntersectCounts[i] = 0;
    }
  }
  *uniqueIntersectEdgeCount = j;
  free(edgeCounts);
  free(nonDegenEdgeCounts);
}

/*
copy data to the linked list
*/
void copyToLinkedList(REAL *intersectedB, REAL *intersectedC, int bListId, int cType,
                      int *neighborsB, int *neighborsC,
                      int *alphaValuesB, int *alphaValuesC,
                      int bIntersectedCount, int cIntersectedCount,
                      int *bInitLabels, int *cInitLabels)
{

  vertex *tmpVertex, *current;
  vertex **bVertexPointers;
  // -------------------------------------------------------------------------------------------
  // Polygon P: (bPolysList)insert intersection vertices and change alpha value in the degenerate cases
  // -------------------------------------------------------------------------------------------
  int i = 0, pi = 0;
  int intersectionPArrayMax = bIntersectedCount * 2;
  bVertexPointers = new vertex *[bIntersectedCount];

  vertex *V = bPolysList[bListId].root;
  for (int ii = 0; ii <= bPolysList[bListId].size; ii++)
  {
    current = V;
    while (*(intersectedB + (i % intersectionPArrayMax)) != V->p.x || *(intersectedB + ((i + 1) % intersectionPArrayMax)) != V->p.y)
    {
      tmpVertex = new vertex(*(intersectedB + i), *(intersectedB + i + 1));
      tmpVertex->label = (IntersectionLabel)(*(bInitLabels + (i / 2)));
      if(tmpVertex->label==-100) tmpVertex->label==0;
      tmpVertex->source = false;
      tmpVertex->intersection = true;
      tmpVertex->next = current;
      current->prev->next = tmpVertex;
      tmpVertex->prev = current->prev;
      current->prev = tmpVertex;
      bVertexPointers[pi++] = tmpVertex;
      i += 2;
    }
    if (ii < bPolysList[bListId].size)
    {
      bVertexPointers[pi] = V;
      pi++;
      V->label = (IntersectionLabel)(*(bInitLabels + (i / 2)));
      if(V->label==-100) V->label==0;
      if (*(alphaValuesB + (i / 2)) != -100)
      {
        V->intersection = true;
      }
    }

    i += 2;
    V = current->next;
  }
  // cout<<"b done"<<endl;
  // -------------------------------------------------------------------------------------------
  // Polygon Q: (cPolysList)insert intersection vertices and change alpha value in the degenerate cases
  // -------------------------------------------------------------------------------------------
  i = 0;
  V = cPolysList[cType].root;
  int intersectionQArrayMax = cIntersectedCount * 2;
  for (int ii = 0; ii <= cPolysList[cType].size; ++ii)
  {
    current = V;
    while (*(intersectedC + (i % intersectionQArrayMax)) != V->p.x || *(intersectedC + ((i + 1) % intersectionQArrayMax)) != V->p.y)
    {
      tmpVertex = new vertex(*(intersectedC + i), *(intersectedC + i + 1));
      tmpVertex->label = (IntersectionLabel)(*(cInitLabels + (i / 2)));
      if(tmpVertex->label==-100) tmpVertex->label==0;
      tmpVertex->source = false;
      tmpVertex->intersection = true;
      tmpVertex->next = current;
      current->prev->next = tmpVertex;
      tmpVertex->prev = current->prev;
      current->prev = tmpVertex;
      tmpVertex->neighbour = bVertexPointers[(*(neighborsC + (i / 2))) - 1];
      bVertexPointers[(*(neighborsC + (i / 2))) - 1]->neighbour = tmpVertex;
      i += 2;
    }
    if (ii < cPolysList[cType].size)
    {
      V->label = (IntersectionLabel)(*(cInitLabels + (i / 2)));
      if(V->label==-100) V->label==0;
      if (*(neighborsC + (i / 2)) != -100)
      {
        V->intersection = true;
        V->neighbour = bVertexPointers[(*(neighborsC + (i / 2))) - 1];
        bVertexPointers[(*(neighborsC + (i / 2))) - 1]->neighbour = V;
      }
    }
    i += 2;
    V = current->next;
  }
  //  if(DEBUG) printf("Copying completed");
}

/*
calculate MBR of a polygon
*/
vector<double> getMBR(REAL *poly, int pSize)
{
  vector<double> mbr;
  REAL maxX, minX, maxY, minY;
  REAL px, py;
  maxX = poly[0];
  maxY = poly[1];
  minX = poly[0];
  minY = poly[1];
  cout << "loop " << endl;
  for (int i = 0; i < pSize * 2; i += 2)
  {
    px = poly[i];
    py = poly[i + 1];
    if (maxX < px)
    {
      maxX = px;
    }
    if (minX > px)
    {
      minX = px;
    }
    if (maxY < py)
    {
      maxY = py;
    }
    if (minY > py)
    {
      minY = py;
    }
    //   cout << setprecision (15) << V->p << endl;
  }
  mbr.push_back(minX);
  mbr.push_back(minY);
  mbr.push_back(maxX);
  mbr.push_back(maxY);
  // uncomment to print in file format
  printf("%f %f %f %f\n", minX, minY, maxX, maxY);
  return mbr;
}
/*
get CMBR for PP and QQ
*/
void getCMBR(REAL *cmbr, REAL *bPoly, REAL *cPoly, int bSize, int cSize)
{
  vector<REAL> PPMBR;
  vector<REAL> QQMBR;
  // double *cmbr; //minx, miny, maxx, maxy
  REAL minX, minY, maxX, maxY;

  PPMBR = getMBR(bPoly, bSize);
  QQMBR = getMBR(cPoly, cSize);

  #ifdef DG3
    cout << "MBR_P [" << PPMBR[0] << ", " << PPMBR[1] << ", " << PPMBR[2] << ", " << PPMBR[3] << endl;
    cout << "MBR_Q [" << QQMBR[0] << ", " << QQMBR[1] << ", " << QQMBR[2] << ", " << QQMBR[3] << endl;
  #endif

  // check intersection between MBRs
  if (PPMBR[0] > QQMBR[2] || PPMBR[2] < QQMBR[0])
  {
    printf("No Overlap between polygons\n");
    exit(0);
  }
  if (PPMBR[1] > QQMBR[3] || PPMBR[3] < QQMBR[1])
  {
    printf("No Overlap between polygons\n");
    exit(0);
  }

  cmbr[0] = max(PPMBR[0], QQMBR[0]);
  cmbr[1] = max(PPMBR[1], QQMBR[1]);
  cmbr[2] = min(PPMBR[2], QQMBR[2]);
  cmbr[3] = min(PPMBR[3], QQMBR[3]);
  #ifdef DG3
    cout << "CMBR [" << cmbr[0] << ", " << cmbr[1] << ", " << cmbr[2] << ", " << cmbr[3] << endl;
  #endif
}

#ifdef GPU
// read data into polygon arrays and interval array
// gpc_read_polygon(bFile, cFile, &bPoly, &cPoly, bSize, cSize, &intervals, ic, minVal, maxVal);
void hybrid_clipping(int argc, char *argv[])
{
  FILE *bFile, **cFile;
  int cPolyCount=argc-2;
  cFile = (FILE **)malloc((cPolyCount) * sizeof(FILE));
  for (int i = 0; i <cPolyCount; ++i)
    *(cFile + i) = fopen(argv[2 + i], "r");

  // FILE *bFile, *cFile;
  // cFile=fopen(argv[2], "r");

  bFile = fopen(argv[1], "r");
  REAL *bPoly, *cPoly, minVal, maxVal;
  Edge *intervals;
  int bSize = 0, cSize = 0, clSize, ic, *cSizes;
  double start_time, run_time;
  cSizes = (int *)malloc((argc - 1) * sizeof(int));

  #ifdef TIME
  high_resolution_clock::time_point cp1, cp2, cp3, cp4, cp6, cp7, cp8, cp9, cp10, cp11, cp12, cp13, cp14, cp15, cp16, cp17, cp18, cp34, cp35, cp19, cp20, cp21, cp22, cp23;
  #endif

  if (argc < 3)
  {
    cout << "3 parameters required\n./seg <base polygon path> <subject polygon path>" << endl;
    exit(0);
  }
  // read data into polygon arrays and interval array
  // gpc_read_polygon(bFile, cFile, &bPoly, &cPoly, bSize, cSize, &intervals, ic, minVal, maxVal);
  // cSizes[0]=0;
  // cSizes[1]=cSize;
  gpc_read_multi_polygon(bFile, cFile, cPolyCount, &bPoly, &cPoly, bSize, clSize, cSizes, &intervals, ic, minVal, maxVal);

  // cout<<"ic="<<ic<<" cal ic="<<bSize+1+clSize+cPolyCount<<endl;

  // duplicates the base polygon allowing the clipping layer polygons to work on clipping independently
  duplicateBasePoly(cPolyCount);

  // print polygons
  // -------------------------------------------------
  // cout << "Base polygon " << bSize << endl;
  // for (int i=0; i<bSize*2; i+=2)
  // {
  //   cout<<i/2<<": "<<*(bPoly+i)<<", "<<*(bPoly+i+1)<<endl;
  // }

  // cout << "\nOverlay polygon " << cSize << endl;
  // for (int i=0; i<cSize*2; i+=2)
  // {
  //   cout<<i/2<<": "<<*(cPoly+i)<<", "<<*(cPoly+i+1)<<endl;
  // }
  // cout << endl;

  // -------------------------------------------------
  #ifdef TIME
  cp1 = high_resolution_clock::now();
  #endif

  // OpenACC initializing 
  #ifdef TIME
  cp19 = high_resolution_clock::now();
  #endif
  #pragma acc kernels
  {

  }
  #ifdef TIME
  cp20 = high_resolution_clock::now();
  #endif

  // create empty tree
  int numEdges = bSize + clSize;
  SegTree sTree(numEdges, ic, minVal, maxVal);

  sTree.copyInputDataToGPU(bPoly, cPoly, bSize, clSize, (clSize+cPolyCount), cPolyCount, ic, intervals); 

  #ifdef TIME
  cp21 = high_resolution_clock::now();
  #endif
  // make intervals a class varaible and use it across all places. Copyin it only one time. Next implement above

  // find CMBR
  // REAL cmbr[4];
  // getCMBR(cmbr, bPoly, cPoly, bSize, cSize);

  // sTree.intervalConstructionParallel
  //   (bPoly, cPoly, bSize, cSize,
  //   &intervals);
  // ----------------------------
  #ifdef DG
  cout << "Empty segment tree with " << numEdges << " edges created" << endl;
  #endif

  // sTree.constructElementaryIntervalsArrayCMBR(intervals, cmbr);
  sTree.constructElementaryIntervalsArrayGPU();
  // sTree.constructElementaryIntervalsArrayGPUParallel(intervals);
  #ifdef TIME
  cp34 = high_resolution_clock::now();
  #endif
  sTree.initHybrid();
  #ifdef TIME 
  cp35 = high_resolution_clock::now();
  #endif

  #ifdef DG
  cout << "Tree initialized. Tree size: " << sTree.treeSize << endl;
  #endif

  sTree.POLY_TYPE_COUNT = argc - 1;
  sTree.constructTreeIntervalArray();
  #ifdef TIME
  cp14 = high_resolution_clock::now();
  #endif

  #ifdef DG
  cout << "Build finished" << endl;
  #endif
  // ----------------------------

  // ----------------------------
  // vector<Edge> *elementaryIntervals = sTree.getElementaryIntervals(intervals);
  // if(DG) cout<<elementaryIntervals->size()<<" elementary intervals copied"<<endl;
  // // copy to array
  // sTree.copyToEleArray(elementaryIntervals);

  // // initialize segment tree
  // sTree.init2(elementaryIntervals);
  // // if(DG) cout<<"Leaf level node count: "<<sTree.numLeaves<<endl;
  // if(DG) cout<<"Tree initialized. Tree size: "<<sTree.treeSize<<endl;

  // sTree.POLY_TYPE_COUNT=argc-1;
  // sTree.constructTreeIntervalArray();
  // if(DG) cout<<"Tree interval array built"<<endl;
  // // build segment tree skeleton
  // sTree.buildSkeleton2(1, 0, elementaryIntervals->size()-1);

  // if(TIME) cp14 = high_resolution_clock::now();

  // if(DG) cout<<"Build finished"<<endl;
  // ----------------------------

  int maxId;
  if(bSize>=(clSize+cPolyCount))
    maxId=bSize;
  else
    maxId=(clSize+cPolyCount);
  

  // sTree.constructCoverList(intervals, bPoly, cPoly, bSize, (clSize + cPolyCount));
  // sTree.insertAllToEdgeListParallel(intervals, maxId);


  // omp_set_nested(1);
  // #pragma omp parallel sections num_threads(2)
  {
    // prepare prefix sum arrays for cover lists sizes and end list sizes
    // #pragma omp section
    {
      #ifdef TIME 
      cp10 = high_resolution_clock::now();
      #endif
      #ifdef DG 
      cout << "End lists starting.." << endl;
      #endif

      sTree.insertAllToEdgeListParallel(maxId);

      #ifdef TIME 
      cp11 = high_resolution_clock::now();
      #endif
      #ifdef DG 
      cout << "End lists done" << endl;
      #endif
    }

    // #pragma omp section
    {
      #ifdef TIME 
      cp8 = high_resolution_clock::now();
      #endif
      #ifdef DG
      cout << "Cover lists starting..." << endl;
      #endif

      // sTree.insertAllParallel(intervals);
      sTree.constructCoverListGPU();
      // sTree.constructCoverListMulticore(intervals);
      // sTree.constructCoverListSerial(intervals);
      // sTree.constructCoverEndLists(bPoly, cPoly, bSize, (clSize + cPolyCount));
      #ifdef DG
      cout << "Cover lists built" << endl;
      #endif

      #ifdef TIME 
      cp9 = high_resolution_clock::now();
      #endif
    }
  }

  // sTree.buildInternalEndListsParallel(bSize);

  // convert tree coverlists  to array and offload data to GPU
  /*#ifdef TIME 
  cp12 = high_resolution_clock::now();
  #endif

  sTree.countEndlistCoverListSizes();
  #ifdef DG
  cout << "End list sizes saved" << endl;
  #endif
  
  sTree.copyNodeCoverListToGPU(bPoly, cPoly, bSize, (clSize + cPolyCount));
  #ifdef TIME 
  cp13 = high_resolution_clock::now();
  #endif

  #ifdef DG
  cout << "Cover lists copied to GPU" << endl;
  #endif



  // int lts = sTree.treeSize;
  // for (int i = 0; i < lts; ++i) {
  //   // omp_destroy_lock(&lock[i] );
  //   omp_destroy_lock(&sTree.lock[i]);
  // }

  #ifdef TIME 
  cp2 = high_resolution_clock::now();
  #endif

  #ifdef TIME 
  cp15 = high_resolution_clock::now();
  #endif*/

  int ncount = sTree.treeSize; // node 0 is invalid. True node count is treeSeze-2
  int bType=0; //bType is base polygon type. Always 0. bListId is the base poly id from the linked list

  // ====================================================================================
  // **************** when b<c swap b and c for that case only. *************
  // ====================================================================================

  // loop over all clipping layer polygons
  for (int bListId=0, cListId=0, cType=1; cType<=cPolyCount; cType++, bListId++, cListId++){ 
    int *iCount;
    iCount = (int *)malloc((ncount) * sizeof(int));

    cSize=cSizes[cType]-cSizes[cListId]-1;

    sTree.countIntersectionsParallel(bPoly, cPoly, bSize, cSize, cPolyCount,
                                    bType, cType,
                                    iCount);
    // #pragma acc update self(iCount[0 : ncount])
    // for(int ii=0; ii<ncount; ++ii){
    //   if(iCount[ii]!=0) cout<<iCount[ii]<<", ";
    // }
    // cout<<endl;

    #ifdef TIME 
    cp16 = high_resolution_clock::now();
    #endif
    #ifdef DG
    cout << "Count intersections" << endl;
    #endif

    // thrust::exclusive_scan(iCount, iCount + ncount, iCount);
    // #pragma acc update device(iCount[0 : ncount])

    #pragma acc host_data use_device(iCount)
    thrustDeviceExclusiveScan(iCount, ncount);

    #ifdef TIME 
    cp22 = high_resolution_clock::now();
    #endif
    #pragma acc update self(iCount[ncount-1])
    #ifdef TIME 
    cp23 = high_resolution_clock::now();
    #endif

    // --------------------------------------------------------
    #ifdef DG
    cout << "Prefix sum: intersection count at nodes " << iCount[ncount - 1] << endl;
    #endif


    int intersectCountDup = iCount[ncount - 1];
    int *bPolLineIdsDup, *cPolLineIdsDup, *bPolLineIdsSortedIndexDup, *intersectTypesDup;

    bPolLineIdsDup = (int *)malloc((intersectCountDup) * sizeof(int));
    cPolLineIdsDup = (int *)malloc((intersectCountDup) * sizeof(int));
    bPolLineIdsSortedIndexDup = (int *)malloc((intersectCountDup) * sizeof(int));
    intersectTypesDup = (int *)malloc((intersectCountDup) * sizeof(int));
      
    #ifdef DG
    cout << "Clip start " << intersectCountDup << endl;
    #endif

    #ifdef TIME 
    cp17 = high_resolution_clock::now();
    #endif
    // sTree.polyClipArrayParallel(bPoly, cPoly, bSize, cSize,
    //                             bType, cType,
    //                             iCount, bPolLineIdsDup, cPolLineIdsDup, intersectTypesDup);
    sTree.polyClipArrayGPUParallel
          (bPoly, cPoly, bSize, cSize,
          bType, cType,
          iCount, bPolLineIdsDup, cPolLineIdsDup,intersectTypesDup);
    #ifdef DG
    cout << "Clip end" << endl;
    #endif
    #ifdef TIME 
    cp18 = high_resolution_clock::now();
    #endif

    radixsort(bPolLineIdsDup, bPolLineIdsSortedIndexDup, intersectCountDup, bSize * 2);

    int *bPolLineIds, *cPolLineIds;
    int *intersectTypes, *bPolLineIdsSortedIndex, intersectCount;
    bPolLineIds = (int *)malloc((intersectCountDup) * sizeof(int));
    cPolLineIds = (int *)malloc((intersectCountDup) * sizeof(int));
    intersectTypes = (int *)malloc((intersectCountDup) * sizeof(int));
    bPolLineIdsSortedIndex = (int *)malloc((intersectCountDup) * sizeof(int));

    // remove intersection duplicates
    removeDuplicatesAtEdge(bPolLineIdsDup, cPolLineIdsDup, intersectTypesDup, bPolLineIdsSortedIndexDup, intersectCountDup,
                          bPolLineIds, cPolLineIds, intersectTypes, bPolLineIdsSortedIndex, intersectCount);

    #ifdef DG
    cout << intersectCountDup << " intersections found. " << intersectCountDup - intersectCount << " duplicate(s) removed. " << intersectCount << " remaining." << endl;
    #endif

    int *cPolLineIdsSortedIndex;
    cPolLineIdsSortedIndex = (int *)malloc((intersectCountDup) * sizeof(int));

    // sort cPolLineIds
    radixsort(cPolLineIds, cPolLineIdsSortedIndex, intersectCount, cSize * 2);
    #ifdef DG
    cout << "Sorted candidate edges by ID" << endl;
    #endif

    // get intersection counts at each edge for both polygons
    int *bAllIntersectCounts, *cAllIntersectCounts;
    int *bNonDegenIntersectCounts, *cNonDegenIntersectCounts;
    bAllIntersectCounts = (int *)malloc(((bSize + 1)) * sizeof(int));
    cAllIntersectCounts = (int *)malloc(((cSize + 1)) * sizeof(int));
    bNonDegenIntersectCounts = (int *)malloc(((bSize + 1)) * sizeof(int));
    cNonDegenIntersectCounts = (int *)malloc(((cSize + 1)) * sizeof(int));

    // default values for the extra location to use for prefix sum
    bAllIntersectCounts[bSize] = 0;
    cAllIntersectCounts[cSize] = 0;
    bNonDegenIntersectCounts[bSize] = 0;
    cNonDegenIntersectCounts[cSize] = 0;

    int *bPolLineUniqueIds;
    int *cPolLineUniqueIds;
    bPolLineUniqueIds = (int *)malloc((intersectCount) * sizeof(int));
    cPolLineUniqueIds = (int *)malloc((intersectCount) * sizeof(int));

    int uniqueIntersectEdgeCountB, uniqueIntersectEdgeCountC;
    countEdgeIntersections(bPolLineIds, bPolLineIdsSortedIndex, 0, intersectTypes, intersectCount,
                          'b', bAllIntersectCounts, bNonDegenIntersectCounts, bPolLineUniqueIds, bSize,
                          &uniqueIntersectEdgeCountB);

    #ifdef DG
    cout << "Base polygon edges counted" << endl;
    #endif
    countEdgeIntersections(cPolLineIds, cPolLineIdsSortedIndex, cSizes[cListId], intersectTypes, intersectCount,
                          'c', cAllIntersectCounts, cNonDegenIntersectCounts, cPolLineUniqueIds, cSize,
                          &uniqueIntersectEdgeCountC);

    #ifdef DG
    cout << "Clipping polygon edges counted" << endl;
    #endif

    // prefix sum of the edge counts
    // #pragma acc routine seq
    // thrust::exclusive_scan(thrust::host, bAllIntersectCounts, bAllIntersectCounts+bSize+1, bAllIntersectCounts);
    thrust::exclusive_scan(bAllIntersectCounts, bAllIntersectCounts + bSize + 1, bAllIntersectCounts);
    thrust::exclusive_scan(cAllIntersectCounts, cAllIntersectCounts + cSize + 1, cAllIntersectCounts);
    thrust::exclusive_scan(bNonDegenIntersectCounts, bNonDegenIntersectCounts + bSize + 1, bNonDegenIntersectCounts);
    thrust::exclusive_scan(cNonDegenIntersectCounts, cNonDegenIntersectCounts + cSize + 1, cNonDegenIntersectCounts);

    #ifdef DG
    cout << "Edge prefix sums calculated" << endl;
    #endif

    int bNonDegenCount = bNonDegenIntersectCounts[bSize];
    int cNonDegenCount = cNonDegenIntersectCounts[cSize];
    int bIntersectedCount = bSize + bNonDegenCount;
    int cIntersectedCount = cSize + cNonDegenCount;

    #ifdef DG
    cout << "Poly B: non-degen count=" << bNonDegenCount << " total=" << bIntersectedCount;
    #endif
    #ifdef DG
    cout << " Poly C: non-degen count=" << cNonDegenCount << " total=" << cIntersectedCount << endl;
    #endif

    // cout<<bSize<<" "<<cSize<<" "<<bSize+cSize<<" "<<sTree.treeSize<<" "<<intersectCount<<" ";

    // save only intersections in new arrays
    // remove intersectType saving here. It is completed at polyClipArrayParallel()
    REAL *intersections;
    intersections = (REAL *)malloc((2 * intersectCount) * sizeof(REAL));

    int *intersectAlphaValuesB, *intersectAlphaValuesC;
    intersectAlphaValuesB = (int *)malloc((intersectCount) * sizeof(int));
    intersectAlphaValuesC = (int *)malloc((intersectCount) * sizeof(int));

    sTree.saveIntersectionsParallel(bPoly, cPoly,
                                    bPolLineIds, cPolLineIds, intersectCount,
                                    intersections, intersectTypes,
                                    intersectAlphaValuesB, intersectAlphaValuesC);
    #ifdef DG
    cout << "Intersections saved" << endl;
    #endif

    // sort saved edges in the common intersections array
    int *sortedIndiciesB, *sortedIndiciesC;
    sortedIndiciesB = (int *)malloc((intersectCount) * sizeof(int));
    sortedIndiciesC = (int *)malloc((intersectCount) * sizeof(int));

    sTree.sortIntersectionsAtEdgesParallel(intersectAlphaValuesB, bPolLineIdsSortedIndex, 0,
                                          bAllIntersectCounts, bNonDegenIntersectCounts, bPolLineUniqueIds, uniqueIntersectEdgeCountB, intersectCount,
                                          sortedIndiciesB);

    sTree.sortIntersectionsAtEdgesParallel(intersectAlphaValuesC, cPolLineIdsSortedIndex, cSizes[cListId],
                                          cAllIntersectCounts, cNonDegenIntersectCounts, cPolLineUniqueIds, uniqueIntersectEdgeCountC, intersectCount,
                                          sortedIndiciesC);
    #ifdef DG
    cout << "Sort at edge completed" << endl;
    #endif

    REAL *intersectedB, *intersectedC;
    intersectedB = (REAL *)malloc((2 * bIntersectedCount) * sizeof(REAL));
    intersectedC = (REAL *)malloc((2 * cIntersectedCount) * sizeof(REAL));

    int *alphaValuesB, *alphaValuesC;
    int *neighborsB, *neighborsC;
    alphaValuesB = (int *)malloc((bIntersectedCount) * sizeof(int));
    alphaValuesC = (int *)malloc((cIntersectedCount) * sizeof(int));
    neighborsB = (int *)malloc((bIntersectedCount) * sizeof(int));
    neighborsC = (int *)malloc((cIntersectedCount) * sizeof(int));

    // data coping not rquired at this point. Existing CPU data update can be done with the available data
    sTree.copyIntersectionsParallel(bPoly, cPoly, bSize, cSize, cSizes[cListId],
                                    bPolLineIds, cPolLineIds, intersectCount,
                                    bPolLineIdsSortedIndex, cPolLineIdsSortedIndex,
                                    sortedIndiciesB, sortedIndiciesC,
                                    intersections, intersectTypes,
                                    bAllIntersectCounts, cAllIntersectCounts,
                                    bNonDegenIntersectCounts, cNonDegenIntersectCounts,
                                    intersectAlphaValuesB, intersectAlphaValuesC,
                                    intersectedB, intersectedC,
                                    neighborsB, neighborsC,
                                    alphaValuesB, alphaValuesC);

    #ifdef DG
    cout << "Intersections copied to source arrays" << endl;
    #endif

    #ifdef TIME 
    cp3 = high_resolution_clock::now();
    #endif

    // initial labeling
    int *bInitLabels, *cInitLabels;
    bInitLabels = (int *)malloc((bIntersectedCount) * sizeof(int));
    cInitLabels = (int *)malloc((cIntersectedCount) * sizeof(int));

    sTree.calculateInitLabelParallel(bSize, bNonDegenIntersectCounts,
                                    intersectedB, intersectedC, alphaValuesB,
                                    neighborsB,
                                    bIntersectedCount, cIntersectedCount, bInitLabels, cInitLabels);

    #ifdef DG
    cout << "Initial labels calculated" << endl;
    #endif

    // copy data to the linked lists
    copyToLinkedList(intersectedB, intersectedC, bListId, cListId,
                    neighborsB, neighborsC,
                    alphaValuesB, alphaValuesC,
                    bIntersectedCount, cIntersectedCount,
                    bInitLabels, cInitLabels);
    #ifdef DG
    cout << "Copied to linked lists" << endl;
    #endif
    // PHASE: 3
    labelIntersections2(bListId, cListId);

    // PHASE: 4
    createResult2();
    #ifdef TIME 
    cp4 = high_resolution_clock::now();
    #endif

    // if(TIME) end = high_resolution_clock::now();
    // post-processing
    cleanUpResult2();

    // write output polygon
    #ifdef SAVE
      cout << "R ";
      savePolygon(rPolyList, outputFile);
    #endif
    // savePolygon2(rPolyList, outputFile);
    

    //clean arrays
    rPolyList.clear();

    // cout << "Clipping complete" << endl;
  }

  // clean arrays
  free(intervals);
  free(bPoly);
  free(cPoly);

  // free(intersectedB);
  // free(intersectedC);
  // free(alphaValuesB);
  // free(alphaValuesC);
  // free(neighborsB);
  // free(neighborsC);

  // free(intersections);
  // free(intersectAlphaValuesB);
  // free(intersectAlphaValuesC);

  // free(bPolLineIdsDup);
  // free(cPolLineIdsDup);
  // free(bPolLineIdsSortedIndexDup);
  // free(intersectTypesDup);
  // free(bPolLineIds);
  // free(cPolLineIds);
  // free(intersectTypes);
  // free(bPolLineIdsSortedIndex);
  // free(cPolLineIdsSortedIndex);
  // free(bAllIntersectCounts);
  // free(cAllIntersectCounts);
  // free(bNonDegenIntersectCounts);
  // free(cNonDegenIntersectCounts);
  // free(bPolLineUniqueIds);
  // free(cPolLineUniqueIds);
  // free(sortedIndiciesB);
  // free(sortedIndiciesC);
  // free(bInitLabels);
  // free(cInitLabels);

  bPolysList.clear();
  cPolysList.clear();
  // rPolyList.clear();

  #ifdef DT_INFO 
  // if (DT_INFO) {
    // Calculating total time taken by the program.
    auto dTotalTime = duration_cast<microseconds>(cp4 - cp1);
    auto dOpenAccInit = duration_cast<microseconds>(cp20 - cp19);
    
    auto dSegTreeConstruct = duration_cast<microseconds>(cp9 - cp1);
    auto dCopyInputData = duration_cast<microseconds>(cp21 - cp20);
    auto dElementaryArrayGpu = duration_cast<microseconds>(cp34 - cp21);
    auto dInit = duration_cast<microseconds>(cp35 - cp34);
    auto dTreeIntervalArray = duration_cast<microseconds>(cp14 - cp35);
    auto dSkeletonConstruct = duration_cast<microseconds>(cp14 - cp1);
    auto dEndlistConstruct = duration_cast<microseconds>(cp11 - cp10);
    auto dCoverListConstruct = duration_cast<microseconds>(cp9 - cp8);
    
    auto dCandidateEdgeFind = duration_cast<microseconds>(cp3 - cp9);
    auto dCandidateEdgeCount = duration_cast<microseconds>(cp16 - cp9);
    auto dCandidateEdgeCountMaxMemcpy = duration_cast<microseconds>(cp23 - cp22);
    auto dCandidateEdgeSave = duration_cast<microseconds>(cp18 - cp16);
    auto dLabeling = duration_cast<microseconds>(cp4 - cp3);


    cout << "PARALLEL TIMING\nAll time in microseconds\nTime: Total: " << fixed << dTotalTime.count() << setprecision(10) << endl;
    cout<<"|-Time: OpenACC initializing: "<<fixed<<dOpenAccInit.count()<<setprecision(10)<<endl;
    cout << "|-Time: Segment tree construction (ALL): " << fixed << dSegTreeConstruct.count() << setprecision(10) << endl;
    // (input data copyin + skeleton + init + tree intervals)
    cout << " |-Time: Segment tree skeleton build: " << fixed << dSkeletonConstruct.count() << setprecision(10) << endl;
    cout<<"  |-Time: copy input data to GPU: "<<fixed<<dCopyInputData.count()<<setprecision(10)<<endl;
    cout << "  |-Time: elementaryArrayGPU: " << fixed << dElementaryArrayGpu.count() << setprecision(10) << endl;
    cout << "  |-Time: Segment tree:init " << fixed << dInit.count() << setprecision(10) << endl;
    cout << "  |-Time: tree interval array " << fixed << dTreeIntervalArray.count() << setprecision(10) << endl;

    cout << " |-Time: end list cosntruct (GPU+CPU): " << fixed << dEndlistConstruct.count() << setprecision(10) << endl;
    cout << " |-Time: cover list costruct: " << fixed << dCoverListConstruct.count() << setprecision(10) << endl;

    cout << "|-Time: Candidate edge find and save: " << fixed << dCandidateEdgeFind.count() << setprecision(10) << endl;
    cout << " |-Time: Candidate edge find count: " << fixed << dCandidateEdgeCount.count() << setprecision(10) << endl;
    cout << " |-Time: Candidate edge count max val copy to CPU: " << fixed << dCandidateEdgeCountMaxMemcpy.count() << setprecision(10) << endl;
    cout << " |-Time: Candidate edge find save ids: " << fixed << (dCandidateEdgeSave.count()-dCandidateEdgeCountMaxMemcpy.count()) << setprecision(10) << endl;

    cout << "|-Time: Labeling and stitch results: " << fixed << dLabeling.count() << setprecision(10) << endl;
  
  #elif TIME
  // } else if (TIME) {
    // Calculating total time taken by the program.
    auto dTotalTime = duration_cast<microseconds>(cp4 - cp1);
    auto dOpenAccInit = duration_cast<microseconds>(cp20 - cp19);
    
    auto dSegTreeConstruct = duration_cast<microseconds>(cp9 - cp1);
    auto dCopyInputData = duration_cast<microseconds>(cp21 - cp20);
    auto dElementaryArrayGpu = duration_cast<microseconds>(cp34 - cp21);
    auto dInit = duration_cast<microseconds>(cp35 - cp34);
    auto dTreeIntervalArray = duration_cast<microseconds>(cp14 - cp35);
    auto dSkeletonConstruct = duration_cast<microseconds>(cp14 - cp1);
    auto dEndlistConstruct = duration_cast<microseconds>(cp11 - cp10);
    auto dCoverListConstruct = duration_cast<microseconds>(cp9 - cp8);
    
    auto dCandidateEdgeFind = duration_cast<microseconds>(cp3 - cp9);
    auto dCandidateEdgeCount = duration_cast<microseconds>(cp16 - cp9);
    auto dCandidateEdgeCountMaxMemcpy = duration_cast<microseconds>(cp23 - cp22);
    auto dCandidateEdgeSave = duration_cast<microseconds>(cp18 - cp16);
    auto dLabeling = duration_cast<microseconds>(cp4 - cp3);


    // input data and intervals copy to GPU, max elementaryArrayIndex copy to CPU, end list partial copy to CPU, end list copy to GPU, cover list max copy to CPU, intersecting candidate edge IDs copy to CPU, candidate edge count max copy to CPU, OpenACC overhead, construct tree skeletone , construct end list , construct cover list, candidate edges find and save, labeling, total time

    cout << fixed << dCandidateEdgeCountMaxMemcpy.count() << setprecision(10) << " " <<dOpenAccInit.count()<<setprecision(10)<<" " << dSkeletonConstruct.count() << setprecision(10) << " " << dEndlistConstruct.count() << setprecision(10) << " " << dCoverListConstruct.count() << setprecision(10) << " " << dCandidateEdgeFind.count() << setprecision(10) << " " << dLabeling.count() << setprecision(10) << " " << dTotalTime.count() << setprecision(10) << endl;

  #endif
  // }
  // add this to destroy gpu data elements
  // #pragma acc exit data delete(coverListSizes, nodeIntervalArray, coverArray, this)
}
#endif


void multicore_clipping(int argc, char *argv[])
{
  // omp_set_num_threads(1);
  // cout<<"*****************Hard-coded for only thread 1*****************"<<endl;
  FILE *bFile, **cFile;
  int cPolyCount=argc-2;
  cFile = (FILE **)malloc((cPolyCount) * sizeof(FILE));
  for (int i = 0; i <cPolyCount; ++i)
    *(cFile + i) = fopen(argv[2 + i], "r");

  // FILE *bFile, *cFile;
  // cFile=fopen(argv[2], "r");

  bFile = fopen(argv[1], "r");
  REAL *bPoly, *cPoly, minVal, maxVal;
  Edge *intervals2, *intervals;
  int bSize = 0, cSize = 0, clSize, ic, *cSizes;
  double start_time, run_time;
  cSizes = (int *)malloc((argc - 1) * sizeof(int));

  #ifdef TIME
  high_resolution_clock::time_point cp1, cp2, cp3, cp4, cp6, cp7, cp8, cp9, cp10, cp11, cp12, cp13, cp14, cp15, cp16, cp17, cp18, cp34, cp35;
  #endif

  if (argc < 3)
  {
    cout << "3 parameters required\n./seg <base polygon path> <subject polygon path>" << endl;
    exit(0);
  }
  // read data into polygon arrays and interval array
  // gpc_read_polygon(bFile, cFile, &bPoly, &cPoly, bSize, cSize, &intervals, ic, minVal, maxVal);
  gpc_read_multi_polygon(bFile, cFile, cPolyCount, &bPoly, &cPoly, bSize, clSize, cSizes, &intervals, ic, minVal, maxVal);

  // duplicates the base polygon allowing the clipping layer polygons to work on clipping independently
  duplicateBasePoly(cPolyCount);
  // print polygons
  // -------------------------------------------------
  // cout << "Base polygon " << bSize << endl;
  // for (int i=0; i<(bSize+1)*2; i+=2)
  // {
  //   cout<<i/2<<": "<<*(bPoly+i)<<", "<<*(bPoly+i+1)<<endl;
  // }

  // cout << "\nOverlay polygon " << clSize << endl;
  // for (int i=0; i<(clSize+cPolyCount)*2; i+=2)
  // {
  //   cout<<i/2<<": "<<*(cPoly+i)<<", "<<*(cPoly+i+1)<<endl;
  // }
  // cout<<"ic="<<ic<<endl;

  // // cSizes[] contains exclusive prefix sum of the individual polygon sizes
  // for (int i = 1; i <=cPolyCount; ++i)
  // {
  //   cout << "Clipping polygon " << i << "-> Size=" << cSizes[i] - cSizes[i - 1]-1 << endl;
  // }

  // cout<<"\nIntervals"<<endl;
  // for(int i=0; i<ic; ++i){
  //   cout<<"i="<<i<<" id="<<intervals[i].id<<" start="<<intervals[i].start<<" end="<<intervals[i].end<<endl;
  // }
  // cout<<endl;


  // // print base and clipping polygons in their linked list form
  // auto next=bPolysList[0].root;
  // cout<<"Base poly "<<bSize<<endl;
  // for(int i=0; i<bSize; ++i){
  //   cout<<"x="<<next->p.x<<", y="<<next->p.y<<" label="<<next->label<<" enx="<<next->enex<<endl;
  //   next=next->next;
  // }
  // next=cPolysList[0].root;
  // cSize=cSizes[1]-1;
  // cout<<"\nClip poly "<<cSize<<endl;
  // for(int i=0; i<cSize; ++i){
  //   cout<<"x="<<next->p.x<<", y="<<next->p.y<<" label="<<next->label<<" enx="<<next->enex<<endl;
  //   next=next->next;
  // }
  // cout<<endl;

  // -------------------------------------------------
  #ifdef TIME
  cp1 = high_resolution_clock::now();
  #endif
  start_time = omp_get_wtime();
  // create empty tree
  int numEdges=bSize+clSize;
  SegTree sTree(numEdges, ic, minVal, maxVal);

  // find CMBR
  // REAL cmbr[4];
  // getCMBR(cmbr, bPoly, cPoly, bSize, cSize);

  // sTree.intervalConstructionParallel
  //   (bPoly, cPoly, bSize, cSize,
  //   &intervals);

  // ----------------------------    
  #ifdef DG 
  cout << "Empty segment tree with " << numEdges << " edges created" << endl;
  #endif
  // sTree.constructElementaryIntervalsArrayCMBR(intervals, cmbr);
  // sTree.constructElementaryIntervalsArrayParallel(intervals);
  sTree.constructElementaryIntervalsArrayMulticore(intervals);
  // sTree.constructElementaryIntervalsArrayGPUParallel(intervals);
  #ifdef TIME
  cp34 = high_resolution_clock::now();
  #endif

  sTree.init();
  #ifdef TIME
  cp35 = high_resolution_clock::now();
  #endif
      
  #ifdef DG 
  cout << "Tree initialized. Tree size: " << sTree.treeSize << endl;
  #endif

  sTree.POLY_TYPE_COUNT = argc - 1;
  sTree.constructTreeIntervalArrayMulticore();
  #ifdef TIME
  cp14 = high_resolution_clock::now();
  #endif

  #ifdef DG
      // Calculating total time taken by the program.
      auto d11 = duration_cast<microseconds>(cp35 - cp34);
      auto d12 = duration_cast<microseconds>(cp14 - cp35);

      cout << "init() " << fixed << d11.count() << setprecision(10) << endl;
      cout << "constructTreeIntervalArrayMulticore() " << fixed << d12.count() << setprecision(10) << endl;
  #endif
    
  #ifdef DG 
  cout << "Build finished" << endl;
  #endif
  // omp_set_nested(1);
  // // omp_set_num_threads(20);
  // // omp_set_teams_thread_limit(20);
  // // #pragma omp parallel
  // #pragma omp parallel sections num_threads(2)
  // {
  //   #pragma omp section
  //   {
  //   if(DEBUG_TIME) cp10 = high_resolution_clock::now();
  //   if(DEBUG) cout<<"End lists starting.."<<endl;
  //   sTree.insertAllToEdgeListMulticore(intervals, bSize);
  //   if(DEBUG_TIME) cp11 = high_resolution_clock::now();
  //   if(DEBUG) cout<<"End lists done"<<endl;
  //   }
  //   #pragma omp section
  //   {
  //   if(DEBUG_TIME) cp8 = high_resolution_clock::now();
  //   if(DEBUG) cout<<"Cover lists starting..."<<endl;
  //   sTree.insertAllParallel(intervals);
  //   if(DEBUG) cout<<"Cover lists built"<<endl;
  //   if(DEBUG_TIME) cp9 = high_resolution_clock::now();

  //   // convert tree coverlists  to array and offload data to GPU
  //   if(DEBUG_TIME) cp12 = high_resolution_clock::now();
  //   sTree.countEndlistCoverListSizesMulticore();
  //   sTree.copyNodeCoverListMulticore();
  //   if(DEBUG_TIME) cp13 = high_resolution_clock::now();
  //   if(DEBUG) cout<<"Cover lists copied to GPU"<<endl;
  //   }
  // }

  int maxId;
  if(bSize>=(clSize+cPolyCount))
    maxId=bSize;
  else
    maxId=(clSize+cPolyCount);

  #ifdef TIME
  cp10 = high_resolution_clock::now();
  #endif

  #ifdef DG 
  cout << "End lists starting.." << endl;
  #endif
  sTree.insertAllToEdgeListMulticore(intervals, maxId);
  #ifdef TIME
  cp11 = high_resolution_clock::now();
  #endif

  #ifdef DG 
  cout << "End lists done" << endl;
  #endif

  #ifdef TIME
  cp8 = high_resolution_clock::now();
  #endif

  #ifdef DG 
  cout << "Cover lists starting..." << endl;
  #endif

  sTree.insertAllParallel(intervals);
    
  #ifdef DG 
  cout << "Cover lists built" << endl;
  #endif

  #ifdef TIME
  cp9 = high_resolution_clock::now();
  #endif

  int lts = sTree.treeSize;
  for (int i = 0; i < lts; ++i)
  {
    // omp_destroy_lock(&lock[i] );
    omp_destroy_lock(&sTree.lock[i]);
  }

  // convert tree coverlists  to array and offload data to GPU
  #ifdef TIME
  cp12 = high_resolution_clock::now();
  #endif
  // sTree.countEndlistCoverListSizesMulticore();
  // sTree.copyNodeCoverListMulticore();
  // if(DEBUG) cout<<"Cover lists copied to GPU"<<endl;
  #ifdef TIME
  cp13 = high_resolution_clock::now();
  #endif

  #ifdef TIME
  cp2 = high_resolution_clock::now();
  #endif
  // -----------------------------------------------------------------------
  // -----------------------------------------------------------------------
  // count the actual node count which are non empty
  // int nEmptyNodeCount=0;
  // for (int nid = 1; nid < sTree.treeSize; nid++){
  //   if(sTree.coverListCounts[nid]!=0 || sTree.endListCounts[nid]!=0){
  //     nEmptyNodeCount++;
  //   }
  // }
  // cout<<"Non empty tree node count: "<<nEmptyNodeCount<<endl;
  // -----------------------------------------------------------------------

  // if(DEBUG_TIME) cp15 = high_resolution_clock::now();
  // int ncount=sTree.treeSize; //node 0 is invalid. True node count is treeSeze-2
  // int *iCount;
  // iCount=(int *)malloc((ncount)*sizeof(int));
  // int bType=0, cType=1;
  // // cSize=16205;

  // sTree.countIntersectionsMulticore
  //     (bPoly, cPoly, bSize, cSize,
  //     bType, cType,
  //     iCount);

  // if(DEBUG_TIME) cp16 = high_resolution_clock::now();
  // if(DEBUG) cout<<"Count intersections"<<endl;

  // thrust::exclusive_scan(iCount, iCount + ncount, iCount);
  // if(DEBUG) cout<<"Prefix sum: intersection count at nodes "<<iCount[ncount-1]<<endl;

  // int intersectCountDup=iCount[ncount-1];
  // int *bPolLineIdsDup, *cPolLineIdsDup, *bPolLineIdsSortedIndexDup, *intersectTypesDup;

  // bPolLineIdsDup=(int *)malloc((intersectCountDup)*sizeof(int));
  // cPolLineIdsDup=(int *)malloc((intersectCountDup)*sizeof(int));
  // bPolLineIdsSortedIndexDup=(int *)malloc((intersectCountDup)*sizeof(int));
  // intersectTypesDup=(int *)malloc((intersectCountDup)*sizeof(int));
  // if(DEBUG) cout<<"Clip start "<<intersectCountDup<<endl;

  // if(DEBUG_TIME) cp17 = high_resolution_clock::now();
  // sTree.polyClipArrayParallel
  //       (bPoly, cPoly, bSize, cSize,
  //       bType, cType,
  //       iCount, bPolLineIdsDup, cPolLineIdsDup,intersectTypesDup);
  // // sTree.polyClipArrayGPUParallel
  // //       (bPoly, cPoly, bSize, cSize,
  // //       bType, cType,
  // //       iCount, bPolLineIdsDup, cPolLineIdsDup,intersectTypesDup);
  // if(DEBUG) cout<<"Clip end"<<endl;
  //-----------------------------------------------------------------------
  //-----------------------------------------------------------------------
  // using STL vector
  int ncount = sTree.treeSize; // node 0 is invalid. True node count is treeSeze-2
  int bType=0; //bType is base polygon type. Always 0. bListId is the base poly id from the linked list

// ====================================================================================
  // **************** when b<c swap b and c for that case only. *************
// ====================================================================================

  // loop over all clipping layer polygons
  for (int bListId=0, cListId=0, cType=1; cType<=cPolyCount; cType++, bListId++, cListId++){ 
  // for (int bListId=0, cListId=0, cType=1; cType<=1; cType++, bListId++, cListId++){ 
    vector<int> bPolLineIdsDup, cPolLineIdsDup, intersectTypesDup;
    // base bType is always 0. cType>=1. 
    // When accessing cSize[] or cPoly[] use (cListId) as id of current clipping polygon
    int *bPolLineIdsSortedIndexDup;

    cSize=cSizes[cType]-cSizes[cListId]-1;
    // cSize=clSize;
    // use pointer vector

    sTree.saveIntersectionsIdsMulticore // for case 2
    // sTree.saveIntersectionsIdsMulticore2
        (bPoly, cPoly, bSize, cSize,
        bType, cType,
        bPolLineIdsDup, cPolLineIdsDup, intersectTypesDup);
    
    // printIntArray(&bPolLineIdsDup[0], bPolLineIdsDup.size(), "bPolLineIdsDup");
    // printIntArray(&cPolLineIdsDup[0], cPolLineIdsDup.size(), "cPolLineIdsDup");

    int intersectCountDup = bPolLineIdsDup.size();
    bPolLineIdsSortedIndexDup = (int *)malloc((intersectCountDup) * sizeof(int));
    #ifdef DG
    cout << "Intersecting IDs saved" << endl;
    #endif

    #ifdef TIME
    cp18 = high_resolution_clock::now();
    #endif
    radixsort(&bPolLineIdsDup[0], bPolLineIdsSortedIndexDup, intersectCountDup, bSize * 2);

    int *bPolLineIds, *cPolLineIds;
    int *intersectTypes, *bPolLineIdsSortedIndex, intersectCount;
    bPolLineIds = (int *)malloc((intersectCountDup) * sizeof(int));
    cPolLineIds = (int *)malloc((intersectCountDup) * sizeof(int));
    intersectTypes = (int *)malloc((intersectCountDup) * sizeof(int));
    bPolLineIdsSortedIndex = (int *)malloc((intersectCountDup) * sizeof(int));

    // remove intersection duplicates
    removeDuplicatesAtEdge(&bPolLineIdsDup[0], &cPolLineIdsDup[0], &intersectTypesDup[0], bPolLineIdsSortedIndexDup, intersectCountDup,
                          bPolLineIds, cPolLineIds, intersectTypes, bPolLineIdsSortedIndex, intersectCount);
    
    // printIntArray(intersectTypes, intersectCountDup, "intersectTypes");
    // removeDuplicatesAtEdge
    //   (bPolLineIdsDup, cPolLineIdsDup, intersectTypesDup, bPolLineIdsSortedIndexDup, intersectCountDup,
    //   bPolLineIds, cPolLineIds, intersectTypes, bPolLineIdsSortedIndex, intersectCount);
    //-----------------------------------------------------------------------
    //-----------------------------------------------------------------------
    // using linked list
    // int ncount=sTree.treeSize; //node 0 is invalid. True node count is treeSeze-2
    // list<int> bPolLineIdsDupList, cPolLineIdsDupList, intersectTypesDupList;
    // int bType=0, cType=1;
    // int *bPolLineIdsSortedIndexDup;

    // sTree.saveIntersectionsIdsMulticore   //for case 2
    // // sTree.saveIntersectionsIdsMulticore2
    //     (bPoly, cPoly, bSize, cSize,
    //     bType, cType,
    //     bPolLineIdsDupList, cPolLineIdsDupList, intersectTypesDupList);

    // int intersectCountDup=bPolLineIdsDupList.size();
    // bPolLineIdsSortedIndexDup=(int *)malloc((intersectCountDup)*sizeof(int));
    // if(DEBUG) cout<<"Intersecting IDs saved"<<endl;

    // // copy to arrays
    // int *bPolLineIdsDup, *cPolLineIdsDup, *intersectTypesDup;
    // bPolLineIdsDup=(int *)malloc((intersectCountDup)*sizeof(int));
    // cPolLineIdsDup=(int *)malloc((intersectCountDup)*sizeof(int));
    // intersectTypesDup=(int *)malloc((intersectCountDup)*sizeof(int));
    // auto itc=cPolLineIdsDupList->begin();
    // auto iti=intersectTypesDup->begin();
    // int i=0;
    // for(auto itb=bPolLineIdsDupList->begin(); itb!=bPolLineIdsDupList->end(); ++itb, ++itc, ++iti){
    //   bPolLineIdsDup[i]=itb;
    //   cPolLineIdsDup[i]=itc;
    //   intersectCountDup[i]=iti;
    //   ++i;
    // }

    // if(DEBUG_TIME) cp18 = high_resolution_clock::now();

    // radixsort(bPolLineIdsDup, bPolLineIdsSortedIndexDup, intersectCountDup, bSize*2);

    // int *bPolLineIds, *cPolLineIds;
    // int *intersectTypes, *bPolLineIdsSortedIndex, intersectCount;
    // bPolLineIds=(int *)malloc((intersectCountDup)*sizeof(int));
    // cPolLineIds=(int *)malloc((intersectCountDup)*sizeof(int));
    // intersectTypes=(int *)malloc((intersectCountDup)*sizeof(int));
    // bPolLineIdsSortedIndex=(int *)malloc((intersectCountDup)*sizeof(int));

    // // remove intersection duplicates
    // removeDuplicatesAtEdge
    //   (bPolLineIdsDup, cPolLineIdsDup, intersectTypesDup, bPolLineIdsSortedIndexDup, intersectCountDup,
    //   bPolLineIds, cPolLineIds, intersectTypes, bPolLineIdsSortedIndex, intersectCount);
    //-----------------------------------------------------------------------      
    #ifdef DG
    cout << intersectCountDup << " intersections found. " << intersectCountDup - intersectCount << " duplicate(s) removed. " << intersectCount << " remaining." << endl;
    #endif

    int *cPolLineIdsSortedIndex;
    cPolLineIdsSortedIndex = (int *)malloc((intersectCountDup) * sizeof(int));

    // sort cPolLineIds
    radixsort(cPolLineIds, cPolLineIdsSortedIndex, intersectCount, clSize * 2);
      
    #ifdef DG
    cout << "Sorted candidate edges by ID" << endl;
    #endif
    // printIntArray(bPolLineIds, intersectCountDup, "bPolLineIds");
    // printIntArray(cPolLineIds, intersectCountDup, "cPolLineIds");

    // get intersection counts at each edge for both polygons
    int *bAllIntersectCounts, *cAllIntersectCounts;
    int *bNonDegenIntersectCounts, *cNonDegenIntersectCounts;
    bAllIntersectCounts = (int *)malloc(((bSize + 1)) * sizeof(int));
    cAllIntersectCounts = (int *)malloc(((cSize + 1)) * sizeof(int));
    bNonDegenIntersectCounts = (int *)malloc(((bSize + 1)) * sizeof(int));
    cNonDegenIntersectCounts = (int *)malloc(((cSize + 1)) * sizeof(int));

    // default values for the extra location to use for prefix sum
    bAllIntersectCounts[bSize] = 0;
    cAllIntersectCounts[cSize] = 0;
    bNonDegenIntersectCounts[bSize] = 0;
    cNonDegenIntersectCounts[cSize] = 0;

    int *bPolLineUniqueIds;
    int *cPolLineUniqueIds;
    bPolLineUniqueIds = (int *)malloc((intersectCount) * sizeof(int));
    cPolLineUniqueIds = (int *)malloc((intersectCount) * sizeof(int));

    int uniqueIntersectEdgeCountB, uniqueIntersectEdgeCountC;
    countEdgeIntersections(bPolLineIds, bPolLineIdsSortedIndex, 0, intersectTypes, intersectCount,
                          'b', bAllIntersectCounts, bNonDegenIntersectCounts, bPolLineUniqueIds, bSize,
                          &uniqueIntersectEdgeCountB);

    // printIntArray(bAllIntersectCounts, (bSize+1), "bAllIntersectCounts");
    // printIntArray(bNonDegenIntersectCounts, (bSize+1), "bNonDegenIntersectCounts");
    // printIntArray(bPolLineUniqueIds, intersectCount, "bPolLineUniqueIds");

    #ifdef DG
    cout << "Base polygon edges counted" << endl;
    #endif
    countEdgeIntersections(cPolLineIds, cPolLineIdsSortedIndex, cSizes[cListId], intersectTypes, intersectCount,
                          'c', cAllIntersectCounts, cNonDegenIntersectCounts, cPolLineUniqueIds, cSize,
                          &uniqueIntersectEdgeCountC);

    // printIntArray(cAllIntersectCounts, (cSize+1), "cAllIntersectCounts");
    // printIntArray(cNonDegenIntersectCounts, (cSize+1), "cNonDegenIntersectCounts");
    // printIntArray(cPolLineUniqueIds, intersectCount, "cPolLineUniqueIds");

    #ifdef DG
    cout<<"Clipping polygon edges counted. b unique edge count="<<uniqueIntersectEdgeCountB<<". c unique edge count="<<uniqueIntersectEdgeCountC<<endl;
    #endif
    // prefix sum of the edge counts
    // #pragma acc routine seq
    // thrust::exclusive_scan(thrust::host, bAllIntersectCounts, bAllIntersectCounts+bSize+1, bAllIntersectCounts);
    thrust::exclusive_scan(bAllIntersectCounts, bAllIntersectCounts + bSize + 1, bAllIntersectCounts);
    thrust::exclusive_scan(cAllIntersectCounts, cAllIntersectCounts + cSize + 1, cAllIntersectCounts);
    thrust::exclusive_scan(bNonDegenIntersectCounts, bNonDegenIntersectCounts + bSize + 1, bNonDegenIntersectCounts);
    thrust::exclusive_scan(cNonDegenIntersectCounts, cNonDegenIntersectCounts + cSize + 1, cNonDegenIntersectCounts);

    // printIntArray(bNonDegenIntersectCounts, (bSize+1), "prefix sum bNonDegenIntersectCounts");
    // printIntArray(cAllIntersectCounts, (cSize+1), "prefix sum cAllIntersectCounts");
    // printIntArray(cNonDegenIntersectCounts, (cSize+1), "prefix sum cNonDegenIntersectCounts");

    #ifdef DG
      cout << "Edge prefix sums calculated" << endl;
    #endif

    int bNonDegenCount = bNonDegenIntersectCounts[bSize];
    int cNonDegenCount = cNonDegenIntersectCounts[cSize];
    int bIntersectedCount = bSize + bNonDegenCount;
    int cIntersectedCount = cSize + cNonDegenCount;

    #ifdef DG
    cout << "Poly B: non-degen count=" << bNonDegenCount << " total=" << bIntersectedCount;
    #endif
    #ifdef DG
    cout << " Poly C: non-degen count=" << cNonDegenCount << " total=" << cIntersectedCount << endl;
    #endif
    // cout<<bSize<<" "<<cSize<<" "<<bSize+cSize<<" "<<sTree.treeSize<<" "<<intersectCount<<" ";

    // save only intersections in new arrays
    // remove intersectType saving here. It is completed at polyClipArrayParallel()
    REAL *intersections;
    intersections = (REAL *)malloc((2 * intersectCount) * sizeof(REAL));

    int *intersectAlphaValuesB, *intersectAlphaValuesC;
    intersectAlphaValuesB = (int *)malloc((intersectCount) * sizeof(int));
    intersectAlphaValuesC = (int *)malloc((intersectCount) * sizeof(int));

    sTree.saveIntersectionsParallel(bPoly, cPoly,
                                    bPolLineIds, cPolLineIds, intersectCount,
                                    intersections, intersectTypes,
                                    intersectAlphaValuesB, intersectAlphaValuesC);
    
    // printDoubleArray(intersections, intersectCount*2, "intersections");
    // printIntArray(intersectTypes, intersectCount, "intersectTypes");
    // printIntArray(intersectAlphaValuesB, intersectCount, "intersectAlphaValuesB");
    // printIntArray(intersectAlphaValuesC, intersectCount, "intersectAlphaValuesC");
    
    #ifdef DG
    cout << "Intersections saved" << endl;
    #endif
    // printIntArray(intersectTypes, intersectCount, "intersectTypes");
    
    // printDoubleArray(intersections, intersectCount*2, "intersections");

    // sort saved edges in the common intersections array
    int *sortedIndiciesB, *sortedIndiciesC;
    sortedIndiciesB = (int *)malloc((intersectCount) * sizeof(int));
    sortedIndiciesC = (int *)malloc((intersectCount) * sizeof(int));

    sTree.sortIntersectionsAtEdgesParallel(intersectAlphaValuesB, bPolLineIdsSortedIndex, 0,
                                          bAllIntersectCounts, bNonDegenIntersectCounts, bPolLineUniqueIds, uniqueIntersectEdgeCountB, intersectCount,
                                          sortedIndiciesB);


    sTree.sortIntersectionsAtEdgesParallel(intersectAlphaValuesC, cPolLineIdsSortedIndex, cSizes[cListId],
                                          cAllIntersectCounts, cNonDegenIntersectCounts, cPolLineUniqueIds, uniqueIntersectEdgeCountC, intersectCount,
                                          sortedIndiciesC);      
    #ifdef DG
    cout << "Sort at edge completed" << endl;
    #endif
    // printDoubleArray(intersections, intersectCount*2, "intersections");
    // printIntArray(intersectTypes, intersectCount, "intersectTypes");
    // printIntArray(intersectAlphaValuesB, intersectCount, "intersectAlphaValuesB");
    // printIntArray(intersectAlphaValuesC, intersectCount, "intersectAlphaValuesC");
    // printDoubleArray(intersections, intersectCount*2, "intersections");

    REAL *intersectedB, *intersectedC;
    intersectedB = (REAL *)malloc((2 * bIntersectedCount) * sizeof(REAL));
    intersectedC = (REAL *)malloc((2 * cIntersectedCount) * sizeof(REAL));

    int *alphaValuesB, *alphaValuesC;
    int *neighborsB, *neighborsC;
    alphaValuesB = (int *)malloc((bIntersectedCount) * sizeof(int));
    alphaValuesC = (int *)malloc((cIntersectedCount) * sizeof(int));
    neighborsB = (int *)malloc((bIntersectedCount) * sizeof(int));
    neighborsC = (int *)malloc((cIntersectedCount) * sizeof(int));

    // data coping not rquired at this point. Existing CPU data update can be done with the available data
    sTree.copyIntersectionsParallel(bPoly, cPoly, bSize, cSize, cSizes[cListId],
                                    bPolLineIds, cPolLineIds, intersectCount,
                                    bPolLineIdsSortedIndex, cPolLineIdsSortedIndex,
                                    sortedIndiciesB, sortedIndiciesC,
                                    intersections, intersectTypes,
                                    bAllIntersectCounts, cAllIntersectCounts,
                                    bNonDegenIntersectCounts, cNonDegenIntersectCounts,
                                    intersectAlphaValuesB, intersectAlphaValuesC,
                                    intersectedB, intersectedC,
                                    neighborsB, neighborsC,
                                    alphaValuesB, alphaValuesC);

    // printDoubleArray(intersectedB, bIntersectedCount*2, "intersectedB");
    // printDoubleArray(intersectedC, cIntersectedCount*2, "intersectedC");
    // printIntArray(neighborsB, bIntersectedCount, "neighborsB");
    // printIntArray(neighborsC, cIntersectedCount, "neighborsC");
    // printIntArray(alphaValuesB, bIntersectedCount, "alphaValuesB");
    // printIntArray(alphaValuesC, cIntersectedCount, "alphaValuesC");
      
    #ifdef DG
   cout << "Intersections copied to source arrays" << endl;
    #endif
    #ifdef TIME
    cp3 = high_resolution_clock::now();
    #endif

    // initial labeling
    int *bInitLabels, *cInitLabels;
    bInitLabels = (int *)malloc((bIntersectedCount) * sizeof(int));
    cInitLabels = (int *)malloc((cIntersectedCount) * sizeof(int));

    // intialize c initLables
    for(int i=0; i<cIntersectedCount; ++i) cInitLabels[i]=0;

    sTree.calculateInitLabelParallel(bSize, bNonDegenIntersectCounts,
                                    intersectedB, intersectedC, alphaValuesB,
                                    neighborsB,
                                    bIntersectedCount, cIntersectedCount, bInitLabels, cInitLabels);
    // printIntArray(bInitLabels, bIntersectedCount, "bInitLabels");
    // printIntArray(cInitLabels, cIntersectedCount, "cInitLabels");

      
    #ifdef DG
    cout << "Initial labels calculated" << endl;
    #endif
    // copy data to the linked lists
    copyToLinkedList(intersectedB, intersectedC, bListId, cListId,
                    neighborsB, neighborsC,
                    alphaValuesB, alphaValuesC,
                    bIntersectedCount, cIntersectedCount,
                    bInitLabels, cInitLabels);

      
    #ifdef DG
    cout << "Copied to linked lists" << endl;
    #endif
      // print base and clipping polygons in their linked list form
    // next=bPolysList[0].root;
    // next=bPolysList[0].root;
    // int pid=0;
    // cout<<"Base poly "<<bIntersectedCount<<endl;
    // for(int i=0; i<bIntersectedCount; ++i){
    //   cout<<"pid="<<pid<<" x="<<next->p.x<<", y="<<next->p.y<<" label="<<next->label<<" enx="<<next->enex<<endl;
    //   next=next->next;
    // }
    // next=cPolysList[cListId].root;
    // cout<<"\nClip poly "<<cIntersectedCount<<endl;
    // for(int i=0; i<cIntersectedCount; ++i){
    //   cout<<"pid="<<pid<<" x="<<next->p.x<<", y="<<next->p.y<<" label="<<next->label<<" enx="<<next->enex<<endl;
    //   next=next->next;
    // }


    // PHASE: 3
    labelIntersections2(bListId, cListId);

    // next=bPolysList[0].root;
    // // int pid=0;
    // pid=0;
    // cout<<"Base poly "<<bIntersectedCount<<endl;
    // for(int i=0; i<bIntersectedCount; ++i){
    //   cout<<"pid="<<pid<<" x="<<next->p.x<<", y="<<next->p.y<<" label="<<next->label<<" enx="<<next->enex<<endl;
    //   next=next->next;
    // }
    // next=cPolysList[cListId].root;
    // cout<<"\nClip poly "<<cIntersectedCount<<endl;
    // for(int i=0; i<cIntersectedCount; ++i){
    //   cout<<"pid="<<pid<<" x="<<next->p.x<<", y="<<next->p.y<<" label="<<next->label<<" enx="<<next->enex<<endl;
    //   next=next->next;
    // }

    // PHASE: 4
    createResult2();
    run_time = omp_get_wtime() - start_time;
    #ifdef TIME
    cp4 = high_resolution_clock::now();
    #endif

    // if(DEBUG_TIME) end = high_resolution_clock::now();
    // post-processing
    cleanUpResult2();


    // write output polygon
    #ifdef SAVE
      cout << "R ";
      savePolygon(rPolyList, outputFile1+to_string(cType)+".txt");
    #endif
    // savePolygon2(rPolyList, outputFile);

    //clean arrays
    bPolLineIdsDup.clear();
    cPolLineIdsDup.clear();
    intersectTypesDup.clear();

    rPolyList.clear();
    

}
  // clean arrays

  free(intervals);
  free(bPoly);
  free(cPoly);


  bPolysList.clear();
  cPolysList.clear();

  // int lts=1;
  // int lts=log2(treeSize);
  // int lts = sTree.treeSize;
  // for (int i = 0; i < lts; ++i)
  // {
  //   // omp_destroy_lock(&lock[i] );
  //   omp_destroy_lock(&sTree.lock[i]);
  // }

  #ifdef DT_INFO
    // Calculating total time taken by the program.
    auto d1 = duration_cast<microseconds>(cp4 - cp1);
    auto d3 = duration_cast<microseconds>(cp3 - cp2);
    auto d4 = duration_cast<microseconds>(cp4 - cp3);
    auto d5 = duration_cast<microseconds>(cp9 - cp8);
    auto d7 = duration_cast<microseconds>(cp11 - cp10);
    auto d8 = duration_cast<microseconds>(cp13 - cp12);
    auto d9 = duration_cast<microseconds>(cp14 - cp1);

    cout << "PARALLEL TIMING\nAll time in microseconds\nTime: Segment tree skeleton construction: " << fixed << d9.count() << setprecision(10) << endl;
    cout << "Time: Copy linked list cover list to array: " << fixed << d8.count() << setprecision(10) << endl;
    cout << "Time: Cover list construction: " << fixed << d5.count() << setprecision(10) << endl;
    cout << "Time: End list construction: " << fixed << d7.count() << setprecision(10) << endl;
    cout << "Time: Candidate edge find and save: " << fixed << d3.count() << setprecision(10) << endl;
    cout << "Time: Labeling and stich results: " << fixed << d4.count() << setprecision(10) << endl;
    cout << "Time: Total : " << fixed << d1.count() << setprecision(10) << endl;

    printf("\n Execution time was %lf seconds\n ", run_time);
  #elif TIME
    // Calculating total time taken by the program.
    auto d1 = duration_cast<microseconds>(cp4 - cp1);
    auto d3 = duration_cast<microseconds>(cp3 - cp2);
    auto d4 = duration_cast<microseconds>(cp4 - cp3);
    auto d5 = duration_cast<microseconds>(cp9 - cp8);
    auto d7 = duration_cast<microseconds>(cp11 - cp10);
    auto d8 = duration_cast<microseconds>(cp13 - cp12);
    auto d9 = duration_cast<microseconds>(cp14 - cp1);

    cout << " " << fixed << d9.count() << setprecision(10) << " " << d7.count() << setprecision(10) << " " << d5.count() << setprecision(10) << " " << d3.count() << setprecision(10) << " "<< d4.count() << setprecision(10) << " " << d1.count() << setprecision(10) << endl;
  #endif
  // add this to destroy gpu data elements
  // #pragma acc exit data delete(coverListSizes, nodeIntervalArray, coverArray, this)
}

int main(int argc, char *argv[])
{
  #ifdef GPU
    cout<<"Hybrid algorithm"<<endl;
    hybrid_clipping(argc, argv); 
  #endif

  #ifdef CPU
    cout<<"Multi-core algorithm"<<endl;
    multicore_clipping(argc, argv);
  #endif

  return 0;
}