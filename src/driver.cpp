#include "../p_seg_tree/segment_tree/segment_tree.h"

int main()
{
  int numEdges = 7;
  //int numEdges = 5;
  
  SegTree sTree(numEdges);
  
  //  Edge edgeArr[7] = {Edge(1,2), Edge(3,4), Edge(2,5), Edge(5,6), Edge(7,8), Edge(9,10), Edge(11,12)};
   // Edge edgeArr[7] = {Edge("s1", 1,2), Edge("s2", 3,4), Edge("s3", 2,5), Edge("s4", 5,6), Edge("s5",7,8), Edge("s6", 9,10), Edge("s7", 11,12)};
   Edge edgeArr[7] = {Edge("s1", 1,2), Edge("s2", 3,8), Edge("s3", 1,5), Edge("s4", 5,6), Edge("s5",3,8), Edge("s6", 2,10), Edge("s7", 11,12)};
  
  //Edge edgeArr[5] = {Edge(1,4), Edge(3,5), Edge(6,8), Edge(7,8), Edge(3,7)};
  // Edge edgeArr[5] = {Edge("s1",1,4), Edge("s2",3,5), Edge("s3", 6,8), Edge("s4",7,8), Edge("s5",3,7)};
  // This DS holds the original edges
  vector<Edge> intervals;
  for(int i = 0; i<numEdges; i++)
  {
     intervals.push_back(edgeArr[i]);
  }
  
  vector<Edge> *elementaryIntervals = sTree.getElementaryIntervals(&intervals);
 
  // The code assumes that the elementary intervals are power of 2.
  sTree.init(elementaryIntervals);
  
     cout<<"Elementary Intervals "<<elementaryIntervals->size() <<endl;
  
  sTree.buildSkeleton(1, 0, elementaryIntervals->size() - 1);
  
     cout<<"Build finished "<<endl;

  sTree.display();
  
  sTree.insertAll(&intervals);

  sTree.display();

  cout<<endl;
  
  vector<Edge> *answers = sTree.querySegTree(4);
  
  cout<<"answers"<<endl;
  for(Edge e : *answers)
     cout<<e.toString()<<" "; 
  
  cout<<"\n-----"<<endl;

  sTree.insertAllToEdgeList(&intervals);
  
  sTree.displayEndLists();
     cout<<"Done "<<endl;
  return 0;
}