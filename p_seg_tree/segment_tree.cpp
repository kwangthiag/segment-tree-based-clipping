#include "segment_tree/segment_tree.h"

//compile: g++ -std=c++0x -o seg segment_tree.cpp

// original source link:
//https://www.hackerearth.com/practice/data-structures/advanced-data-structures/segment-trees/tutorial/

void SegTree :: buildSkeleton(int node, int start, int end)
{
 	   if(start == end)
  	   {
    	    // Leaf node will have a single element
    	    treeNode->at(node)->interval = elementaryVec->at(start);
    	}
   		else
    	{
        	int mid = (start + end) / 2;
        	// Recurse on the left child
	        buildSkeleton(2*node, start, mid);
    	    // Recurse on the right child
        	buildSkeleton(2*node+1, mid+1, end);
	        // Internal node will have the sum of both of its children
    	    treeNode->at(node)->interval = seg_union(treeNode->at(2*node)->interval, treeNode->at(2*node+1)->interval);
    	}
}

void SegTree :: buildInternalEndLists()
{
	int firstLeafIndex = treeSize/2; 
    int lastLeafIndex = treeSize-1; 
    
    buildInternalEndLists(1, firstLeafIndex, lastLeafIndex);
}
  
void SegTree :: buildInternalEndLists(int node, int start, int end)
{
   if(start == end){
         // Leaf node will have a single element
         //treeNode->at(node)->interval = elementaryVec->at(start);
         //set<Edge,edge_compare> *endSet = treeNode->at(node)->endEdgeSet;
         //vector<Edge> *lEndList = new vector<Edge>( endSet->begin(), endSet->end());
         
   }else{
      int mid = (start + end) / 2;
      // Recurse on the left child
      buildInternalEndLists(2*node, start, mid);
      // Recurse on the right child
      buildInternalEndLists(2*node+1, mid+1, end);
         
      // Internal node will have the concatenation of both of its children
      set<Edge, edge_compare> *left = treeNode->at(2*node)->endEdgeSet;
      set<Edge, edge_compare> *right = treeNode->at(2*node + 1)->endEdgeSet;

      set<Edge, edge_compare> *root = treeNode->at(node)->endEdgeSet;
      root->insert(left->begin(), left->end());
      root->insert(right->begin(), right->end());    	    
   }
}

// displays interval and cover list for each node
void SegTree :: display()
{
	   float numLevels = log2(treeSize);
	   int nLevels = static_cast<int>(numLevels);
      if(DEBUG) cout<<"float_#Levels "<<numLevels<<" int_#Levels = "<<nLevels<<endl;
	   
	   for(int l = 1; l < treeSize; l++)
	   {
	      //cout<<l<<" "<<treeNode[l].interval.toString()<<endl;
	   if(DEBUG) cout<<" ["<<l<<"] Interval: "<<treeNode->at(l)->interval.toString()<<", size: "<<treeNode->at(l)->coverList->size()<<", ";//<<endl;
	      vector<Edge>* coverList = treeNode->at(l)->coverList;
         if(coverList->size() > 0)
         {
            for(Edge e : *coverList)
               if(DEBUG) cout<<" "<<e.toString()<<", ";
               
            if(DEBUG) cout<<endl;  
         }
    	   else
    	   if(DEBUG) cout<<endl;     
	   }
}

void SegTree :: insertAllToEdgeList(vector<Edge> *edges)
{
   insertToLeafLevelEdgeList(edges);
   buildInternalEndLists();
}

void SegTree :: insertToLeafLevelEdgeList(vector<Edge> *edges)
{
   map<REAL, pair<Node*, Node*> >* ptToNodeMap = preparePointsToNodeMap();
   
   if(ptToNodeMap == NULL)
   {
   if(DEBUG) cout<<"BUG: map empty"<<endl; 
      exit(1);
   }
   else
   {
      if(DEBUG) cout<<"Map Size "<<ptToNodeMap->size()<<endl;
   }
   
   std::map<REAL, pair<Node*, Node*> >::iterator it;
   		// cout<<"before for loop all edges "<<endl;
   for(Edge e : *edges)
   {
      // start point
              // cout<<"Start point "<<e.start<<endl;
      it = ptToNodeMap->find(e.start);
      if(it != ptToNodeMap->end())
      {
      	pair<Node*, Node*> nodePair = ptToNodeMap->at(e.start);
      
      	Node* after = nodePair.second;
      	//after->addToEndList(e);
      	after->addToEndEdgeSet(e);
      }
      
      // end point
      it = ptToNodeMap->find(e.end);
      if(it != ptToNodeMap->end())
      {
      	pair<Node*, Node*> nodePair = ptToNodeMap->at(e.end);
      	Node* before =nodePair.first;

      	//before->addToEndList(e);
      	before->addToEndEdgeSet(e);
      }
   }
   /*
   // convert set ds to vector ds
   int firstLeafIndex = treeSize/2; 
   int lastLeafIndex = treeSize-1;  
   
   for(int nodeIdx = firstLeafIndex; nodeIdx <= lastLeafIndex; nodeIdx++)
   {
      Node *leaf = treeNode->at(nodeIdx);
      //leaf->endEdgeSet
      set<Edge,edge_compare> *endSet = leaf->endEdgeSet;
      leaf->endList = new vector<Edge>( endSet->begin(), endSet->end());
   }
   */
}

void SegTree :: displayEndLists()
{
   //int firstLeafIndex = treeSize/2; 
   int lastLeafIndex = treeSize-1;  
   
   //for(int nodeIdx = firstLeafIndex; nodeIdx <= lastLeafIndex; nodeIdx++)
   for(int nodeIdx = 1; nodeIdx <= lastLeafIndex; nodeIdx++)
   {
      Node *leaf = treeNode->at(nodeIdx);
      if(DEBUG)  cout<<" [ "<<nodeIdx<<" ] "<<leaf->interval.toString()<<" : ";
      //for(Edge e: *leaf->endList)
      for(Edge e: *leaf->endEdgeSet)
      {
         if(DEBUG) cout<<e.toString()<<", ";
      }
      if(DEBUG) cout<<endl;
   }
}

map<REAL, pair<Node*, Node*> >* SegTree :: preparePointsToNodeMap()
{
   if(m_elementaryPoints == NULL && false == m_elementaryPoints->empty())
   {
     if(DEBUG) cout<< "Error: m_elementaryPoints in insertToLeafLevelEndList"<<endl;
     return NULL;
   }
   		
   // special handling for 1st point 
   REAL pt0 = m_elementaryPoints->at(0);
   		//cout<<"pt0 is "<<pt0<<endl;
   
   map<REAL, pair<Node*, Node*> > *ptToAdjNodeMap = new map<REAL, pair<Node*, Node*> >();

   int firstLeafIndex = treeSize/2;   
   Node *firstLeaf = treeNode->at(firstLeafIndex);
   pair <Node*, Node*> nodePair1 = make_pair(firstLeaf,firstLeaf);
   
   pair<REAL, pair<Node*, Node*> > kv1 = make_pair(pt0, nodePair1);
   ptToAdjNodeMap->insert(kv1);

   int lastLeafIndex = treeSize-1;
   if(DEBUG) cout<<endl<<" firstLeafIndex "<<firstLeafIndex<<" lastLeafIndex "<<lastLeafIndex<<" m_elementaryPoints->size() "<<m_elementaryPoints->size()<<endl;
   
   int leafIndex = firstLeafIndex;
   int i;
   for( i =1; i<m_elementaryPoints->size()-1; i++, leafIndex++)
   {
      Node *leftLeaf = treeNode->at(leafIndex);
      Node *rightLeaf = treeNode->at(leafIndex+1);
      
      pair <Node*, Node*> ithNodePair = make_pair(leftLeaf, rightLeaf);
      
      REAL ithPt = m_elementaryPoints->at(i);
      pair<REAL, pair<Node*, Node*> > kv = make_pair(ithPt, ithNodePair);
      ptToAdjNodeMap->insert(kv);
      
   }
   // special handling for the last point 
   //cout<<endl<<"i = "<<i<<endl;
   
   Node *lastLeaf = treeNode->at(lastLeafIndex);
   pair <Node*, Node*> lastNodePair = make_pair(lastLeaf,lastLeaf);

   REAL lastPt = m_elementaryPoints->at(m_elementaryPoints->size()-1);
   pair<REAL, pair<Node*, Node*> > lastKV = make_pair(lastPt, lastNodePair);
   ptToAdjNodeMap->insert(lastKV);
   
   return ptToAdjNodeMap;
}

void SegTree :: insertAll(vector<Edge> *edges)
{
   //int counter = 0;
   
   for(Edge e : *edges)
   {
      insert(e);
  	    //  cout<<"inserted "<<counter<<" "<<endl;
		//  counter++;
   } 
   
   // cout<<endl;
   //insert(Edge(1,4));
 //      insert(Edge(1,2));
}

void SegTree :: insertEdge(Edge x, Node *root, int nodeIdx)
{
    //cout<<"N# "<<nodeIdx<<endl;
    
	if(x.contains(root->interval)) 
    {
       //cout<<" Root  "<<nodeIdx<<" "<<endl;
       root->addToCoverList(x);
    }
    else 
    {
       int leftIdx = 2*nodeIdx;
       
       if(leftIdx < treeSize)
       {
         Node* leftChild = treeNode->at(leftIdx);
       
         if(x.intersects(leftChild->interval))
         {
            insertEdge(x, leftChild, leftIdx);
         }
       } 
       
       int rightIdx = 2*nodeIdx+1;

       if(rightIdx < treeSize)
       {
         Node* rightChild = treeNode->at(rightIdx);
       
         if(x.intersects(rightChild->interval))
         {
            insertEdge(x, rightChild, rightIdx);
         }
       }
     }  
}

void SegTree :: insert(Edge x)
{
    insertEdge(x, treeNode->at(1), 1);
}

vector<REAL>* SegTree :: getPointsFromEdges(vector<Edge> *intervals)
{
   set<REAL>* ptList = new set<REAL>();
   for(Edge e : *intervals)
   {
      ptList->insert(e.start);
      ptList->insert(e.end);
   }
   vector<REAL> *ptVect = new vector<REAL>(ptList->size());
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

vector<Edge>* SegTree :: getElementaryIntervals(vector<Edge> *intervals)
{
	vector<REAL>* pts = getPointsFromEdges(intervals);
	
	m_elementaryPoints = pts;
	
		// for(REAL coord : *pts)
//   		{
//        		cout<<coord<<" ";
//   		}
// 		cout<<endl;
	
	vector<Edge> *elementaryIntervals = new vector<Edge>();
	
	for(int i = 0; i < pts->size()-1; i++)
	{ 
	  REAL start = pts->at(i);
	  REAL end = pts->at(i+1);
	  Edge e(start, end);
	//   if(DEBUG) cout<<start<<"->"<<end<<" ";
	  elementaryIntervals->push_back(e);
	}
	
	// find upper bound on #intervals that is a 2 exponent
	int numIntervals = elementaryIntervals->size();
	int upperBound = Util::getNPowOf2(numIntervals);
	// if(DEBUG) cout<<endl<<"numIntervals "<<numIntervals<<" upperBound "<<upperBound<<endl;

	// normalize : Make the number of edges power of 2 for convenience 
	int i = pts->size() - 1;
	for(int j = i; j < upperBound; j++)
	{ 
	  REAL start = pts->at(i);
	  REAL end = pts->at(i);
	  Edge e(start, end);
	//   if(DEBUG) cout<<start<<"->"<<end<<" ";
	  elementaryIntervals->push_back(e);
	  
	  // this part is useful for endlist for leaves. We want to normalize the points 2^integer
	  m_elementaryPoints->push_back(start);
    }
	
	// if(DEBUG) cout<<endl;
	return elementaryIntervals;
}

vector<Edge>* SegTree :: querySegTree(REAL qx)
{ 
   // Input: the root of a (subtree of a) segment tree and query point qx.
   // Output: All intervals in the tree containing qx.
      vector<Edge> *edges = new vector<Edge>();
      
      recQueryTree(qx, treeNode->at(1), 1, edges);
      
      return edges;
}

void SegTree :: recQueryTree(REAL qx, Node *root, int nodeIdx, vector<Edge> *edges)
{
    //edges->push_back(root.coverList);
    string edge = root->interval.toString();
    
    //cout<<"cover list size "<<root->coverList->size()<<" "<<nodeIdx<<" "<<edge<<endl;
    
    if(root->coverList != NULL && root->coverList->size() > 0) {
       //cout<<"cover list size "<<root.coverList->size()<<endl;
       edges->insert(edges->end(), root->coverList->begin(), root->coverList->end());
    }
    
    // if root is not a leaf
    int leftIdx = 2*nodeIdx;
    
    if(leftIdx < treeSize)
    {
         Node *leftChild = treeNode->at(leftIdx);
       
         if(leftChild->interval.contains(qx))
         {
            recQueryTree(qx, leftChild, leftIdx, edges);
         }
    } 
    
    int rightIdx = 2*nodeIdx+1; 
    
    if(rightIdx < treeSize)
    {
         Node *rightChild = treeNode->at(rightIdx);
       
         if(rightChild->interval.contains(qx))
         {
            recQueryTree(qx, rightChild, rightIdx, edges);
         }
    }
    
}



/*
    	    vector<Edge> *left = treeNode->at(2*node)->endList;
    	    vector<Edge> *right = treeNode->at(2*node +1)->endList;
    	    
    	    treeNode->at(node)->endList->reserve(left->size() + right->size());
    	    treeNode->at(node)->endList->insert(treeNode->at(node)->endList->end(), left->begin(), left->end());
    	    treeNode->at(node)->endList->insert(treeNode->at(node)->endList->end(), right->begin(), right->end());
    	    */