#ifndef NODE_H
#define NODE_H

#include <vector>
#include <set>
#include <list>
#include "edge.h"

// struct edge_compare {
//     bool operator() (const Edge& lhs, const Edge& rhs) const 
//     {
//         if(lhs.id == rhs.id)
//            return true;
//         else
//            return false;
//     }
// };

struct edge_compare {
   bool operator() (const Edge& lhs, const Edge& rhs) const 
   {
      if(lhs.type==rhs.type && lhs.id==rhs.id)
         return false;
      else
         return lhs.id<=rhs.id;
         // return true;
   }
};

class Node
{
   public:
   Edge interval;
   list<Edge> *coverList;  //Chanage to list 
   vector<Edge> *endList;
   set<Edge, edge_compare> *endEdgeSet;
   
   Node()
   {
      coverList = new list<Edge>();
      endList = new vector<Edge>();
      endEdgeSet = new set<Edge, edge_compare>();
   }
   
   void addToCoverList(Edge e)
   {
     coverList->push_back(e);
   }
   
   void addToEndList(Edge e)
   {
     endList->push_back(e);
   }
   
   void addToEndEdgeSet(Edge e)
   {
     endEdgeSet->insert(e);
   }
};

#endif
