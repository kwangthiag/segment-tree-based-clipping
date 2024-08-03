#ifndef EDGE_H
#define EDGE_H

#include <string>

#include "data_types.h"

using namespace std;

typedef struct point{
  REAL x;
  REAL y;
}point;

typedef struct line{
  point p1;
  point p2;
}line;

class Edge{
  public:
    int id;
    int type;
    REAL start;
    REAL end;
    REAL first;
    
    Edge(){ }

    Edge(int i){
      id=i;
    }
    
    Edge(REAL aStart, REAL aEnd){
      start = aStart;
      first = aStart;
      end = aEnd;
    }
    
    Edge(int aId, int tp, REAL aStart, REAL aEnd){
      id = aId;
      type = tp;
      start = aStart;
      first = aStart;
      end = aEnd;
    }
    
    string toString(){
      string s = std::to_string(start);
      string e = std::to_string(end);
      string x = std::to_string(id);
      string t = std::to_string(type);
      return "{" + x + "}[" + t + "]("+ s + " " + e +")";
    }
    
    // returns true if this edge contains arg. e
    bool contains(Edge e){
      // if(start <= e.start && end >= e.end)
      if(start <= e.start && end >= e.end)
        return true;
      // else if(start == e.start && end == e.end)
      else
        return false;
    }
    
    bool contains(REAL x1){
      if(x1 >= start && x1 <= end)
        return true;
      else
        return false;
    }
    
    bool contains(REAL estart, REAL eend){
      // if(start <= e.start && end >= e.end)
      if(start <= estart && end >= eend)
        return true;
      // else if(start == e.start && end == e.end)
      else
        return false;
    }
    
    bool intersects(Edge e){
      // returns true if this edge contains arg. e or vice versa
      if((start <= e.start && end >= e.end) || (e.start <= start && e.end >= end))
        return true;
      // partial overlap
      if((start <= e.start && end >= e.start) ||  (e.end <= end && e.end >= start))
        return true;
      else
        return false;
    }

    bool intersects(REAL estart, REAL eend){
      // returns true if this edge contains arg. e or vice versa
      if((start <= estart && end >= eend) || (estart <= start && eend >= end))
        return true;
      // partial overlap
      if((start <= estart && end >= estart) ||  (eend <= end && eend >= start))
        return true;
      else
        return false;
    }

    void getData(int aId, int tp, REAL aStart, REAL aEnd){
      id = aId;
      type = tp;
      first = aStart;
      if(aStart<aEnd){
        start = aStart;
        end = aEnd;
      }else{
        end = aStart;
        start = aEnd;
      }
    }
};
#endif
