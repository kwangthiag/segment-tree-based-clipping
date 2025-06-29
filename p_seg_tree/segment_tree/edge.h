// #ifndef EDGE_H
// #define EDGE_H

// #include <string>

// #include "data_types.h"

// using namespace std;

// typedef struct point{
//   REAL x;
//   REAL y;
// }point;

// typedef struct line{
//   point p1;
//   point p2;
// }line;

// class Edge{
//   public:
//     int id;
//     int type;
//     REAL start;
//     REAL end;
//     REAL first;
    
//     Edge(){ }

//     Edge(int i){
//       id=i;
//     }
    
//     Edge(REAL aStart, REAL aEnd){
//       start = aStart;
//       first = aStart;
//       end = aEnd;
//     }
    
//     Edge(int aId, int tp, REAL aStart, REAL aEnd){
//       id = aId;
//       type = tp;
//       start = aStart;
//       first = aStart;
//       end = aEnd;
//     }
    
//     string toString(){
//       string s = std::to_string(start);
//       string e = std::to_string(end);
//       string x = std::to_string(id);
//       string t = std::to_string(type);
//       return "{" + x + "}[" + t + "]("+ s + " " + e +")";
//     }
    
//     // returns true if this edge contains arg. e
//     bool contains(Edge e) const {
//       // if(start <= e.start && end >= e.end)
//       if(start <= e.start && end >= e.end)
//         return true;
//       // else if(start == e.start && end == e.end)
//       else
//         return false;
//     }
    
//     bool contains(REAL x1) const {
//       if(x1 >= start && x1 <= end)
//         return true;
//       else
//         return false;
//     }
    
//     bool contains(REAL estart, REAL eend) const {
//       // if(start <= e.start && end >= e.end)
//       if(start <= estart && end >= eend)
//         return true;
//       // else if(start == e.start && end == e.end)
//       else
//         return false;
//     }
    
//     bool intersects(const Edge e) const {
//       // returns true if this edge contains arg. e or vice versa
//       if((start <= e.start && end >= e.end) || (e.start <= start && e.end >= end))
//         return true;
//       // partial overlap
//       if((start <= e.start && end >= e.start) ||  (e.end <= end && e.end >= start))
//         return true;
//       else
//         return false;
//     }

//     bool intersects(REAL estart, REAL eend) const {
//       // returns true if this edge contains arg. e or vice versa
//       if((start <= estart && end >= eend) || (estart <= start && eend >= end))
//         return true;
//       // partial overlap
//       if((start <= estart && end >= estart) ||  (eend <= end && eend >= start))
//         return true;
//       else
//         return false;
//     }

//     void getData(int aId, int tp, REAL aStart, REAL aEnd){
//       id = aId;
//       type = tp;
//       first = aStart;
//       if(aStart<aEnd){
//         start = aStart;
//         end = aEnd;
//       }else{
//         end = aStart;
//         start = aEnd;
//       }
//     }
// };
// #endif


#ifndef EDGE_H
#define EDGE_H

#include <string>
#include "data_types.h"

using namespace std;

// Forward declaration of the Edge class
class Edge;

typedef struct point {
    REAL x;
    REAL y;
} point;

typedef struct line {
    point p1;
    point p2;
} line;

class Edge {
public:
    int id;
    int type;
    REAL start;
    REAL end;
    REAL y; // Assuming this member exists for vertical ordering

    Edge() : id(0), type(0), start(0), end(0), y(0) {}

    Edge(int i) : id(i), type(0), start(0), end(0), y(0) {}

    Edge(REAL aStart, REAL aEnd) : id(0), type(0), start(aStart), end(aEnd), y(aStart) {}

    Edge(int aId, int tp, REAL aStart, REAL aEnd) : id(aId), type(tp), start(aStart), end(aEnd), y(aStart) {}

    string toString() const { // Added const
        string s = std::to_string(start);
        string e = std::to_string(end);
        string x = std::to_string(id);
        string t = std::to_string(type);
        return "{" + x + "}[" + t + "](" + s + " " + e + ")";
    }

    // CORRECTED: Added 'const' to all member functions that do not modify the object.
    bool contains(const Edge& e) const {
        return (start <= e.start && end >= e.end);
    }
    
    bool contains(REAL x1) const {
        return (x1 >= start && x1 <= end);
    }
    
    bool contains(REAL estart, REAL eend) const {
        return (start <= estart && end >= eend);
    }
    
    bool intersects(const Edge& e) const {
        if ((start <= e.start && end >= e.end) || (e.start <= start && e.end >= end)) return true;
        if ((start <= e.start && end >= e.start) || (e.end <= end && e.end >= start)) return true;
        return false;
    }

    bool intersects(REAL estart, REAL eend) const {
        if ((start <= estart && end >= eend) || (estart <= start && eend >= end)) return true;
        if ((start <= estart && end >= estart) || (eend <= end && eend >= start)) return true;
        return false;
    }

    void getData(int aId, int tp, REAL aStart, REAL aEnd) {
        id = aId;
        type = tp;
        y = aStart; // Assuming the y-coordinate is the start
        if (aStart < aEnd) {
            start = aStart;
            end = aEnd;
        } else {
            end = aStart;
            start = aEnd;
        }
    }
};
#endif
