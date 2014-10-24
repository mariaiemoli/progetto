/**
 * TriangleHandler.h
 * 
 * Classe che costruisce il triangolo 2d dell'intersezione
 * 
 */

#ifndef __TRIANGLEHANDLER_H__
#define __TRIANGLEHANDLER_H__1

#include <iosfwd>
#include <Eigen/Dense>
#include "Core.h"
#include <assert.h>
#include "UsefulFunctions.h"
#include "PointHandler.h"

typedef Eigen::Vector2d Vector;

const int ndim=2;


class TriangleHandler {
	
public:
    static int const myDim=2;
    static const int numVertices=3;
    static const int numSides=3;

    Triangle( const FracturePtrContainer_Type& M_fractures ); //Constructs an empty triangle
    
    Triangle(Point&,Point&,Point&); //Points are given (by reference)
    
    Triangle(const Triangle&);
    
    Triangle & operator=(const Triangle&);
    
    // We get the points by operator [] (defined in-class for inlining)
    void setPoint(int i, Point const &  p);
    
    // Can be used ONLY if empty()==false
    Point const & operator[](int i) const {return M_points[i];}
    
    // Can be used ONLY if empty()==false
    Point & operator[](int i){return M_points[i];}
    
    double measure() const; // Triangle area
    
    Point& edgePoint(int edgenum,int endnum); // The point on an edge
    
    Point const & edgePoint(int edgenum,int endnum) const; // The const version
    
    //!Get baricenter
    Point baricenter()const;
    
    //!Get edge baricenter
    Point edgeBaricenter(int edgeNum) const;
    
    //!Get vector connecting edge baricenter with baricenter
    Vector c(int edgeNum) const;
    
    // Outward normal to edge times edge length;
    Vector unscaledNormal(int edgeNum)const;
    
    static int edge(int edgeNum, int endNum); // The edge numbering
    
    friend std::ostream & operator <<(std::ostream &, Triangle const &);

    
	
private:
    PointHandler  M_points[numVertices];
    static int const M_edge[numSides][2];

}; // TriangleHandler


typedef TriangleHandler TriangleHandler_Type;
typedef std::vector < TriangleHandler_Type > TriangleHandlerContainer_Type;
typedef boost::shared_ptr < TriangleHandler_Type > TriangleHandlerPtr_Type;
typedef std::vector < TriangleHandlerPtr_Type> TriangleHandlerPtrContainer_Type;


#endif