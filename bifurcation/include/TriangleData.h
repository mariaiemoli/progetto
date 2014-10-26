/**
 * TriangleData.h
 * 
 * Classe che contiene i dati principali per costrire il triangolo
 * 
 */

#ifndef __TRIANGLEDATA_H__
#define __TRIANGLEDATA_H__


#include <iosfwd>
#include <Eigen/Dense>

#include <assert.h>
#include "UsefulFunctions.h"
#include "PointHandler.h"
#include "FractureHandler.h"
#include "TriangleData.h"

class TriangleData {
	
public:
    static size_type const myDim=2;
    static const size_type numVertices=3;
    static const size_type numSides=3;

    TriangleData( const FracturePtrContainer_Type& M_fractures ); //Constructs an empty triangle
    
    TriangleData(PointHandler&,PointHandler&,PointHandler&); //PointHandlers are given (by reference)
    
    TriangleData(const TriangleData&);
    
    TriangleData& operator=(const TriangleData&);
    
    // We get the point by operator [] (defined in-class for inlining)
    void setPoint(size_type i, PointHandler const &  p);
    
    // Can be used ONLY if empty()==false
    PointHandler const & operator[](size_type i) const {return M_point[i];}
    
    // Can be used ONLY if empty()==false
    PointHandler & operator[](size_type i){return M_point[i];}
    
    scalar_type measure() const; // Triangle area
    
    PointHandler& edgePoint(size_type edgenum,size_type endnum); // The posize_type on an edge
    
    PointHandler const & edgePoint(size_type edgenum,size_type endnum) const; // The const version
    
    //!Get baricenter
    PointHandler baricenter()const;
    
    //!Get edge baricenter
    PointHandler edgeBaricenter(size_type edgeNum) const;
    
    //!Get vector connecting edge baricenter with baricenter
    Vector2d c(size_type edgeNum) const;
    
    // Outward normal to edge times edge length;
    Vector2d unscaledNormal(size_type edgeNum)const;
    
    static size_type edge(size_type edgeNum, size_type endNum); // The edge numbering
    
    friend std::ostream & operator <<(std::ostream &, TriangleData const &);
	
	
private:
    PointHandler  M_point[numVertices];
    static size_type const M_edge[numSides][2];
	

}; // TriangleData


typedef TriangleData TriangleData_Type;
typedef std::vector < TriangleData_Type > TriangleDataContainer_Type;
typedef boost::shared_ptr < TriangleData_Type > TriangleDataPtr_Type;
typedef std::vector < TriangleDataPtr_Type> TriangleDataPtrContainer_Type;


#endif