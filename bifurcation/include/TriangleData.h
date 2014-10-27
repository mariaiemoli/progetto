/**
 * TriangleData.h
 * 
 * Classe che contiene i dati principali per costrire il triangolo
 * 
 */

#ifndef __TRIANGLEDATA_H__
#define __TRIANGLEDATA_H__

#include <eigen3/Eigen/Dense>

#include <assert.h>
#include "PointData.h"
#include "StringUtility.h"

class TriangleData {
	
	public:
		static int const myDim=2;
		static const int numVertices=3;
		static const int numSides=3;
		
		//Costruttori
		TriangleData(); //vuoto
		
		TriangleData(PointData&,PointData&,PointData&); //Da tre punti presi per refernza
		
		TriangleData(const TriangleData&); // Da un altro triangolo == Copia
		
		//Mi permette di accedere al punto i-esimo del triangolo
		PointData& getPoint( size_type i ) 
		{
			assert ( i < M_point.size() );
			
			return M_point[i];
		}
		
		// We get the point by operator [] (defined in-class for inlining)
		void setPoint(size_type i, PointData const &  p);
		
		//misura l'area del tringolo
		scalar_type measure() const; 
		
		PointData& edgePoint(size_type edgenum,size_type endnum); // The posize_type on an edge
		
		PointData const & edgePoint(size_type edgenum,size_type endnum) const; // The const version
		
		//!Get baricenter
		PointData baricenter()const;
		
		static size_type edge(size_type edgeNum, size_type endNum); // The edge numbering
		
		//!Get edge baricenter
		PointData edgeBaricenter(size_type edgeNum) const;
		
		//!Get vector connecting edge baricenter with baricenter
		Vector2d c(size_type edgeNum) const;
		
		// Outward normal to edge times edge length;
		Vector2d unscaledNormal(size_type edgeNum)const;
		
		//Operatori
		TriangleData& operator=(const TriangleData&);
		
		friend std::ostream & operator <<(std::ostream &, TriangleData const &);
		
		// Can be used ONLY if empty()==false
		PointData const & operator[](size_type i) const 
		{
			return M_point[i];
		}
		
		// Can be used ONLY if empty()==false
		PointData & operator[](size_type i)
		{
			return M_point[i];
		}
		
		
	private:
		PointDataContainer_Type  M_point;
		static size_type const M_edge[numSides][2];
		
		//static sizeVectorContainer_Type const M_edge = {{0, 1}, {1, 2}, {2, 0}};
		

}; // TriangleData


typedef TriangleData TriangleData_Type;
typedef std::vector < TriangleData_Type > TriangleDataContainer_Type;
typedef boost::shared_ptr < TriangleData_Type > TriangleDataPtr_Type;
typedef std::vector < TriangleDataPtr_Type> TriangleDataPtrContainer_Type;


#endif