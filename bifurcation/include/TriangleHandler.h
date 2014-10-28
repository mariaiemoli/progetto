/**
 * TriangleHandler.h
 * 
 * Classe che costruisce il triangolo 2d dell'intersezione
 * 
 */

#ifndef __TRIANGLEHANDLER_H__
#define __TRIANGLEHANDLER_H__

#include "TriangleData.h"
#include "FractureHandler.h"


/**
 * Classe che rappresenta la fine di una frattura
 */ 
class FractureEnd
{
	public:
		
		//Costruttore
	    FractureEnd( PointData& , scalar_type ); 
		
	    inline PointData getPoint ( ) const
	    {
	        return M_endPoint;
	    }
		
	    inline scalar_type getThickness ( ) const
	    {
	        return M_thickness;
	    }
		
		//Operatori
		FractureEnd & operator=(const FractureEnd& t)
		{
		  	if (this !=&t)
		   	{
		  		M_endPoint = t.getPoint();
		  		M_thickness = t.getThickness();
		    }
	
			return *this;
		}

	private:
	  
		// Punto finale della frattura
		PointData  M_endPoint;
		// Spessore della frattura
		scalar_type M_thickness;
		
};

typedef FractureEnd FractureEnd_Type;
typedef std::vector < FractureEnd_Type > FractureEndContainer_Type;
typedef boost::shared_ptr < FractureEnd_Type > FractureEndPtr_Type;
typedef std::vector < FractureEndPtr_Type> FractureEndPtrContainer_Type;

/**
 * Classe che rappresenta l'intersezione
 */
class Intersection
{
	public:
		
		
		//Costruttori
		/**
		* Le fratture devono essere date in senso orario
		*/ 
		Intersection(FractureEnd const & gamma0, FractureEnd const & gamma1, FractureEnd const & gamma2, PointData const & intersectionPoint);
		
		Intersection()
		{};
		
		//Troviamo il triangolo di intersezione partendo dal vettore delle fratture che si intersecano
		void setIntersection( const FracturePtrContainer_Type& M_FracturesSet );
		
		
		//Metodi
		TriangleData const & computeIntersectionTriangle();
		
		void setTriangle ( const TriangleData& triangle)
		{
			M_intersectionTriangle = triangle;
		}
		
		TriangleData const & intersectionTriangle()const 
		{
			return M_intersectionTriangle;
		}
		
		PointData getPointIntersection( ) const
		{
			return M_intersection;
		}
		
		
		Intersection & operator =(const Intersection & Int)
		{
			this-> M_fractures = Int.M_fractures;
			this-> M_intersection = Int.M_intersection;
			
			this-> M_tangents [ 0 ] = Int.M_tangents [ 0 ];
			this-> M_tangents [ 1 ] = Int.M_tangents [ 1 ];
			this-> M_tangents [ 2 ] = Int.M_tangents [ 2 ];
			
			this-> M_normals [ 0 ] = Int.M_normals [ 0 ];
			this-> M_normals [ 1 ] = Int.M_normals [ 1 ];
			this-> M_normals [ 2 ] = Int.M_normals [ 2 ];
			
			this-> M_intersectionTriangle = Int.M_intersectionTriangle; 
			
		}
	  

	private:
		  
		  FractureEndContainer_Type M_fractures;
		  PointData M_intersection;
		  Vector2d M_tangents [ 3 ];
		  Vector2d M_normals [ 3 ];
	  
		  TriangleData M_intersectionTriangle;
  
};

typedef Intersection Intersection_Type;
typedef std::vector < Intersection_Type > IntersectionContainer_Type;
typedef boost::shared_ptr < Intersection_Type > IntersectionPtr_Type;
typedef std::vector < IntersectionPtr_Type> IntersectionPtrContainer_Type;

#endif