/**
 * TriangleHandler.h
 * 
 * Classe che costruisce il triangolo 2d dell'size_typeersezione
 * 
 */

#ifndef __TRIANGLEHANDLER_H__
#define __TRIANGLEHANDLER_H__

#include "TriangleData.h"
#include "FractureHandler.h"


//! Class that represents the end of a fracture
class FractureEnd
{
	public:
		//Costruttore
	    FractureEnd( PointData& , scalar_type ); //Da tre punti presi per refernza
		
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
	  //! The end point of the fracture
	  PointData  M_endPoint;
	  // The thickness of the fracture
	  scalar_type M_thickness;
};

typedef FractureEnd FractureEnd_Type;
typedef std::vector < FractureEnd_Type > FractureEndContainer_Type;
typedef boost::shared_ptr < FractureEnd_Type > FractureEndPtr_Type;
typedef std::vector < FractureEndPtr_Type> FractureEndPtrContainer_Type;

//! Class to represent the intersection
class Intersection
{
	public:
		
	  //Costruttori
	  //! Fractures must be set and given in clockwise order.
	  Intersection(FractureEnd const & gamma0, FractureEnd const & gamma1, FractureEnd const & gamma2, PointData const & intersectionPoint);
	  
	  //Troviamo il triangolo di intersezione partendo dal vettore delle fratture che si intersecano
	  Intersection( const FracturePtrContainer_Type& M_fractures );
	  
	  //Metodi
	  TriangleData const & computeIntersectionTriangle();
  
	  TriangleData const & intersectionTriangle()const 
		  {
			  return intersectionTriangle_;
		  }

	private:
	  FractureEndContainer_Type fractures;
	  PointData intersection_;
	  Vector2d tangents;
	  Vector2d normals;
  
	  TriangleData intersectionTriangle_;
  
};

typedef Intersection Intersection_Type;
typedef std::vector < Intersection_Type > IntersectionContainer_Type;
typedef boost::shared_ptr < Intersection_Type > IntersectionPtr_Type;
typedef std::vector < IntersectionPtr_Type> IntersectionPtrContainer_Type;

#endif