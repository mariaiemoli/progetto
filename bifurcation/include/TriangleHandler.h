/**
 * TriangleHandler.h
 * 
 * Classe che costruisce il triangolo 2d dell'size_typeersezione
 * 
 */

#ifndef __TRIANGLEHANDLER_H__
#define __TRIANGLEHANDLER_H__

#include "TriangleData.h"


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
		

	private:
	  //! The end point of the fracture
	  PointData  M_endPoint;
	  // The thickness of the fracture
	  scalar_type M_thickness;
};

//! Class to represent the intersection
class Intersection
{
	public:
	  //! Fractures must be set and given in clockwise order.
	  Intersection(FractureEnd const & gamma0, FractureEnd const & gamma1, FractureEnd const & gamma2, PointData const & intersectionPoint);

	  TriangleData const & computeIntersectionTriangle();
  
	  TriangleData const & intersectionTriangle()const 
		  {
			  return intersectionTriangle_;
		  }

	private:
	  FractureEnd fractures[3];
	  PointData intersection_;
	  Vector2d tangents[3];
	  Vector2d normals[3];
  
	  TriangleData intersectionTriangle_;
  
};

#endif