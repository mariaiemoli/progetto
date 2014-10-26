/**
 * TriangleHandler.h
 * 
 * Classe che costruisce il triangolo 2d dell'size_typeersezione
 * 
 */

#ifndef __TRIANGLEHANDLER_H__
#define __TRIANGLEHANDLER_H__


//! Class that represents the end of a fracture
struct FractureEnd
{
  //! The end point of the fracture
  PointHandler  endPoint;
  // The thickness of the fracture
  scalar_type thickness;
};

//! Class to represent the intersection
class Intersection
{
public:
  //! Fractures must be set and given in clockwise order.
  Intersection(FractureEnd const & gamma0, FractureEnd const & gamma1, FractureEnd const & gamma2, PointHandler const & intersectionPoint);

  TriangleData const & computeIntersectionTriangle();
  
  TriangleData const & intersectionTriangle()const {return intersectionTriangle_;}

private:
  FractureEnd fractures[3];
  PointHandler intersection_;
  Vector2d tangents[3];
  Vector2d normals[3];
  
  TriangleData intersectionTriangle_;
  
};

#endif