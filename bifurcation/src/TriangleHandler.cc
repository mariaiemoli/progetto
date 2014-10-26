#include "../include/TriangleHandler.h"

 
TriangleHandler::Intersection(FractureEnd const & gamma0, FractureEnd const & gamma1, FractureEnd const & gamma2, PointHandler const & intersectionPoint):
  fractures{gamma0,gamma1,gamma2}
{
  Vector2d tmp;
  for (unsigned int i=0;i<3;++i)
    {
		// get tangent at the end
		tmp(0) = intersectionPoint.x()-fractures[i].endPoint.x();
		tmp(1) = intersectionPoint.y()-fractures[i].endPoint.y();
		tmp.normalize();
		tangents[i]   = tmp;
		normals[i](0) = -tmp(1);
		normals[i](1) =  tmp(0);
    }
}

