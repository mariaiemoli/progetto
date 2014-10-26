#include "../include/TriangleHandler.h"


TriangleHandler::FractureEnd(PointData& a , scalar_type t)
{
  M_point = a;
  M_endPoint = t;
}//costruttore dati i punti
 
TriangleHandler::Intersection(FractureEnd const & gamma0, FractureEnd const & gamma1, FractureEnd const & gamma2, 
								PointData const & intersectionPoint //: fractures{gamma0,gamma1,gamma2} 
{
	this->fracture[0] = gamma0;
	this->fracture[1] = gamma1;
	this->fracture[2] = gamma2;
	
	Vector2d tmp;
	
	for (unsigned int i=0;i<3;++i)
	{
		// get tangent at the end
		tmp(0) = intersectionPoint.x()-fractures[i].getPoint.x();
		tmp(1) = intersectionPoint.y()-fractures[i].getPoint.y();
		tmp.normalize();
		tangents[i]   = tmp;
		normals[i](0) = -tmp(1);
		normals[i](1) =  tmp(0);
	}
}

