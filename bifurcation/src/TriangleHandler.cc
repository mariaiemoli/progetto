#include "../include/TriangleHandler.h"


FractureEnd::FractureEnd(PointData& a , scalar_type t)
{
	M_endPoint = a;

	M_thickness = t;
}//costruttore FractureEnd
 
Intersection::Intersection(FractureEnd const & gamma0, FractureEnd const & gamma1, FractureEnd const & gamma2, 
								PointData const & intersectionPoint) //: fractures{gamma0,gamma1,gamma2} 
{
	//QUESTO é DA RISCRIVERE NELLA FORMA : 
	fractures.clear();
  	fractures.push_back( gamma0 ); 
  	fractures.push_back( gamma1 );
  	fractures.push_back( gamma2 );
	
	intersection_ = intersectionPoint;
	
	Vector2d ttmp;
	Vector2d ntmp;
	
//	tangents.clear();
//	normals.clear();
	
	for (unsigned int i=0;i<3;++i)
	{
		// get tangent at the end
		ttmp(0) = intersectionPoint.x()-fractures[i].getPoint ().x();
		ttmp(1) = intersectionPoint.y()-fractures[i].getPoint ().y();
		ttmp.normalize();
//		tangents.push_back( ttmp );
		tangents [ i ] = ttmp;
		ntmp(0) = -ttmp(1); 
		ntmp(1) = ttmp(0); 
//		normals.push_back( ntmp );
		normals [ i ] = ntmp;
	}
	
	intersectionTriangle_ = this->computeIntersectionTriangle();
}// costruttore intersezione


Intersection::Intersection( const FracturePtrContainer_Type& M_fractures )
{
	//M_point.clear();
	
	assert ( M_fractures.size() == 3 );
	
	LevelSetHandlerPtr_Type l0 = M_fractures [ 0 ]-> getLevelSet ();
	LevelSetHandlerPtr_Type l1 = M_fractures [ 1 ]-> getLevelSet ();
	LevelSetHandlerPtr_Type l2 = M_fractures [ 2 ]-> getLevelSet ();
	
	const size_type nbDof0 = l0 ->getLevelSet().get_mesh_fem().nb_basic_dof();
	const size_type nbDof1 = l1 ->getLevelSet().get_mesh_fem().nb_basic_dof();
	const size_type nbDof2 = l2 ->getLevelSet().get_mesh_fem().nb_basic_dof();
	
	base_node node0 = l0-> getLevelSet().get_mesh_fem().point_of_basic_dof( 0 );
	base_node node1 = l1-> getLevelSet().get_mesh_fem().point_of_basic_dof( 0 );
	base_node node2 = l2-> getLevelSet().get_mesh_fem().point_of_basic_dof( 0 );
	base_node nodeI = l0-> getLevelSet().get_mesh_fem().point_of_basic_dof( nbDof0-1 );

	if ( l1->getData()->levelSetFunction ( node0 ) == 0 )
	{
		node0 = l0-> getLevelSet().get_mesh_fem().point_of_basic_dof( nbDof0-1 );
		nodeI = l0-> getLevelSet().get_mesh_fem().point_of_basic_dof( 0 );
	}

	if ( l0->getData()->levelSetFunction ( node1 ) == 0 )
	{
		node1 = l1-> getLevelSet().get_mesh_fem().point_of_basic_dof( nbDof1-1 );
	}

	if ( l1->getData()->levelSetFunction ( node2 ) == 0 )
	{
		node2 = l2->  getLevelSet().get_mesh_fem().point_of_basic_dof( nbDof2-1 );
	}

	//Estremi liberi delle fratture +  punto di intersezione
	PointData p0( node0[0], node0[1] );
	PointData p1( node1[0], node1[1] );
	PointData p2( node2[0], node2[1] );
	PointData pi( nodeI[0], nodeI[1] );
	
	//Spessori delle mie fratture
	scalar_type t0,t1,t2;
	t0 = M_fractures [ 0 ] -> getData().getThickness();
	t1 = M_fractures [ 1 ] -> getData().getThickness();
	t2 = M_fractures [ 2 ] -> getData().getThickness();
	
    FractureEnd f0( p0, t0 );
    FractureEnd f1( p1, t1 );
    FractureEnd f2( p2, t2 );
	
    Intersection intersection(f0,f1,f2,pi);
    
}// costruttore intersezione


TriangleData const & Intersection::computeIntersectionTriangle()
{
	PointData point;
	const double tol=1.e-5;
	double s(0.);
	
	for (unsigned int i=0; i<3;++i)
	{
		unsigned int j = (i +1 ) % 3;
		double ninj = normals[i].dot(normals[j] );
		double nitj = normals[i].dot(tangents[j]);
		
		if(std::fabs(nitj)<tol)
		{
			point = intersection_ + 
					PointData ( normals [ i ] ( 0 ), normals [ i ] ( 1 ) ) * ( 0.5 * ( fractures [ j ].getThickness() + fractures [ i ].getThickness() ));
		}
		else
		{
			// parametric coordinate
			s    = 0.5 * ( fractures [ j ].getThickness() * ninj + fractures [ i ].getThickness() )/nitj;
		
			// the ith point is Pji in the note
			point = intersection_ + ( PointData ( tangents [ j ] ( 0 ), tangents [ j ] ( 1 ) ) * s )- 
									( PointData ( normals [ j ] ( 0 ), normals [ j ] ( 1 ) ) * ( 0.5 * fractures [ j ].getThickness() ) );
		}
	
		intersectionTriangle_.setPoint(i,point);
	}
	
	return intersectionTriangle_;

}

Vector3d pressureCoeff(const Matrix3d &T)
{
  double factor = (1.0/T.sum());
  Vector3d tmp;
  tmp(0)=factor*(T.col(0).sum());
  tmp(1)=factor*(T.col(1).sum());
  tmp(2)=factor*(T.col(2).sum());
  return tmp;
}