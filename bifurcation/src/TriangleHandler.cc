#include "../include/TriangleHandler.h"


FractureEnd::FractureEnd(PointData& a , scalar_type t)
{
	M_endPoint = a;

	M_thickness = t;
}//costruttore FractureEnd
 
Intersection::Intersection(FractureEnd const & gamma0, FractureEnd const & gamma1, FractureEnd const & gamma2, 
								PointData const & intersectionPoint) : intersection_ (intersectionPoint)
{
	//QUESTO é DA RISCRIVERE NELLA FORMA : 
	fractures.clear();
  	fractures.push_back( gamma0 ); 
  	fractures.push_back( gamma1 );
  	fractures.push_back( gamma2 );
	
	Vector2d ttmp;
	Vector2d ntmp;
		
	for (size_type i=0;i<3;++i)
	{
		// get tangent at the end
		ttmp(0) = intersectionPoint.x()-fractures[i].getPoint ().x();

		ttmp(1) = intersectionPoint.y()-fractures[i].getPoint ().y();
		ttmp.normalize();

		tangents [ i ] = ttmp;
		ntmp(0) = -ttmp(1); 
		ntmp(1) = ttmp(0); 

		normals [ i ] = ntmp;
	}
	
	intersectionTriangle_ = this->computeIntersectionTriangle();
	
}// costruttore intersezione


void Intersection::setIntersection( const FracturePtrContainer_Type& M_fractures )
{	
	assert ( M_fractures.size() == 3 );
	
	const size_type nbDof0 =  M_fractures [ 0 ]-> getMeshFEMVelocity().nb_basic_dof();
	const size_type nbDof1 =  M_fractures [ 1 ]-> getMeshFEMVelocity().nb_basic_dof();
	const size_type nbDof2 =  M_fractures [ 2 ]-> getMeshFEMVelocity().nb_basic_dof();
	
	base_node node0(2);
	base_node node1(2);
	base_node node2(2);
	base_node nodeI(2);
	
	//truccetto xk y_map vuole i base_node
	base_node tmp0( 1 );
	base_node tmp1( 1 );
	tmp0[ 0 ]= 0.;
	tmp1[ 0 ]= 1.;
	
	node0[ 0 ] = M_fractures [ 0 ]-> getMeshFEMVelocity().point_of_basic_dof( 0 )[ 0 ];
	node1[ 0 ] = M_fractures [ 1 ]-> getMeshFEMVelocity().point_of_basic_dof( 0 )[ 0 ];
	node2[ 0 ] = M_fractures [ 2 ]-> getMeshFEMVelocity().point_of_basic_dof( 0 )[ 0 ];
	nodeI[ 0 ] = M_fractures [ 0 ]-> getMeshFEMVelocity().point_of_basic_dof( 0 )[ 0 ];
	
	node0[ 1 ] = M_fractures [ 0 ]-> getLevelSet()->getData()->y_map( tmp0 );
	node1[ 1 ] = M_fractures [ 1 ]-> getLevelSet()->getData()->y_map( tmp0 );
	node2[ 1 ] = M_fractures [ 2 ]-> getLevelSet()->getData()->y_map( tmp0 );
	nodeI[ 1 ] = M_fractures [ 0 ]-> getLevelSet()->getData()->y_map( tmp1 );
	
	if ( M_fractures [ 1 ]-> getLevelSet()->getData()->ylevelSetFunction ( node0 ) == 0 )
	{
		node0 [ 0 ] = M_fractures [ 0 ]-> getMeshFEMVelocity().point_of_basic_dof( nbDof0-1 )[ 0 ];
		nodeI [ 0 ] = M_fractures [ 0 ]-> getMeshFEMVelocity().point_of_basic_dof( 0 )[ 0 ];
		node0 [ 1 ] = M_fractures [ 0 ]-> getLevelSet()->getData()->y_map( tmp1 );
		nodeI [ 1 ] = M_fractures [ 0 ]-> getLevelSet()->getData()->y_map( tmp0 );
	}

	if ( M_fractures [ 0 ]-> getLevelSet()->getData()->ylevelSetFunction ( node1 ) == 0 )
	{
		node1 [ 0 ] = M_fractures [ 1 ]-> getMeshFEMVelocity().point_of_basic_dof( nbDof1-1 )[ 0 ];
		node1 [ 1 ] = M_fractures [ 1 ]-> getLevelSet()->getData()->y_map( tmp1 );
	}

	if ( M_fractures [ 1 ]-> getLevelSet()->getData()->ylevelSetFunction ( node2 ) == 0)
	{
		node2 [ 0 ] = M_fractures [ 2 ]-> getMeshFEMVelocity().point_of_basic_dof( nbDof2-1 )[ 0 ];
		node2 [ 1 ] = M_fractures [ 2 ]-> getLevelSet()->getData()->y_map( tmp1 );
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
	
    Intersection tmp(f0,f1,f2,pi);
    
    this-> fractures = tmp.fractures;
    this-> intersection_ = tmp.intersection_;
    
    this-> tangents [ 0 ] = tmp.tangents [ 0 ];
    this-> tangents [ 1 ] = tmp.tangents [ 1 ];
    this-> tangents [ 2 ] = tmp.tangents [ 2 ];
    
    this-> normals [ 0 ] = tmp.normals [ 0 ];
    this-> normals [ 1 ] = tmp.normals [ 1 ];
    this-> normals [ 2 ] = tmp.normals [ 2 ];
    
    this-> intersectionTriangle_ = tmp.intersectionTriangle_; 
	
	std::cout << tmp.intersectionTriangle_ << std::endl;
    
}// costruttore intersezione


TriangleData const & Intersection::computeIntersectionTriangle()
{
	PointData point;
	const scalar_type tol=1.e-5;
	scalar_type s(0.);
	
	for (size_type i=0; i<3;++i)
	{
		size_type j = (i +1 ) % 3;
		scalar_type ninj = normals[i].dot(normals[j] );
		scalar_type nitj = normals[i].dot(tangents[j]);
		
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

}// computeIntersectionTriangle

