
#include "../include/TriangleHandler.h"

/**************************************************************************/
/*  TriangleHandler.h													  */
/*  Libreria che costruisce il triangolo 2d dell'intersezione             */
/**************************************************************************/


FractureEnd::FractureEnd(PointData& a , scalar_type t)
{
	M_endPoint = a;

	M_thickness = t;
}//costruttore FractureEnd
 
Intersection::Intersection(FractureEnd const & gamma0, FractureEnd const & gamma1, FractureEnd const & gamma2, 
								PointData const & intersectionPoint) : M_intersection (intersectionPoint)
{
	//QUESTO é DA RISCRIVERE NELLA FORMA : 
	M_fractures.clear();
  	M_fractures.push_back( gamma0 ); 
  	M_fractures.push_back( gamma1 );
  	M_fractures.push_back( gamma2 );
	
	Vector2d ttmp;
	Vector2d ntmp;
		
	for (size_type i=0;i<3;++i)
	{
		// get tangent at the end
		ttmp(0) = intersectionPoint.x()-M_fractures[i].getPoint ().x();

		ttmp(1) = intersectionPoint.y()-M_fractures[i].getPoint ().y();
		ttmp.normalize();

		M_tangents [ i ] = ttmp;
		ntmp(0) = -ttmp(1); 
		ntmp(1) = ttmp(0); 

		M_normals [ i ] = ntmp;
	}
	
	M_intersectionTriangle = this->computeIntersectionTriangle();

	
}// costruttore intersezione


void Intersection::setIntersection( FracturePtrContainer_Type& M_FracturesSet )
{	
	if( M_FracturesSet.size() == 2 )
	{
		setIntersectionWith2Fractures( M_FracturesSet );
	}
	else
	{
		setIntersectionWith3Fractures( M_FracturesSet );
	}
	
    
}// costruttore intersezione

void Intersection::setIntersectionWith2Fractures( FracturePtrContainer_Type& M_FracturesSet )
{	
	const size_type nbDof0 =  M_FracturesSet [ 0 ]-> getMeshFEMPressure().nb_basic_dof();
	const size_type nbDof1 =  M_FracturesSet [ 1 ]-> getMeshFEMPressure().nb_basic_dof();
					std::cout << "Node " <<  M_FracturesSet[ 0 ]->getDOFBifurcation() << std::endl;
	size_type DOFsplit = M_FracturesSet[ 0 ]->getDOFBifurcation()[ 0 ];
						std::cout << "FIne Node " << std::endl;
	
	base_node node0(2);
	base_node node1(2);
	base_node node2(2);
	base_node nodeI(2);
	
	//truccetto xk y_map vuole i base_node
	base_node tmp0( 1 );
	base_node tmp1( 1 );
	base_node tmp2( 1 );
	base_node tmp3( 1 );
	tmp0[ 0 ]= DOFsplit/( nbDof0 -1 );
	tmp1[ 0 ]= (DOFsplit + 1 )/( nbDof0 -1 );
	tmp2[ 0 ]= 0.;
	tmp3[ 0 ]= 1.;
	

	
	node0[ 0 ] = M_FracturesSet [ 0 ]-> getMeshFEMPressure().point_of_basic_dof( DOFsplit )[ 0 ];
	node1[ 0 ] = M_FracturesSet [ 0 ]-> getMeshFEMPressure().point_of_basic_dof( DOFsplit + 1 )[ 0 ];
	node2[ 0 ] = M_FracturesSet [ 1 ]-> getMeshFEMPressure().point_of_basic_dof( 0 )[ 0 ];
	
	nodeI[ 0 ] = M_FracturesSet [ 1 ]-> getMeshFEMPressure().point_of_basic_dof( 0 )[ 0 ];
	
	node0[ 1 ] = M_FracturesSet [ 0 ]-> getLevelSet()->getData()->y_map( tmp0 );
	node1[ 1 ] = M_FracturesSet [ 0 ]-> getLevelSet()->getData()->y_map( tmp1);
	node2[ 1 ] = M_FracturesSet [ 1 ]-> getLevelSet()->getData()->y_map( tmp2 );
	
	nodeI[ 1 ] = M_FracturesSet [ 1 ]-> getLevelSet()->getData()->y_map( tmp2 );
	

	
	//Estremi liberi delle fratture +  punto di intersezione
	PointData p0( node0[0], node0[1] );
	PointData p1( node1[0], node1[1] );
	PointData p2( node2[0], node2[1] );
	
	PointData pi( nodeI[0], nodeI[1] );
	
	//Spessori delle mie fratture
	scalar_type t0,t1;
	t0 = M_FracturesSet [ 0 ] -> getData().getThickness();
	t1 = M_FracturesSet [ 1 ] -> getData().getThickness();

	FractureEnd f0( p0, t0 );
    FractureEnd f1( p1, t0 );
	FractureEnd f2( p2, t1 );
	
    Intersection tmp(f0,f1,f2,pi);
    
    this-> M_fractures = tmp.M_fractures;
    this-> M_intersection = tmp.M_intersection;
    
    this-> M_tangents [ 0 ] = tmp.M_tangents [ 0 ];
    this-> M_tangents [ 1 ] = tmp.M_tangents [ 1 ];
    this-> M_tangents [ 2 ] = tmp.M_tangents [ 2 ];
    
    this-> M_normals [ 0 ] = tmp.M_normals [ 0 ];
    this-> M_normals [ 1 ] = tmp.M_normals [ 1 ];
    this-> M_normals [ 2 ] = tmp.M_normals [ 2 ];
    
    this-> M_intersectionTriangle = tmp.M_intersectionTriangle; 
	
	std::cout << tmp.M_intersectionTriangle << std::endl;
	
	return;
}

void Intersection::setIntersectionWith3Fractures( FracturePtrContainer_Type& M_FracturesSet )
{
	// ordino le fratture in senso orario
	sortFractures ( M_FracturesSet );

	
	const size_type nbDof0 =  M_FracturesSet [ 0 ]-> getMeshFEMPressure().nb_basic_dof();
	const size_type nbDof1 =  M_FracturesSet [ 1 ]-> getMeshFEMPressure().nb_basic_dof();
	const size_type nbDof2 =  M_FracturesSet [ 2 ]-> getMeshFEMPressure().nb_basic_dof();
	
	base_node node0(2);
	base_node node1(2);
	base_node node2(2);
	base_node nodeI(2);
	
	//truccetto xk y_map vuole i base_node
	base_node tmp0( 1 );
	base_node tmp1( 1 );
	tmp0[ 0 ]= 0.;
	tmp1[ 0 ]= 1.;
	
	node0[ 0 ] = M_FracturesSet [ 0 ]-> getMeshFEMPressure().point_of_basic_dof( 0 )[ 0 ];
	node1[ 0 ] = M_FracturesSet [ 1 ]-> getMeshFEMPressure().point_of_basic_dof( 0 )[ 0 ];
	node2[ 0 ] = M_FracturesSet [ 2 ]-> getMeshFEMPressure().point_of_basic_dof( 0 )[ 0 ];
	nodeI[ 0 ] = M_FracturesSet [ 0 ]-> getMeshFEMPressure().point_of_basic_dof( nbDof0-1 )[ 0 ];
	
	node0[ 1 ] = M_FracturesSet [ 0 ]-> getLevelSet()->getData()->y_map( tmp0 );
	node1[ 1 ] = M_FracturesSet [ 1 ]-> getLevelSet()->getData()->y_map( tmp0 );
	node2[ 1 ] = M_FracturesSet [ 2 ]-> getLevelSet()->getData()->y_map( tmp0 );
	nodeI[ 1 ] = M_FracturesSet [ 0 ]-> getLevelSet()->getData()->y_map( tmp1 );
	
	if ( gmm::abs( M_FracturesSet [ 1 ]-> getLevelSet()->getData()->ylevelSetFunction ( node0 ) ) < 1.0E-2 )
	{
		node0 [ 0 ] = M_FracturesSet [ 0 ]-> getMeshFEMPressure().point_of_basic_dof( nbDof0-1 )[ 0 ];
		node0 [ 0 ] = M_FracturesSet [ 0 ]-> getMeshFEMPressure().point_of_basic_dof( nbDof0-1 )[ 0 ];
		nodeI [ 0 ] = M_FracturesSet [ 0 ]-> getMeshFEMPressure().point_of_basic_dof( 0 )[ 0 ];
		node0 [ 1 ] = M_FracturesSet [ 0 ]-> getLevelSet()->getData()->y_map( tmp1 );
		nodeI [ 1 ] = M_FracturesSet [ 0 ]-> getLevelSet()->getData()->y_map( tmp0 );
	}

	if ( gmm::abs( M_FracturesSet [ 0 ]-> getLevelSet()->getData()->ylevelSetFunction ( node1 ) ) < 1.0E-2 )
	{
		
		node1 [ 0 ] = M_FracturesSet [ 1 ]-> getMeshFEMPressure().point_of_basic_dof( nbDof1-1 )[ 0 ];
		node1 [ 1 ] = M_FracturesSet [ 1 ]-> getLevelSet()->getData()->y_map( tmp1 );
	}

	if ( gmm::abs( M_FracturesSet [ 1 ]-> getLevelSet()->getData()->ylevelSetFunction ( node2 ) ) < 1.0E-2 )
	{
		node2 [ 0 ] = M_FracturesSet [ 2 ]-> getMeshFEMPressure().point_of_basic_dof( nbDof2-1 )[ 0 ];
		node2 [ 1 ] = M_FracturesSet [ 2 ]-> getLevelSet()->getData()->y_map( tmp1 );
	}

	//Estremi liberi delle fratture +  punto di intersezione
	PointData p0( node0[0], node0[1] );
	PointData p1( node1[0], node1[1] );
	PointData p2( node2[0], node2[1] );
	PointData pi( nodeI[0], nodeI[1] );
	
	//Spessori delle mie fratture
	scalar_type t0,t1,t2;
	t0 = M_FracturesSet [ 0 ] -> getData().getThickness();
	t1 = M_FracturesSet [ 1 ] -> getData().getThickness();
	t2 = M_FracturesSet [ 2 ] -> getData().getThickness();
	

	FractureEnd f0( p0, t0 );
    FractureEnd f1( p1, t1 );
	FractureEnd f2( p2, t2 );
	
    Intersection tmp(f0,f1,f2,pi);
    
    this-> M_fractures = tmp.M_fractures;
    this-> M_intersection = tmp.M_intersection;
    
    this-> M_tangents [ 0 ] = tmp.M_tangents [ 0 ];
    this-> M_tangents [ 1 ] = tmp.M_tangents [ 1 ];
    this-> M_tangents [ 2 ] = tmp.M_tangents [ 2 ];
    
    this-> M_normals [ 0 ] = tmp.M_normals [ 0 ];
    this-> M_normals [ 1 ] = tmp.M_normals [ 1 ];
    this-> M_normals [ 2 ] = tmp.M_normals [ 2 ];
    
    this-> M_intersectionTriangle = tmp.M_intersectionTriangle; 
	
	std::cout << tmp.M_intersectionTriangle << std::endl;
	
	return;
}


TriangleData const & Intersection::computeIntersectionTriangle()
{
	PointData point;
	const scalar_type tol=1.e-5;
	scalar_type s(0.);
	
	for (size_type i=0; i<3;++i)
	{
		size_type j = (i +1 ) % 3;
		scalar_type ninj = M_normals[i].dot(M_normals[j] );
		scalar_type nitj = M_normals[i].dot(M_tangents[j]);
		
		if(std::fabs(nitj)<tol)
		{
			point = M_intersection + 
					PointData ( M_normals [ i ] ( 0 ), M_normals [ i ] ( 1 ) ) * ( 0.5 * ( M_fractures [ j ].getThickness() + M_fractures [ i ].getThickness() ));
		}
		else
		{
			// parametric coordinate
			s    = 0.5 * ( M_fractures [ j ].getThickness() * ninj + M_fractures [ i ].getThickness() )/nitj;
		
			// the ith point is Pji in the note
			point = M_intersection + ( PointData ( M_tangents [ j ] ( 0 ), M_tangents [ j ] ( 1 ) ) * s )- 
									( PointData ( M_normals [ j ] ( 0 ), M_normals [ j ] ( 1 ) ) * ( 0.5 * M_fractures [ j ].getThickness() ) );
		}
	
		M_intersectionTriangle.setPoint(i,point);
		
	}
	
	return M_intersectionTriangle;

}// computeIntersectionTriangle


void Intersection::sortFractures ( FracturePtrContainer_Type& M_FracturesSet )
{

	FractureHandlerPtr_Type f0 = M_FracturesSet [ 0 ];
	FractureHandlerPtr_Type f1 = M_FracturesSet [ 1 ];
	FractureHandlerPtr_Type f2 = M_FracturesSet [ 1 ];
	
	if ( isPos ( f0, f1) )
	{	
		if ( !isPos ( f0, f2 ) )
			return;
		else if ( isPos ( f1, f2 ) )
			return;
		else
		{	
			// devo invertire la due e la uno
			FractureHandlerPtr_Type tmp = f1;
			M_FracturesSet [ 1 ] = M_FracturesSet [ 2 ];
			M_FracturesSet [ 2 ] = tmp;
			
			return;
		}
	}
		
	else 
	{	
		if ( isPos ( f0, f2 ) )
		{
			// devo invertire la 0 e la 1
			FractureHandlerPtr_Type tmp = f1;
			M_FracturesSet [ 1 ] = M_FracturesSet [ 0 ];
			M_FracturesSet [ 0 ] = tmp;
			
			return;

		}
		else
		{
			if ( isPos ( f1, f2 ) )
			{	
				// devo invertire la 2 e la 1
				FractureHandlerPtr_Type tmp = f1;
				M_FracturesSet [ 1 ] = M_FracturesSet [ 2 ];
				M_FracturesSet [ 2 ] = tmp;
				
				return;

			}

		}
	}
	
	return;
}// sortFractures



bool Intersection::isPos ( FractureHandlerPtr_Type& fracture, FractureHandlerPtr_Type& otherfracture)
{
	base_node n0(2);
	n0 [ 0 ] = 0;
	n0 [ 1 ] = otherfracture->getLevelSet()->getData()->y_map ( n0 );
	n0 [ 0 ] = otherfracture->getLevelSet()->getData()->x_map ( n0 );

	base_node n1(2);
	n1 [ 0 ] = 1;
	n1 [ 1 ] = otherfracture->getLevelSet()->getData()->y_map ( n0 );
	n1 [ 0 ] = otherfracture->getLevelSet()->getData()->x_map ( n0 );

	if ( fracture->getLevelSet()->getData()->ylevelSetFunction ( n0 ) >= 0 && fracture->getLevelSet()->getData()->ylevelSetFunction ( n1 ) >= 0 )
		return true;
	else 
		return false;
	
}// isPos
