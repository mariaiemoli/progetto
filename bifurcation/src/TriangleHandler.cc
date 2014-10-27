#include "../include/TriangleHandler.h"


TriangleHandler::FractureEnd(PointData& a , scalar_type t)
{
  M_point = a;
  M_endPoint = t;
}//costruttore dati i punti
 
TriangleHandler::Intersection(FractureEnd const & gamma0, FractureEnd const & gamma1, FractureEnd const & gamma2, 
								PointData const & intersectionPoint) //: fractures{gamma0,gamma1,gamma2} 
{
	//QUESTO é DA RISCRIVERE NELLA FORMA : 
	fracture.clear();
  	fracture.push_back( gamma0 ); 
  	fracture.push_back( gamma1 );
  	fracture.push_back( gamma2 );
	
	intersection_ = intersectionPoint;
	
	Vector2d ttmp;
	Vector2d ntmp;
	
	tangents.clear();
	normals.clear();
	
	for (unsigned int i=0;i<3;++i)
	{
		// get tangent at the end
		ttmp(0) = intersectionPoint.x()-fractures[i].getPoint.x();
		ttmp(1) = intersectionPoint.y()-fractures[i].getPoint.y();
		tmp.normalize();
		tangents.push_back( ttmp );
		ntmp(0) = -ttmp(1); 
		ntmp(1) = ttmp(0); 
		normals.push_back( ntmp );
	}
	
	intersectionTriangle_ = intersection.computeIntersectionTriangle();
}


TriangleHandler::Intersection( const FracturePtrContainer_Type& M_fractures )
{
	M_point.clear();
	
	assert ( M_fractures.size() == 3 );
	
	LevelSetHandlerPtr_Type l0 = M_fractures [ 0 ]-> getLevelSet ();
	LevelSetHandlerPtr_Type l1 = M_fractures [ 1 ]-> getLevelSet ();
	LevelSetHandlerPtr_Type l2 = M_fractures [ 2 ]-> getLevelSet ();
	
	const size_type nbDof0 = l0 ->getLevelSet()-> get_mesh_fem().nb_basic_dof();
	const size_type nbDof1 = l1 ->getLevelSet()-> get_mesh_fem().nb_basic_dof();
	const size_type nbDof2 = l2 ->getLevelSet()-> get_mesh_fem().nb_basic_dof();
	
	base_node node0 = l0-> getLevelSet()-> get_mesh_fem().point_of_basic_dof( 0 );
	base_node node1 = l1-> getLevelSet()-> get_mesh_fem().point_of_basic_dof( 0 );
	base_node node2 = l2-> getLevelSet()-> get_mesh_fem().point_of_basic_dof( 0 );
	base_node nodeI = l0-> getLevelSet()-> get_mesh_fem().point_of_basic_dof( nbDof0-1 );

	if ( l1->getData()->levelSetFunction ( node0 ) == 0 )
	{
		node0 = l0-> getLevelSet()-> get_mesh_fem().point_of_basic_dof( nbDof0-1 );
		nodeI = l0-> getLevelSet()-> get_mesh_fem().point_of_basic_dof( 0 );
	}

	if ( l0->getData()->levelSetFunction ( node1 ) == 0 )
	{
		node1 = l1-> getLevelSet()-> get_mesh_fem().point_of_basic_dof( nbDof1-1 );
	}

	if ( l1->getData()->levelSetFunction ( node2 ) == 0 )
	{
		node2 = l2->  getLevelSet()-> get_mesh_fem().point_of_basic_dof( nbDof2-1 );
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
	
    Intersection intersection(f0,f1,f2,pI);
    
}

