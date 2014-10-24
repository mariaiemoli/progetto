#include <iostream>
#include "../include/TriangleHandler.h"


// Definition of the static variable edge

const size_type TriangleHandler::M_edge[3][2] = {{0, 1}, {1, 2}, {2, 0}};

//COSTRUTTORI	
TriangleHandler::TriangleHandler()
	{}// costruttore vuoro


TriangleHandler::TriangleHandler(PointHandler & a, PointHandler & b, PointHandler & c)
{
  M_point[0] = a;
  M_point[1] = b;
  M_point[2] = c;
}//costruttore dati i punti

TriangleHandler::TriangleHandler( const FracturePtrContainer_Type& M_fractures )
{
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
	PointHandler p0( node0[0], node0[1] );
	PointHandler p1( node1[0], node1[1] );
	PointHandler p2( node2[0], node2[1] );
	PointHandler pi( nodeI[0], nodeI[1] );
	
	//Spessori delle mie fratture
	scalar_type t0,t1,t2;
	t0 = M_fractures [ 0 ] -> getData().getThickness();
	t1 = M_fractures [ 1 ] -> getData().getThickness();
	t2 = M_fractures [ 2 ] -> getData().getThickness();
	
    FractureEnd f0{ p0, t0 };
    FractureEnd f1{ p1, t1 };
    FractureEnd f2{ p2, t2 };
	
    Intersection intersection(f0,f1,f2,pI);
    
	this->intersection.computeIntersectionTriangle();
    
}//Costruttore date le fratture

void TriangleHandler::setPoint(size_type i, PointHandler const & p)
{
M_point[i] = p;
}

PointHandler & TriangleHandler::edgePoint(size_type edgenum, size_type endnum)
{
return M_point[M_edge[edgenum][endnum]];
}

const PointHandler & TriangleHandler::edgePoint(size_type edgenum, size_type endnum) const
{
return M_point[M_edge[edgenum][endnum]];
}

TriangleHandler::TriangleHandler(const TriangleHandler & t)
{
M_point[0]=t.M_point[0];
M_point[1]=t.M_point[1];
	M_point[2]=t.M_point[2];
}

TriangleHandler & TriangleHandler::operator =(const TriangleHandler & t)
{
if (this !=&t){
  M_point[0]=t.M_point[0];
  M_point[1]=t.M_point[1];
  M_point[2]=t.M_point[2];
}
return *this;
}

size_type TriangleHandler::edge(size_type i, size_type j)
{
return M_edge[i][j];
}

scalar_type TriangleHandler::measure() const
{
const TriangleHandler & t = *this;
return -0.5 *
  (t[1][0] * (t[2][1] - t[0][1]) +
   t[2][0] * (t[0][1] - t[1][1]) +
   t[0][0] * (t[1][1] - t[2][1]));
}

PointHandler TriangleHandler::baricenter()const
{
PointHandler tmp(M_point[0]);
tmp += M_point[1];
tmp += M_point[2];
return tmp * (1./3.0);
}

PointHandler TriangleHandler::edgeBaricenter(size_type edgeNum)const
{
PointHandler tmp(M_point[edge(edgeNum,0)]+M_point[edge(edgeNum,1)]);
return tmp *0.5;
}

Vector2d TriangleHandler::c(size_type edgeNum) const
{
PointHandler baric = this->baricenter();
PointHandler eBaric= this->edgeBaricenter(edgeNum);
return eBaric.asVector()-baric.asVector();
}

Vector2d TriangleHandler::unscaledNormal(size_type edgeNum) const
{
PointHandler const & a=this->edgePoint(edgeNum,0);
PointHandler const & b=this->edgePoint(edgeNum,1);
return Vector2d(a.y()-b.y(),b.x()-a.x());
}

std::ostream & operator <<(std::ostream & stream, TriangleHandler const & t)
{
stream<<" ************* TRIANGLE  POINTS ********"<<std::endl;
for (size_type i=0;i<3;++i)
  stream<< t.M_point[i]<<std::endl;
stream<<"****************************************"<<std::endl;
return stream;
}
