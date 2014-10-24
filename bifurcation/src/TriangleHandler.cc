#include <iostream>
#include "../include/TriangleHandler.h"


// Definition of the static variable edge

const int Triangle::M_edge[3][2] = {{0, 1}, {1, 2}, {2, 0}};

Triangle::Triangle( const FracturePtrContainer_Type& M_fractures )
{
	assert ( M_fractures.size() == 3 );
	
	FractureHandlerPtr_Type f0 = M_fractures [ 0 ];
	FractureHandlerPtr_Type f1 = M_fractures [ 1 ];
	FractureHandlerPtr_Type f2 = M_fractures [ 2 ];

	LevelSetHandlerPtr_Type l0 = M_fractures [ 0 ]-> getLevelSet ();
	LevelSetHandlerPtr_Type l1 = M_fractures [ 1 ]-> getLevelSet ();
	LevelSetHandlerPtr_Type l2 = M_fractures [ 2 ]-> getLevelSet ();

	
	std::string val0;
	std::string val1;
	
	/*
	stringContainer_Type M_subRegion;
		
	M_subRegion [ 0 ] = "++";
	M_subRegion [ 1 ] = "--";
	M_subRegion [ 2 ] = "+-";
	M_subRegion [ 3 ] = "-+";

    sizeVectorContainer_Type listOfLevelSet0 (2);
    sizeVectorContainer_Type listOfLevelSet1 (2);
    
    listOfLevelSet0 [ 0 ] = M_fractures [ 0 ]->getID();
    listOfLevelSet1 [ 0 ] = M_fractures [ 0 ]->getID();
    
    listOfLevelSet0 [ 1 ] = M_fractures [ 1 ]->getID();
    listOfLevelSet1 [ 1 ] = M_fractures [ 2 ]->getID();
	
    for ( size_type i = 0; i < maxIntersection; ++i )
    {
        booleanOperation [i] = getOperation ( M_subRegion[i], levelSets );

    }

  */
	
	const size_type nbDof0 = l0 ->get_mesh_fem().nb_basic_dof();
	const size_type nbDof1 = l1 ->get_mesh_fem().nb_basic_dof();
	const size_type nbDof2 = l2 ->get_mesh_fem().nb_basic_dof();
	
	base_node p0 = l0-> get_mesh_fem().point_of_basic_dof( 0 );
	base_node p1 = l1-> get_mesh_fem().point_of_basic_dof( 0 );
	base_node p2 = l2-> get_mesh_fem().point_of_basic_dof( 0 );

	if ( l1->getData()->levelSetFunction ( p0 ) == 0 )
	{
		p0 = l0-> get_mesh_fem().point_of_basic_dof( nbDof0-1 )
	}

	if ( l0->getData()->levelSetFunction ( p1 ) == 0 )
	{
		p1 = l1-> get_mesh_fem().point_of_basic_dof( nbDof1-1 )
	}

	if ( l1->getData()->levelSetFunction ( p2 ) == 0 )
	{
		p2 = l2-> get_mesh_fem().point_of_basic_dof( nbDof2-1 )
	}

	
    val0 = l0->getData()->levelSetFunction ( p1 );
    val1 = l0->getData()->levelSetFunction ( p2 );
    
    if( val0 * val1 < 0 )
    {
    	//questa è la prima frattura da usare insieme a quella che sta nella sua parte negativa ( per esempio )
    	
    }
    else
    {
        val0 = l1->getData()->levelSetFunction ( p0 );
        val1 = l1->getData()->levelSetFunction ( p2 );
        
        if ( val0 * val1 < 0 )
        {
        	//questa è la prima frattura da usare insieme a quella che sta nella sua parte negativa ( per esempio )
        }
        else
        {
        	// è sicuramente la terza quella da usare
        }
    }
    

	
}


Triangle::Triangle(Point & a, Point & b, Point & c)
{
	M_points[0] = a;
	M_points[1] = b;
	M_points[2] = c;
}

void Triangle::setPoint(int i, Point const & p)
{
M_points[i] = p;
}

Point & Triangle::edgePoint(int edgenum, int endnum)
{
return M_points[M_edge[edgenum][endnum]];
}

const Point & Triangle::edgePoint(int edgenum, int endnum) const
{
return M_points[M_edge[edgenum][endnum]];
}

Triangle::Triangle(const Triangle & t)
{
M_points[0]=t.M_points[0];
M_points[1]=t.M_points[1];
	M_points[2]=t.M_points[2];
}

Triangle & Triangle::operator =(const Triangle & t)
{
if (this !=&t){
  M_points[0]=t.M_points[0];
  M_points[1]=t.M_points[1];
  M_points[2]=t.M_points[2];
}
return *this;
}

int Triangle::edge(int i, int j)
{
return M_edge[i][j];
}

double Triangle::measure() const
{
const Triangle & t = *this;
return -0.5 *
  (t[1][0] * (t[2][1] - t[0][1]) +
   t[2][0] * (t[0][1] - t[1][1]) +
   t[0][0] * (t[1][1] - t[2][1]));
}

Point Triangle::baricenter()const
{
Point tmp(M_points[0]);
tmp += M_points[1];
tmp += M_points[2];
return tmp * (1./3.0);
}

Point Triangle::edgeBaricenter(int edgeNum)const
{
Point tmp(M_points[edge(edgeNum,0)]+M_points[edge(edgeNum,1)]);
return tmp *0.5;
}

Vector Triangle::c(int edgeNum) const
{
Point baric = this->baricenter();
Point eBaric= this->edgeBaricenter(edgeNum);
return eBaric.asVector()-baric.asVector();
}

Vector Triangle::unscaledNormal(int edgeNum) const
{
Point const & a=this->edgePoint(edgeNum,0);
Point const & b=this->edgePoint(edgeNum,1);
return Vector(a.y()-b.y(),b.x()-a.x());
}

std::ostream & operator <<(std::ostream & stream, Triangle const & t)
{
stream<<" ************* TRIANGLE  POINTS ********"<<std::endl;
for (int i=0;i<3;++i)
  stream<< t.M_points[i]<<std::endl;
stream<<"****************************************"<<std::endl;
return stream;
}
Point Point::operator +=(const Point & rhs)
{
M_coor[0]+=rhs.M_coor[0];
M_coor[1]+=rhs.M_coor[1];
return *this;
}

Point Point::operator -=(const Point & rhs)
{
M_coor[0]-=rhs.M_coor[0];
M_coor[1]-=rhs.M_coor[1];
return *this;
}

Point operator -(const Point & a, const Point & b)
{
return Point(a.M_coor[0]-b.M_coor[0],
	 a.M_coor[1]-b.M_coor[1]);
}

Point operator +(const Point & a, const Point & b)
{
return Point(a.M_coor[0]+b.M_coor[0],
	 a.M_coor[1]+b.M_coor[1]);
}


Point Point::operator *(const double & d) const
{
return Point(d*M_coor[0],d*M_coor[1]);
}

Point::Point(const Point & p)    {
M_coor[0]=p.M_coor[0];
M_coor[1]=p.M_coor[1];
}

Point & Point::operator =(const Point & p)
{
if(this!=&p){
  M_coor[0]=p.M_coor[0];
  M_coor[1]=p.M_coor[1];
}
return *this;
}

Point operator *(const double & d, const Point & p)
{
return p*d;
}

double Point::dot(Point const & p)const
{
return this->x()*p.x()+this->y()*p.y();
}

std::ostream & operator <<(std::ostream &stream , Point const & p)
{
  stream <<"("<<p.M_coor[0]<<","<<p.M_coor[1]<<")"<<std::endl;
  return stream;
}

