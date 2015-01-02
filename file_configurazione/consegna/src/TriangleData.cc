
#include <iostream>
#include "../include/TriangleData.h"

/**************************************************************************/
/*  TriangleData.cc														  */
/*  Classe che contiene i dati principali per costrire il triangolo       */
/**************************************************************************/

// Definition of the static variable edge
const size_type TriangleData::M_edge[3][2] = {{0, 1}, {1, 2}, {2, 0}};


//COSTRUTTORI	
TriangleData::TriangleData()
{
	M_point.clear();
	M_point.resize(3);
	
}// costruttore vuoto


TriangleData::TriangleData(PointData & a, PointData & b, PointData & c)
{
	M_point.clear();
  	M_point.push_back(a); 
  	M_point.push_back(b);
  	M_point.push_back(c);
  	
}//costruttore dati i punti

TriangleData::TriangleData(const TriangleData & t)
{
	M_point[0]=t.M_point[ 0 ];
	M_point[1]=t.M_point[ 1 ];
	M_point[2]=t.M_point[ 2 ];

}//Da triangolo

void TriangleData::setPoint(size_type i, PointData const & p)
{
	M_point[i] = p;
	
	return;
}//setPoint

scalar_type TriangleData::measure() const
{
	const TriangleData & t = *this;
	
	return -0.5 * ( t[1][0] * (t[2][1] - t[0][1]) +
	   				t[2][0] * (t[0][1] - t[1][1]) +
	   				t[0][0] * (t[1][1] - t[2][1])   );
}//measure

size_type TriangleData::size () const
{
	return M_point.size();
}

PointData & TriangleData::edgePoint(size_type edgenum, size_type endnum)
{
	return M_point[M_edge[edgenum][endnum]];
	
}//edgePoint

const PointData & TriangleData::edgePoint(size_type edgenum, size_type endnum) const
{
	return M_point[M_edge[edgenum][endnum]];
	
}//edgePoint

PointData TriangleData::baricenter()const
{
	PointData tmp(M_point[0]);
	
	tmp += M_point[1];
	tmp += M_point[2];
	
	return tmp * (1./3.0);
	
}//baricenter

size_type TriangleData::edge(size_type i, size_type j)
{
	return M_edge[i][j];
	
}//edge

PointData TriangleData::edgeMean(size_type edgeNum)const
{
	PointData tmp(M_point[edge(edgeNum,0)]+M_point[edge(edgeNum,1)]);
	
	return tmp *0.5;
	
}//edgeMedium

Vector2d TriangleData::c(size_type edgeNum) const
{
	PointData baric = this->baricenter();
	PointData eMean= this->edgeMean(edgeNum);
	
	return eMean.asVector()-baric.asVector();
}//c

Vector2d TriangleData::unscaledNormal(size_type edgeNum) const
{
	PointData const & a=this->edgePoint(edgeNum,0);
	PointData const & b=this->edgePoint(edgeNum,1);

	return Vector2d(a.y()-b.y(),b.x()-a.x());
}//unscaledNormal


TriangleData & TriangleData::operator =(const TriangleData & t)
{
  	if (this !=&t)
   	{
  		M_point[0] = t.M_point[0];
  		M_point[1] = t.M_point[1];
  		M_point[2] = t.M_point[2];
    }
	
	return *this;
}// operator =

std::ostream & operator <<(std::ostream & stream, TriangleData const & t)
{
	stream<<" ************* TRIANGLE  POINTS ********"<<std::endl;
	for (size_type i=0;i<3;++i)
	{
		stream<< t.M_point[i]<<std::endl;
	}
	stream<<"****************************************"<<std::endl;
	
	return stream;
}// operator <<
