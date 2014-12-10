
#include <iostream>
#include "../include/QuadrilaterData.h"

/**************************************************************************/
/*  QuadrilaterData.cc													  */
/*  Classe che contiene i dati principali per costrire il triangolo       */
/**************************************************************************/

// Definition of the static variable edge
const size_type QuadrilaterData::M_edge[4][2] = {{0, 1}, {1, 2}, {2, 3}, {3,0}};


//COSTRUTTORI	
QuadrilaterData::QuadrilaterData()
{
	M_point.clear();
	M_point.resize(4);
	
}// costruttore vuoto


QuadrilaterData::QuadrilaterData(PointData & a, PointData & b, PointData & c, PointData & d )
{
	M_point.clear();
  	M_point.push_back(a); 
  	M_point.push_back(b);
  	M_point.push_back(c);
  	M_point.push_back(d);
  	
}//costruttore dati i punti

QuadrilaterData::QuadrilaterData(const QuadrilaterData & t)
{
	M_point[0]=t.M_point[ 0 ];
	M_point[1]=t.M_point[ 1 ];
	M_point[2]=t.M_point[ 2 ];
	M_point[3]=t.M_point[ 3 ];

}//Da quadrilatero


void QuadrilaterData::setPoint(size_type i, PointData const & p)
{
	M_point[i] = p;
	
	return;
}//setPoint


scalar_type QuadrilaterData::measure() const
{
	const QuadrilaterData & t = *this;
	
	scalar_type d1,d2;
	
	d1 = sqrt( pow( t[ 2 ][ 0 ] - t[ 0 ][ 0 ], 2) + pow( t[ 2 ][ 1 ] - t[ 0 ][ 1 ], 2) );
	d2 = sqrt( pow( t[ 1 ][ 0 ] - t[ 3 ][ 0 ], 2) + pow( t[ 1 ][ 1 ] - t[ 3 ][ 1 ], 2) );

	return 0.5*d1*d2;
}//measure


PointData & QuadrilaterData::edgePoint(size_type edgenum, size_type endnum)
{
	return M_point[M_edge[edgenum][endnum]];
	
}//edgePoint

const PointData & QuadrilaterData::edgePoint(size_type edgenum, size_type endnum) const
{
	return M_point[M_edge[edgenum][endnum]];
	
}//edgePoint


PointData QuadrilaterData::baricenter()const
{
	PointData tmp(M_point[0]);
	
	tmp += M_point[1];
	tmp += M_point[2];
	tmp += M_point[3];
	
	return tmp * (1./4.0);
	
}//baricenter


size_type QuadrilaterData::edge(size_type i, size_type j)
{
	return M_edge[i][j];
	
}//edge


PointData QuadrilaterData::edgeMean(size_type edgeNum)const
{
	PointData tmp(M_point[edge(edgeNum,0)]+M_point[edge(edgeNum,1)]);
	
	return tmp *0.5;
	
}//edgeMedium


Vector2d QuadrilaterData::c(size_type edgeNum) const
{
	PointData baric = this->baricenter();
	PointData eMean= this->edgeMean(edgeNum);
	
	return eMean.asVector()-baric.asVector();
}//c



Vector2d QuadrilaterData::unscaledNormal(size_type edgeNum) const
{
	PointData const & a=this->edgePoint(edgeNum,0);
	PointData const & b=this->edgePoint(edgeNum,1);

	return Vector2d(a.y()-b.y(),b.x()-a.x());
}//unscaledNormal


QuadrilaterData & QuadrilaterData::operator =(const QuadrilaterData & t)
{	
  	if (this !=&t)
   	{
  		M_point[0] = t.M_point[0];
  		M_point[1] = t.M_point[1];
  		M_point[2] = t.M_point[2];
  		M_point[3] = t.M_point[3];
    }
	
	return *this;
}// operator =

std::ostream & operator <<(std::ostream & stream, QuadrilaterData const & t)
{
	stream<<" ************* QUADRILATER  POINTS ********"<<std::endl;
	for (size_type i=0;i<4;++i)
	{
		stream<< t.M_point[i]<<std::endl;
	}
	stream<<"****************************************"<<std::endl;
	
	return stream;
}// operator <<
