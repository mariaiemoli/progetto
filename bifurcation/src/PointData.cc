
#include "../include/PointData.h" 

/**************************************************************************/
/*  PointData.cc														  */
/*  Classe che rappresenta un punto geometrico			                  */
/**************************************************************************/

PointData::PointData(scalar_type x, scalar_type y)
{
	M_coor.clear();
  	M_coor.push_back(x); 
  	M_coor.push_back(y);
  	
}// costruttore

PointData::PointData(const PointData & p)   
{
  	M_coor.clear();
	M_coor.push_back(p.x());
  	M_coor.push_back(p.y());
  	
}// costruttore

void PointData::setCoordinates(scalar_type x, scalar_type y)
{
	M_coor.clear();
  	M_coor.push_back(x); 
  	M_coor.push_back(y);
  	
  	return;
  	
}//setCoordinates

void PointData::getCoordinates(scalar_type & x, scalar_type & y) const
{
  	x=this->x();
  	y=this->y();
  	
  	return;
  	
}//getCoordinates

PointData PointData::operator +=(const PointData & rhs)
{
  M_coor[0]+=rhs.x();
  M_coor[1]+=rhs.y();
  
  return *this;
  
}// operator +=

PointData PointData::operator -=(const PointData & rhs)
{
  M_coor[0]-=rhs.x();
  M_coor[1]-=rhs.y();
  
  return *this;
  
}// operator -=

PointData operator -(const PointData & a, const PointData & b)
{
  return PointData( a.x()-b.x(),
	 				a.y()-b.y() );
  
}// operator -

PointData operator +(const PointData & a, const PointData & b)
{
  return PointData( a.x()+b.x(),
	 				a.y()+b.y() );
  
}// operator +


PointData PointData::operator *(const scalar_type & d) const
{
  return PointData( d*this->x(), d*this->y() );

}// operator *


PointData & PointData::operator =(const PointData & p)
{
  if(this!=&p)
  {
    	M_coor.clear();
  		M_coor.push_back(p.x());
    	M_coor.push_back(p.y());
  }
  
  return *this;
  
}// operator =

PointData operator *(const scalar_type & d, const PointData & p)
{
  return p*d;
  
}// operator *

scalar_type PointData::dot(PointData const & p)const
{
  return this->x()*p.x() + this->y()*p.y();
  
}// operator prodotto scalare

std::ostream & operator <<(std::ostream &stream , PointData const & p)
{
	stream << " ( " << p.x() << ", " << p.y() << " ) " << std::endl;
	
	return stream;
	
}// operator <<