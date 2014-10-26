#include "../include/PointData.h" 

PointData::PointData(scalar_type x, scalar_type y)
{
  M_coor[0]=x; 
  M_coor[1]=y;
}

void PointData::setCoordinates(scalar_type x, scalar_type y)
{
  M_coor[0]=x; 
  M_coor[1]=y;
}//setCoordinates

void PointData::getCoordinates(scalar_type & x, scalar_type & y) const
{
  x=M_coor[0];
  y=M_coor[1];
}//getCoordinates

PointData PointData::operator +=(const PointData & rhs)
{
  M_coor[0]+=rhs.M_coor[0];
  M_coor[1]+=rhs.M_coor[1];
  return *this;
}

PointData PointData::operator -=(const PointData & rhs)
{
  M_coor[0]-=rhs.M_coor[0];
  M_coor[1]-=rhs.M_coor[1];
  return *this;
}

PointData operator -(const PointData & a, const PointData & b)
{
  return PointData( a.M_coor[0]-b.M_coor[0],
	 				   a.M_coor[1]-b.M_coor[1] );
}

PointData operator +(const PointData & a, const PointData & b)
{
  return PointData(a.M_coor[0]+b.M_coor[0],
	 a.M_coor[1]+b.M_coor[1]);
}


PointData PointData::operator *(const scalar_type & d) const
{
  return PointData( d*M_coor[0], d*M_coor[1] );
}

PointData::PointData(const PointData & p)   
{
  M_coor[0]=p.M_coor[0];
  M_coor[1]=p.M_coor[1];
}

PointData & PointData::operator =(const PointData & p)
{
  if(this!=&p)
  {
    M_coor[0]=p.M_coor[0];
    M_coor[1]=p.M_coor[1];
  }
  return *this;
}

PointData operator *(const scalar_type & d, const PointData & p)
{
  return p*d;
}

scalar_type PointData::dot(PointData const & p)const
{
  return this->x()*p.x()+this->y()*p.y();
}

std::ostream & operator <<(std::ostream &stream , PointData const & p)
{
	stream <<"("<<p.M_coor[0]<<","<<p.M_coor[1]<<")"<<std::endl;
	return stream;
}