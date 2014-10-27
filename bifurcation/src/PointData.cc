#include "../include/PointData.h" 

PointData::PointData(scalar_type x, scalar_type y)
{
	M_coor.clear();
  	M_coor.push_back(x); 
  	M_coor.push_back(y);
}

PointData::PointData(const PointData & p)   
{
  	M_coor.clear();
	M_coor.push_back(p.x());
  	M_coor.push_back(p.y());
}

void PointData::setCoordinates(scalar_type x, scalar_type y)
{
	M_coor.clear();
  	M_coor.push_back(x); 
  	M_coor.push_back(y);
}//setCoordinates

void PointData::getCoordinates(scalar_type & x, scalar_type & y) const
{
  	x=this->x();
  	y=this->y();
}//getCoordinates

PointData PointData::operator +=(const PointData & rhs)
{
  M_coor[0]+=rhs.x();
  M_coor[1]+=rhs.y();
  return *this;
}

PointData PointData::operator -=(const PointData & rhs)
{
  M_coor[0]-=rhs.x();
  M_coor[1]-=rhs.y();
  return *this;
}

PointData operator -(const PointData & a, const PointData & b)
{
  return PointData( a.x()-b.x(),
	 				a.y()-b.y() );
}

PointData operator +(const PointData & a, const PointData & b)
{
  return PointData( a.x()+b.x(),
	 				a.y()+b.y() );
}


PointData PointData::operator *(const scalar_type & d) const
{
  return PointData( d*this->x(), d*this->y() );
}


PointData & PointData::operator =(const PointData & p)
{
  if(this!=&p)
  {
    	M_coor.clear();
  		M_coor.push_back(p.x());
    	M_coor.push_back(p.y());
  }
  return *this;
}

PointData operator *(const scalar_type & d, const PointData & p)
{
  return p*d;
}

scalar_type PointData::dot(PointData const & p)const
{
  return this->x()*p.x() + this->y()*p.y();
}

std::ostream & operator <<(std::ostream &stream , PointData const & p)
{
	stream << " ( " << p.x() << ", " << p.y() << " ) " << std::endl;
	return stream;
}