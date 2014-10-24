#include "../include/PointHandler.h" 

PointHandler::PointHandler(scalar_type x, scalar_type y){
  M_coor[0]=x; M_coor[1]=y;
}
void PointHandler::setCoordinates(scalar_type x, scalar_type y){
  M_coor[0]=x; M_coor[1]=y;
}
void PointHandler::getCoordinates(scalar_type & x, scalar_type & y) const{
  x=M_coor[0];y=M_coor[1];
}

PointHandler PointHandler::operator +=(const PointHandler & rhs)
{
  M_coor[0]+=rhs.M_coor[0];
  M_coor[1]+=rhs.M_coor[1];
  return *this;
}

PointHandler PointHandler::operator -=(const PointHandler & rhs)
{
  M_coor[0]-=rhs.M_coor[0];
  M_coor[1]-=rhs.M_coor[1];
  return *this;
}

PointHandler operator -(const PointHandler & a, const PointHandler & b)
{
  return PointHandler(a.M_coor[0]-b.M_coor[0],
	 a.M_coor[1]-b.M_coor[1]);
}

PointHandler operator +(const PointHandler & a, const PointHandler & b)
{
  return PointHandler(a.M_coor[0]+b.M_coor[0],
	 a.M_coor[1]+b.M_coor[1]);
}


PointHandler PointHandler::operator *(const scalar_type & d) const
{
  return PointHandler(d*M_coor[0],d*M_coor[1]);
}

PointHandler::PointHandler(const PointHandler & p)    {
  M_coor[0]=p.M_coor[0];
  M_coor[1]=p.M_coor[1];
}

PointHandler & PointHandler::operator =(const PointHandler & p)
{
  if(this!=&p){
    M_coor[0]=p.M_coor[0];
    M_coor[1]=p.M_coor[1];
  }
  return *this;
}

PointHandler operator *(const scalar_type & d, const PointHandler & p)
{
  return p*d;
}

scalar_type PointHandler::dot(PointHandler const & p)const
{
  return this->x()*p.x()+this->y()*p.y();
}

std::ostream & operator <<(std::ostream &stream , PointHandler const & p)
{
	stream <<"("<<p.M_coor[0]<<","<<p.M_coor[1]<<")"<<std::endl;
	return stream;
}