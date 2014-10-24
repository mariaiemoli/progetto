#include "../include/MatrixBifurcationHandler.h" 


Bifurcation::Bifurcation():Pc_(Matrix3d::Constant(1./3.0))

{
}

Bifurcation::Bifurcation(scalar_type K):Pc_(Matrix3d::Constant(1./3.0))
{
  K_(0,0)=K;
  K_(0,1)=K_(1,0)=0.0;
  K_(1,1)=K;
}

Bifurcation::Bifurcation(Matrix2d K):K_(K),Pc_(Matrix3d::Constant(1./3.0)){};

void Bifurcation::computeN()
{
  N_.row(0)=triangle_.unscaledNormal(0).transpose();
  N_.row(1)=triangle_.unscaledNormal(1).transpose();
  N_.row(2)=triangle_.unscaledNormal(2).transpose();
}

void Bifurcation::computeC()
{
  C_.row(0) = triangle_.c(0).transpose();
  C_.row(1) = triangle_.c(1).transpose();
  C_.row(2) = triangle_.c(2).transpose();
}

void Bifurcation::computeQc()
{
  this->computeC();
  Qc_.col(0)  = C_.col(0).normalized();
  Vector3d v0 = Qc_.col(0);
  Vector3d v1 = C_.col(1) - ( v0.dot(C_.col(1)) )*v0;
  Qc_.col(1)  = v1.normalized();
}

void Bifurcation::computePc()
{
  /*
  this->computeQc();
  Pc_= Matrix3d::Identity() - Qc_*Qc_.transpose();
  */
  // For a triangle we know how to do it!
}

void Bifurcation::computeT(scalar_type t)
{
  //this->computePc();
  this->computeQc();
  this->computeN();
  scalar_type area = triangle_.measure();
  Matrix3d Nkn = N_*K_*N_.transpose();
  Eigen::DiagonalMatrix<scalar_type,3,3> Nd (Nkn.diagonal());
  Matrix3d tmp = Pc_ *Nd *Pc_;
  T_=(1./area)*( Nkn + t*tmp );
}

void Bifurcation::computeTsimple(scalar_type t)
{
  //this->computePc();
  this->computeQc();
  this->computeN();
  scalar_type area = triangle_.measure();
  Matrix3d Nkn = N_*K_*N_.transpose();
  scalar_type Nd=Nkn.trace();
  Matrix3d tmp = Nd *Pc_;
  T_=(1./area)*( Nkn + t*tmp );
}
 
Intersection::Intersection(FractureEnd const & gamma0, FractureEnd const & gamma1, FractureEnd const & gamma2, PointHandler const & intersectionPoint):
  fractures{gamma0,gamma1,gamma2}
{
  Vector2d tmp;
  for (unsigned int i=0;i<3;++i)
    {
		// get tangent at the end
		tmp(0) = intersectionPoint.x()-fractures[i].endPoint.x();
		tmp(1) = intersectionPoint.y()-fractures[i].endPoint.y();
		tmp.normalize();
		tangents[i]   = tmp;
		normals[i](0) = -tmp(1);
		normals[i](1) =  tmp(0);
    }
}

TriangleHandler const &
Intersection::computeIntersectionTriangle()
{
  PointHandler point;
  const scalar_type tol=1.e-5;
  scalar_type s(0.);
  for (unsigned int i=0; i<3;++i)
    {
		unsigned int j = (i +1 ) % 3;
		scalar_type ninj = normals[i].dot(normals[j] );
		scalar_type nitj = normals[i].dot(tangents[j]);
		if(std::fabs(nitj)<tol)
  {
    point=intersection_+PointHandler(normals[i](0),normals[i](1))*(0.5*(fractures[j].thickness+fractures[i].thickness));
  }
else
  {
    // parametric coordinate
    s    = 0.5*(
		   fractures[j].thickness*ninj +
		   fractures[i].thickness
		   )/nitj;

    // the ith point is Pji in the note
    point = intersection_+ 
      (PointHandler(tangents[j](0),tangents[j](1))*s)- 
      (PointHandler(normals[j](0),normals[j](1))*(0.5*fractures[j].thickness));
  }
     intersectionTriangle_.setPoint(i,point);
   }
return intersectionTriangle_;
}
Vector3d pressureCoeff(const Matrix3d &T)
{
  scalar_type factor = (1.0/T.sum());
  Vector3d tmp;
  tmp(0)=factor*(T.col(0).sum());
  tmp(1)=factor*(T.col(1).sum());
  tmp(2)=factor*(T.col(2).sum());
  return tmp;
}