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

Vector3d pressureCoeff(const Matrix3d &T)
{
  scalar_type factor = (1.0/T.sum());
  Vector3d tmp;
  tmp(0)=factor*(T.col(0).sum());
  tmp(1)=factor*(T.col(1).sum());
  tmp(2)=factor*(T.col(2).sum());
  return tmp;
}