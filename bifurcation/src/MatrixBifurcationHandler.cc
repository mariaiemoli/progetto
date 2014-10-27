#include "../include/MatrixBifurcationHandler.h" 


MatrixBifurcationHandler::MatrixBifurcationHandler( const GetPot& dataFile,
 												    const std::string& section,
												    const std::string& subsection):Pc_(Matrix3d::Constant(1./3.0))
{
	Matrix2d invk; 
	scalar_type det = dataFile ( (section + subsection + "invK").data(), 1. );
	
	invk(0,0) = dataFile ( ( section + subsection + "invKDist11" ).data (), 1. );
	invk(0,1) = dataFile ( ( section + subsection + "invKDist12" ).data (), 0. );
	invk(1,0) = dataFile ( ( section + subsection + "invKDist12" ).data (), 0. );
	invk(1,1) = dataFile ( ( section + subsection + "invKDist22" ).data (), 1. );
	
	
	this-> inversion2X2 ( invk, det );
	
}// costruttore


MatrixBifurcationHandler::MatrixBifurcationHandler(scalar_type K):Pc_(Matrix3d::Constant(1./3.0))
{
	K_(0,0)=K;
	K_(0,1)=K_(1,0)=0.0;
	K_(1,1)=K;
}// costruttore

MatrixBifurcationHandler::MatrixBifurcationHandler(Matrix2d K):K_(K),Pc_(Matrix3d::Constant(1./3.0)){};


void MatrixBifurcationHandler::setMatrices ( FracturePtrContainer_Type& fractures )
{
	M_intersection.setIntersection ( fractures );
	
	computeT();	
}// setMatrices


void MatrixBifurcationHandler::computeN()
{
	N_.row(0)=M_intersection.intersectionTriangle ().unscaledNormal(0).transpose();
	N_.row(1)=M_intersection.intersectionTriangle ().unscaledNormal(1).transpose();
	N_.row(2)=M_intersection.intersectionTriangle ().unscaledNormal(2).transpose();
}// computeN

void MatrixBifurcationHandler::computeC()
{
	C_.row(0) = M_intersection.intersectionTriangle ().c(0).transpose();
	C_.row(1) = M_intersection.intersectionTriangle ().c(1).transpose();
	C_.row(2) = M_intersection.intersectionTriangle ().c(2).transpose();
}// computeC


void MatrixBifurcationHandler::computeQc()
{
	this->computeC();
	Qc_.col(0)  = C_.col(0).normalized();
	Vector3d v0 = Qc_.col(0);
	Vector3d v1 = C_.col(1) - ( v0.dot(C_.col(1)) )*v0;
	Qc_.col(1)  = v1.normalized();
}// computeQc


void MatrixBifurcationHandler::computeT(scalar_type t)
{
	this->computeQc();
	this->computeN();

	scalar_type area = M_intersection.intersectionTriangle ().measure();
	Matrix3d Nkn = N_*K_*N_.transpose();
	
	Eigen::DiagonalMatrix<scalar_type,3,3> Nd (Nkn.diagonal());
	
	Matrix3d tmp = Pc_ *Nd *Pc_;
	
	T_=(1./area)*( Nkn + t*tmp );
}// computeT


void MatrixBifurcationHandler::computeTsimple(scalar_type t)
{
	this->computeQc();
	this->computeN();

	scalar_type area = M_intersection.intersectionTriangle ().measure();
	Matrix3d Nkn = N_*K_*N_.transpose();
	
	scalar_type Nd=Nkn.trace();
	
	Matrix3d tmp = Nd *Pc_;
	
	T_=(1./area)*( Nkn + t*tmp );
}// computeTsimple


void MatrixBifurcationHandler::inversion2X2 ( const Matrix2d& invk, scalar_type det )
{
	this-> K_ ( 0, 0 ) = invk ( 1, 1 )*1./det;
	this-> K_ ( 0, 1 ) = -invk ( 0, 1 )*1./det;
	this-> K_ ( 1, 0 ) = -invk ( 1, 0 )*1./det;
	this-> K_ ( 1, 1 ) = invk ( 0, 0 )*1./det;
}// inversion2X2

