
#include "../include/MatrixBifurcationHandler.h" 

/**************************************************************************/
/*  MatrixBifurcationHandler.cc											  */
/*  Classe che gestisce la costruzione delle matrici per il triangolo di  */
/*  intersezione                 										  */
/**************************************************************************/

MatrixBifurcationHandler::MatrixBifurcationHandler( const GetPot& dataFile,
													const std::string& type, 
 												    const std::string& section,
												    const std::string& subsection ):
												    M_Pc(Matrix4d::Constant(1./4.0)),
												    M_type( type )
{
	Matrix2d invk; 
	scalar_type det = dataFile ( (section + subsection + "invK").data(), 1. );
	
	invk(0,0) = dataFile ( ( section + subsection + "invKDist11" ).data (), 1. );
	invk(0,1) = dataFile ( ( section + subsection + "invKDist12" ).data (), 0. );
	invk(1,0) = dataFile ( ( section + subsection + "invKDist12" ).data (), 0. );
	invk(1,1) = dataFile ( ( section + subsection + "invKDist22" ).data (), 1. );
	
	this-> inversion2X2 ( invk, det );
	
}// costruttore


MatrixBifurcationHandler::MatrixBifurcationHandler(scalar_type K, const std::string& type):
												   M_Pc(Matrix4d::Constant(1./4.0)),
												   M_type ( type) 
{
	M_K(0,0)=K;
	M_K(0,1)=M_K(1,0)=0.0;
	M_K(1,1)=K;
}// costruttore

MatrixBifurcationHandler::MatrixBifurcationHandler(Matrix2d K, const std::string& type ):
												M_K(K),M_Pc(Matrix4d::Constant(1./4.0)), 
												M_type( type)
{};


void MatrixBifurcationHandler::setMatrices ( FracturePtrContainer_Type& fractures )
{
	M_intersection.setIntersection ( fractures );
	
	
	if ( M_type == "Bifurcation" )
	{
		for ( size_type i = 0; i < 4; i++ )
		{
			for ( size_type j = 0; j < 4; j++ )
			{
				if( i == 3 || j == 3)
				{
					M_Pc(i,j) = 0;
				}
				else
				{
					M_Pc(i,j) = 1./3.0;
				}
			}
		}
	}
	
	computeT();	
	
	return;

}// setMatrices


void MatrixBifurcationHandler::computeN()
{
	if ( M_type == "Bifurcation2")
	{
		M_N.row(0)=M_intersection.intersectionQuadrilater ().unscaledNormal(0).transpose();
		M_N.row(1)=M_intersection.intersectionQuadrilater ().unscaledNormal(1).transpose();
		M_N.row(2)=M_intersection.intersectionQuadrilater ().unscaledNormal(2).transpose();
		M_N.row(3)=M_intersection.intersectionQuadrilater ().unscaledNormal(3).transpose();
		
	}
	else
	{
		M_N.row(0)=M_intersection.intersectionTriangle ().unscaledNormal(0).transpose();
		M_N.row(1)=M_intersection.intersectionTriangle ().unscaledNormal(1).transpose();
		M_N.row(2)=M_intersection.intersectionTriangle ().unscaledNormal(2).transpose();
		M_N.row(3)=0*M_intersection.intersectionTriangle ().unscaledNormal(2).transpose();
	}

	return;
}// computeN

void MatrixBifurcationHandler::computeC()
{
	if ( M_type == "Bifurcation" )
	{	
		M_C.row(0) = M_intersection.intersectionTriangle ().c(0).transpose();
		M_C.row(1) = M_intersection.intersectionTriangle ().c(1).transpose();
		M_C.row(2) = M_intersection.intersectionTriangle ().c(2).transpose();
		M_C.row(3) = 0*M_intersection.intersectionTriangle ().c(2).transpose();
	}
	else
	{
		M_C.row(0) = M_intersection.intersectionQuadrilater ().c(0).transpose();
		M_C.row(1) = M_intersection.intersectionQuadrilater ().c(1).transpose();
		M_C.row(2) = M_intersection.intersectionQuadrilater ().c(2).transpose();
		M_C.row(3) = M_intersection.intersectionQuadrilater ().c(3).transpose();
	}

	return;
}// computeC


void MatrixBifurcationHandler::computeQc()
{
	this->computeC();
	M_Qc.col(0)  = M_C.col(0).normalized();
	Vector4d v0 = M_Qc.col(0);
	Vector4d v1 = M_C.col(1) - ( v0.dot(M_C.col(1)) )*v0;
	M_Qc.col(1)  = v1.normalized();

	return;
}// computeQc


void MatrixBifurcationHandler::computeT(scalar_type t)
{
	this->computeQc();
	this->computeN();

	scalar_type area;
	
	if ( M_type == "Bifurcation")
	{
		area = M_intersection.intersectionTriangle ().measure();
	}
	else 
	{
		area = M_intersection.intersectionQuadrilater ().measure();
	}
	
	Matrix4d Nkn = M_N*M_K*M_N.transpose();
	
	Eigen::DiagonalMatrix<scalar_type,4,4> Nd (Nkn.diagonal());
	
	Matrix4d tmp = M_Pc *Nd *M_Pc;
	
	M_T=(1./area)*( Nkn + t*tmp );
		
	return;
}// computeT


void MatrixBifurcationHandler::computeTsimple(scalar_type t)
{
	this->computeQc();
	this->computeN();

	scalar_type area = M_intersection.intersectionTriangle ().measure();
	Matrix4d Nkn = M_N*M_K*M_N.transpose();
	
	scalar_type Nd=Nkn.trace();
	
	Matrix4d tmp = Nd *M_Pc;
	
	M_T=(1./area)*( Nkn + t*tmp );
	
	return;
}// computeTsimple


void MatrixBifurcationHandler::computeScap( scalar_type& s, scalar_type t )
{
	Matrix4d Nkn = M_N*M_K*M_N.transpose();
	
	Eigen::DiagonalMatrix<scalar_type,4,4> Nd (Nkn.diagonal());
	
	Matrix4d S = t*Nd;
	
	if ( M_intersection.intersectionTriangle ().measure() == 0 )
	{
		s = S.trace()*1./4.0; 
	}
	else
	{
		s = S.trace()*1./3.0;
	}
	
	return; 
}// computeS


void MatrixBifurcationHandler::inversion2X2 ( const Matrix2d& invk, scalar_type det )
{
	this-> M_K ( 0, 0 ) = invk ( 1, 1 )*1./det;
	this-> M_K ( 0, 1 ) = -invk ( 0, 1 )*1./det;
	this-> M_K ( 1, 0 ) = -invk ( 1, 0 )*1./det;
	this-> M_K ( 1, 1 ) = invk ( 0, 0 )*1./det;
	
	return;
}// inversion2X2


void MatrixBifurcationHandler::SetDOFIntersecton( FractureHandlerPtr_Type& fractures, scalar_type& DOF )
{
	const size_type nbDof =  fractures-> getMeshFEMPressure().nb_basic_dof();
	
	base_node node(2);
	
	base_node tmp0( 1 );
	base_node tmp1( 1 );
	tmp0[ 0 ]= 0.;
	tmp1[ 0 ]= 1.;
	
	node[ 0 ] = fractures-> getMeshFEMPressure().point_of_basic_dof( 0 )[ 0 ];
	node[ 1 ] = fractures-> getLevelSet()->getData()->y_map( tmp0 );
		
	PointData temp = M_intersection.getPointIntersection();
	if( fractures->getDOFBifurcation().size() != 0 )
	{
		DOF = fractures->getDOFBifurcation()[ 0 ];
	}
	else
	{
		if( ( gmm::abs( node[ 0 ] - temp.x() ) < 1.0E-2 ) && ( gmm::abs( node[ 1 ] - temp.y() ) < 1.0E-2 ) )
		{
			DOF = 0;
		}
		else 
		{
			DOF= nbDof - 1;
		}
	}
	
	return;
			
}//SetDOFIntersecton