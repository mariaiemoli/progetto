/**
 * MatrixBifurcationHandler.h
 * 
 * Classe che gestisce la costruzione delle matrici per il triangolo di intersezione
 * 
 */

#ifndef __MATRIXBIFURCATIONHANDLER_H__
#define __MATRIXBIFURCATIONHANDLER_H__

#include "Core.h"
#include "TriangleHandler.h"
#include <Eigen/Dense>



class Bifurcation{
public:
  //! Default constructor
  Bifurcation();
  //! Passing the permeability tensor
  Bifurcation(Matrix2d K);
  //! Passing a scalar permeability
  explicit Bifurcation(scalar_type K);
  //! Give the triangle.
  void setTriangle(TriangleHandler const & T){triangle_=T;}
  //! Change the permeability.
  void setK(Matrix2d K){K_=K;}
  //! Get permeability.
  Matrix2d K()const{ return K_;}
  //! Compute matrix N having unscaled normals as rows.
  /*! 
    The internal triangle must have been set before calling this method.
  */
  void computeN();
  //! Get matrix N;
  Matrix32 N()const {return N_;}
  //! compute matrix C.
  /*!
    The matrix containing as rows the vector from triangle baricenter
    to triangle edge baricenter.
  */
  void computeC();
  //! get matrix C;
  Matrix32 C(){return C_;}
  //! compute matric Qc
  /*!
    The columns form an orthonormal basis for the transpose of C.
    It computes also C
  */
  void computeQc();
  //! get Qc.
  Matrix32 Qc()const {return Qc_;}
  //! Compute Projection matrix on the null space of Qc.
  /*
    Empty method since for a triangle PC is given
  */
  void computePc();
  //! Get Pc. I computes also Qc and C.
  Matrix3d Pc()const {return Pc_;}
  //! Computes the Trasmissibility matrix
  /*!
    It computes also C, N e Qc
    It uses the formula
    \f$ T= N K N^T + t P_c \operatorname{diag}(N K N^T) P_c
  */
  void computeT(scalar_type t=6.0);
  //! Computes the Trasmissibility matrix
  /*!
    It computes also C, N e Qc
    It uses the formula
    \f$ T= N K N^T + t \operatorname{trace}(NKN^T)P_c
  */
  void computeTsimple(scalar_type t=6.0);
  Matrix3d T()const {return T_;}
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW 

private:
	  Matrix2d K_;
	  TriangleHandler triangle_;
	  Matrix32 N_;
	  Matrix32 C_;
	  Matrix32 Qc_;
	  Matrix3d Pc_;
	  Matrix3d T_;
};


//! Returns the coefficients that links the pressure to edge pressure
/*!
  It computes the coefficients that links the pressure internal at the 
  intersection with that at the boundary (fracture end pressures).

  @param T transmissibility matrix.
  @return The vector the gives the coefficients.
 */
Vector3d pressureCoeff(const Matrix3d & T);


#endif