/*
 * PROGETTO DI PACS 2014
 *
 * \author Bonomi Claudia
 * 
 * \author Iemoli Maria
 *
 * Problema di Darcy per un network di fratture
 *
 */

#ifndef _DARCY_OPERATORSX_
#define _DARCY_OPERATORSX_ 1

#include "Core.h"
#include "FracturesSet.h"
#include "LevelSetHandler.h"
#include "MeshHandler.h"
#include "FractureHandler.h"
#include "BCHandler.h"
#include "MatrixBifurcationHandler.h"

/**************************************************************************/
/*  XFEMOperators.h													  	  */
/*  Libreria che definisce le forme lineari e bilineari per il problema	  */
/*  di Darcy             												  */
/**************************************************************************/


namespace getfem
{


// Classe che rappresenta la normale ad un levelset
class level_set_unit_normal : public getfem::nonlinear_elem_term
{
public:

    level_set_unit_normal ( const getfem::mesh_fem& mf_,
                            const scalarVector_Type& U_ );

    inline const bgeot::multi_index& sizes ( ) const
    {
        return sizes_;
    }

    virtual void compute ( getfem::fem_interpolation_context& ctx,
                           bgeot::base_tensor& t );
    
private:
    const getfem::mesh_fem& mf;
    scalarVector_Type U;
    size_type N;
    base_matrix gradU;
    bgeot::base_vector coeff;
    bgeot::multi_index sizes_;

    
};// level_set_unit_normal


/**
 * Funzione che costruisce la matrice corrispondente alla forma bilineare a(u,v):  Aij = a(\phi_j, \phi_i) = A11
 */
void darcy_A11F ( sparseMatrixPtr_Type& M,
                  const FractureHandlerPtr_Type& fracture,
                  const scalar_type& gammaU,
                  const scalarVector_Type& invKTangentialInterpolated,
                  const sizeVector_Type &ExtBoundary,
                  const size_type& uncutRegionFlag );

/**
 * Funzione che aggiorna la matrice corrispondente alla forma bilineare a(u,v):  Aij = a(\phi_j, \phi_i) = A11
 * nel caso di frattura con intersezione di tipo " Cross "
 */
void darcy_A11F_Cross ( sparseMatrixPtr_Type& M,
					    const FractureHandlerPtr_Type& fracture,
					    const scalarVector_Type& invKTangentialInterpolated,
					    const FractureHandlerPtr_Type& otherFracture,
					    const size_type& cutRegionFlag );


/**
 * Funzione che costruisce la matrice tau tau corrispondente alla forma bilineare b(u,p):  Bij = b(\phi_j, \omega_i) = A12
 */
void darcy_A12F ( sparseMatrixPtr_Type& M,
                  const FractureHandlerPtr_Type& fracture,
                  const size_type& uncutRegionFlag );

/**
 * Funzione che aggiorna la matrice tau tau corrispondente alla forma bilineare b(u,p):  Bij = b(\phi_j, \omega_i) = A12
 * nel caso di frattura con intersezione di tipo " Cross "
 */
void darcy_A12F_Cross ( sparseMatrixPtr_Type& M,
                  	    const FractureHandlerPtr_Type& fracture,
                  	    const FractureHandlerPtr_Type& otherFracture,
                  	    const size_type& cutRegionFlag );


/**
 * Funzione che calcola il termine noto del sistema F(u)
 */
void darcy_dataF ( scalarVectorPtr_Type &Bstress,
                   scalarVectorPtr_Type &Bvel,
                   const BCHandlerPtr_Type& bcHandler,
                   const FractureHandlerPtr_Type& fracture,
                   const scalar_type& gammaU,
                   const scalar_type& invK,
                   const scalarVectorPtr_Type& Pneumann,
                   const scalarVectorPtr_Type& v_diri );



/**
 * Funzione che calcola il termine noto del sistema Q(p)
 */
void assembling_Source_BoundaryF ( scalarVectorPtr_Type& D,
                                   const scalarVectorPtr_Type& source,
                                   const FractureHandlerPtr_Type& fracture,
                                   const size_type& uncutRegionFlag );


/**
 * Funzione che aggiorna il termine noto del sistema Q(p) nel caso di frattura con intersezione di tipo " Cross "
 */
void assembling_SourceF ( scalarVectorPtr_Type& D,
                          const scalarVectorPtr_Type& source,
                          const FractureHandlerPtr_Type& fracture,
                          const FractureHandlerPtr_Type& otherFracture,
                          const size_type& cutRegionFlag );


/**
 * Funzione che accoppia le variabili corrispondenti alle fratture che si intersecano
 */
void coupleFractures ( sparseMatrixPtr_Type& M, const FracturesSetPtr_Type& fractures, const size_type i );


/**
 * Funzione che calcola il salto di velocità per una frattura con intersezione di tipo " Cross "
 */
void velocityJump_Cross ( sparseMatrixPtr_Type& M,
                    const FractureHandlerPtr_Type& fracture,
                    const FractureHandlerPtr_Type& otherFracture,
                    const size_type& convex );


scalarVector_Type setDOF_v ( scalarVector_Type& DOF, FracturePtrContainer_Type& Fracture, 
							 MatrixBifurcationHandler_Type& Matrix );
							 
void setDOF_v ( scalarVector_Type& DOF_p0, scalarVector_Type& DOF_v0,
				scalar_type& DOF_p1, scalar_type& DOF_v1, 
				FracturePtrContainer_Type& Fracture, 
				MatrixBifurcationHandler_Type& Matrix, 
				size_type shiftIntersect, size_type fractureNumberDOFVelocity, size_type fractureNumberDOFPressure );
							 
void setAup_i ( sparseMatrixPtr_Type& Aup_i, 
 				size_type id, size_type id_i, size_type id_j, size_type id_k, 
 				scalarVector_Type& DOF, scalarVector_Type& DOF_v, 
 				sizeVector_Type& shiftIntersect, sizeVector_Type& fractureNumberGlobalDOFVelocity,  
 				Matrix4d T, const size_type Index,
				scalar_type s = 0 );
				
void setAup_i ( sparseMatrixPtr_Type& Aup_i, size_type id, 
				FracturePtrContainer_Type& Fracture,
				scalarVector_Type& DOF_p0, scalarVector_Type& DOF_v0,
				scalar_type& DOF_p1, scalar_type& DOF_v1, 
				sizeVector_Type& shiftIntersect, sizeVector_Type& fractureNumberGlobalDOFVelocity,  
				const Matrix4d& T, const size_type Index,
				scalar_type s = 0 );
						 


} // namespace getfem

size_type GlobalIndex_Bifurcation( FracturesSetPtr_Type& M_fractures, size_type id0, size_type id1, size_type id2 );

#endif