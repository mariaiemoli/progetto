/** 
 * 
 * XFEMOperators.h
 * 
 * Darcy bilinear and linear forms.
 *
 * REMEMBER: the matrix corresponding to a(u,v) is Aij = a(\phi_j, \phi_i)
 *
 */

#ifndef _DARCY_OPERATORSX_
#define _DARCY_OPERATORSX_ 1

#include "Core.h"
#include "LevelSetHandler.h"
#include "MeshHandler.h"
#include "FractureHandler.h"
#include "BCHandler.h"

namespace getfem
{

class level_set_unit_normal : public getfem::nonlinear_elem_term
{
    const getfem::mesh_fem& mf;
    scalarVector_Type U;
    size_type N;
    base_matrix gradU;
    bgeot::base_vector coeff;
    bgeot::multi_index sizes_;

public:

    level_set_unit_normal ( const getfem::mesh_fem& mf_,
                            const scalarVector_Type& U_ );

    inline const bgeot::multi_index& sizes ( ) const
    {
        return sizes_;
    }

    virtual void compute ( getfem::fem_interpolation_context& ctx,
                           bgeot::base_tensor& t );
};



//matrice tau tau per la frattura
void darcy_A11F ( sparseMatrixPtr_Type& M,
                  const FractureHandlerPtr_Type& fracture,
                  const scalar_type& gammaU,
                  const scalarVector_Type& invKTangentialInterpolated,
                  const sizeVector_Type &ExtBoundary,
                  const size_type& uncutRegionFlag );

//matrice tau tau per la frattura con un'intersezione di tipo Cross
void darcy_A11F_Cross ( sparseMatrixPtr_Type& M,
					    const FractureHandlerPtr_Type& fracture,
					    const scalarVector_Type& invKTangentialInterpolated,
					    const FractureHandlerPtr_Type& otherFracture,
					    const size_type& cutRegionFlag );


//matrice tau tau per la frattura con un'intersezione di tipo Bifurcation
void darcy_A11F_Bifurcation ( sparseMatrixPtr_Type& M,
							  const FractureHandlerPtr_Type& fracture,
							  const scalarVector_Type& invKTangentialInterpolated,
							  const FractureHandlerPtr_Type& otherFracture,
							  const size_type& cutRegionFlag );


//lo stesso per la frattura
void darcy_A12F ( sparseMatrixPtr_Type& M,
                  const FractureHandlerPtr_Type& fracture,
                  const size_type& uncutRegionFlag );


//lo stesso per la frattura intersecata con un'intersezione di tipo Cross
void darcy_A12F_Cross ( sparseMatrixPtr_Type& M,
                  	    const FractureHandlerPtr_Type& fracture,
                  	    const FractureHandlerPtr_Type& otherFracture,
                  	    const size_type& cutRegionFlag );


//lo stesso per la frattura intersecata con un'intersezione di tipo Bifurcation
void darcy_A12F_Bifurcation ( sparseMatrixPtr_Type& M,
                  	  	  	  const FractureHandlerPtr_Type& fracture,
                  	  	  	  const FractureHandlerPtr_Type& otherFracture,
							  const size_type& cutRegionFlag );


//stesso lavoro con la frattura - più semplice
void darcy_dataF ( scalarVectorPtr_Type &Bstress,
                   scalarVectorPtr_Type &Bvel,
                   const BCHandlerPtr_Type& bcHandler,
                   const FractureHandlerPtr_Type& fracture,
                   const scalar_type& gammaU,
                   const scalar_type& invK,
                   const scalarVectorPtr_Type& Pneumann,
                   const scalarVectorPtr_Type& v_diri );


//termine sorgente per la frattura
void assembling_Source_BoundaryF ( scalarVectorPtr_Type& D,
                                   const scalarVectorPtr_Type& source,
                                   const FractureHandlerPtr_Type& fracture,
                                   const size_type& uncutRegionFlag );


//termine sorgente per la frattura
void assembling_SourceF ( scalarVectorPtr_Type& D,
                          const scalarVectorPtr_Type& source,
                          const FractureHandlerPtr_Type& fracture,
                          const FractureHandlerPtr_Type& otherFracture,
                          const size_type& cutRegionFlag );


//funzione che accoppia le variabili corrispondenti alle fratture che si intersecano
void coupleFractures ( sparseMatrixPtr_Type& M, const FracturesSetPtr_Type& fractures );


void velocityJump_Cross ( sparseMatrixPtr_Type& M,
                    const FractureHandlerPtr_Type& fracture,
                    const FractureHandlerPtr_Type& otherFracture,
                    const size_type& convex );


void velocityJump_Bifurcation ( sparseMatrixPtr_Type& M,
                    			const FractureHandlerPtr_Type& fracture,
                    			const FractureHandlerPtr_Type& otherFracture,
                    			const size_type& convex );


} // namespace getfem

#endif