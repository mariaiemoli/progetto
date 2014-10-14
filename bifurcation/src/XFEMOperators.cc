/** XFEMOperators.cc
 *
 * Darcy bilinear and linear forms.
 *
 * REMEMBER: the matrix corresponding to a(u,v) is Aij = a(\phi_j, \phi_i)
 *
 */

#include "../include/XFEMOperators.h"

namespace getfem
{


//matrice tau tau per la frattura
void darcy_A11F ( sparseMatrixPtr_Type& M,
                  const FractureHandlerPtr_Type& fracture,
                  const scalar_type& gammaU,
                  const scalarVector_Type& invKTangentialInterpolated,
                  const sizeVector_Type& ExtBoundary,
                  const size_type& uncutRegionFlag )
{
    const size_type shiftVelocity = fracture->getMeshFEMVelocity().nb_dof();
    const size_type shiftData = fracture->getMeshFEMPressure().nb_dof();

    sparseMatrix_Type M_;
    gmm::resize(M_, shiftVelocity, shiftVelocity);
    gmm::clear(M_);

    scalarVector_Type etaGammaUinvh(shiftData);

    for ( size_type i = 0; i < shiftData; ++i )
    {
        etaGammaUinvh [ i ] = invKTangentialInterpolated [ i ] * gammaU * fracture->getInverseMeshSize(i);
    }

    // Volume integration
    const size_type shiftMapFactor = fracture->getMagnificationMapFactor1().size();
    scalarVector_Type invF(shiftMapFactor, 0.);

    generic_assembly assem, assemGam;
    
    if ( fracture->getMeshFEMVelocity().get_qdim() == 1 )
    {
        for ( size_type i = 0; i < shiftMapFactor; ++i )
        {
            invF [ i ] = 1 / fracture->getMagnificationMapFactor1(i);
        }
        assem.set("w=data$1(#2);" "q=data$2(#2);"
            "a=comp(Base(#1).Base(#1).Base(#2).Base(#2));"
            "M(#1,#1)+=a(:, :, i,k).w(i).q(k);");
    }
    else
    {
        for ( size_type i = 0; i < shiftMapFactor; ++i )
        {
            invF [ i ] = 1 / (fracture->getMagnificationMapFactor1(i)
                    * fracture->getMagnificationMapFactor2(i));
        }
        assem.set("w=data(#2);" "q=data$2(#2);"
            "a=comp(vBase(#1).vBase(#1).Base(#2).Base(#2));"
            "M(#1,#1)+=a(:,i,:,i,j,k).w(j).q(k);");
    }

    // Assign the M_mediumMesh integration method
    assem.push_mi(fracture->getIntegrationMethodVelocity());

    // Assign the M_mediumMesh finite element space
    assem.push_mf(fracture->getMeshFEMVelocity());

    // Assign the M_mediumMesh finite element space for the coefficients
    assem.push_mf(fracture->getMeshFEMPressure());
    assem.push_mf(fracture->getMeshFEMPressure());

    // Assign the coefficients
    assem.push_data(invKTangentialInterpolated);
    assem.push_data(invF);

    // Set the matrices to save the evaluations
    assem.push_mat(M_);

    // Computes the matrices
    assem.assembly ( uncutRegionFlag );

    for ( size_type i = 0; i < shiftVelocity; ++i )
    {
        for ( size_type j = 0; j < shiftVelocity; ++j )
        {
            (*M)(i, j) = M_(i, j);
        }
    };

    cout << "DARCY :: operator a(volume)      [OK]" << endl;

    // Boundary integration for the fracture
    gmm::clear(M_);

    getfem::generic_assembly assem_surf, assem_normx, assem_normy;

    assem_surf.set("gamma=data$1(#2);"
        "t=comp(vBase(#1).Normal().vBase(#1).Normal().Base(#2));"
        "M$1(#1,#1)+=(t(:,i, i, :,j, j, k).gamma(k));");

    // Assign the M_mediumMesh integration method
    assem_surf.push_mi(fracture->getIntegrationMethodVelocity());

    // Assign the M_mediumMesh finite element space
    assem_surf.push_mf(fracture->getMeshFEMVelocity());

    // Assign the M_mediumMesh finite element space for the coefficients
    assem_surf.push_mf(fracture->getMeshFEMPressure());

    // Assign the coefficients
    assem_surf.push_data(etaGammaUinvh);

    // Set the matrices to save the evaluations
    assem_surf.push_mat(M_);

    // Assemble in each sub region
    for ( size_type bndID = 0; bndID < ExtBoundary.size(); bndID++ )
    {
        assem_surf.assembly(
                fracture->getMeshFEMVelocity().linked_mesh().get_mpi_sub_region(
                        ExtBoundary [ bndID ]));
    }

    for ( size_type i = 0; i < shiftVelocity; ++i )
    {
        for ( size_type j = 0; j < shiftVelocity; ++j )
        {
            (*M)(i, j) += M_(i, j);
        }
    }
    cout << "DARCY :: operator a(surface)     [OK]" << endl;

} // darcy_A11F


//lo stesso per la frattura
void darcy_A12F ( sparseMatrixPtr_Type& M,
                  const FractureHandlerPtr_Type& fracture,
                  const size_type& uncutRegionFlag )
{
    const size_type velocityShift = fracture->getMeshFEMVelocity().nb_dof();
    const size_type pressureShift = fracture->getMeshFEMPressure().nb_dof();

    sparseMatrix_Type M_;

    gmm::resize(M_, velocityShift, pressureShift);
    gmm::clear(M_);

    // Volume integration
    getfem::generic_assembly assem;

    if ( fracture->getMeshFEMVelocity().get_qdim() == 1 )
    {
        assem.set("M(#1,#2)+=-comp(vGrad(#1).Base(#2))" "(:, i,i,:);");
    }
    else
    {
        assem.set("M(#1,#2)+=-comp(vGrad(#1).Base(#2))" "(:,i,i, :);");
    }

    // Assign the M_mediumMesh integration method
    assem.push_mi(fracture->getIntegrationMethodVelocity());

    // Assign the M_mediumMesh finite element space
    assem.push_mf(fracture->getMeshFEMVelocity());
    assem.push_mf(fracture->getMeshFEMPressure());

    // Set the matrices to save the evaluations
    assem.push_mat(M_);

    // Computes the matrices
    assem.assembly ( uncutRegionFlag );

    for ( size_type i = 0; i < velocityShift; ++i )
    {
        for ( size_type j = 0; j < pressureShift; ++j )
        {
            (*M)(i, j) = M_(i, j);
        }
    };

    cout << "DARCY :: operator b(volumic)     [OK]" << endl;

} // darcy_A12F


//stesso lavoro con la frattura - più semplice
void darcy_dataF ( scalarVectorPtr_Type &Bstress,
                   scalarVectorPtr_Type &Bvel,
                   const BCHandlerPtr_Type& bcHandler,
                   const FractureHandlerPtr_Type& fracture,
                   const scalar_type& gammaU,
                   const scalar_type& invK,
                   const scalarVectorPtr_Type& Pneumann,
                   const scalarVectorPtr_Type& v_diri )
{

    // ----------------- Penalty    ---------------

    const scalar_type fractureID = fracture->getId();
    const LevelSetHandlerPtr_Type& levelSet = fracture->getLevelSet();
    const size_type shiftVelocity = fracture->getMeshFEMVelocity().nb_dof();
    const size_type shiftCoefficinents =
            fracture->getMeshFEMPressure().nb_dof();

    scalarVector_Type etaGammaUinvh(shiftCoefficinents, 0.), Bvel_tot(
            shiftVelocity, 0.);

    for ( size_type i = 0; i < shiftCoefficinents; ++i )
    {
        etaGammaUinvh [ i ] = invK * gammaU * fracture->getInverseMeshSize(i);
    }

    getfem::generic_assembly assem2;

    assem2.set("gamma=data$1(#2);" "vel=data$2(#2);"
        "t=comp(Base(#2).vBase(#1).Normal().Base(#2));"
        "V(#1)+=(t(m, :,j, j, k).gamma(k).vel(m));");

    // Assign the M_mediumMesh integration method
    assem2.push_mi(fracture->getIntegrationMethodVelocity());

    // Assign the M_mediumMesh finite element space
    assem2.push_mf(fracture->getMeshFEMVelocity());

    // Assign the M_mediumMesh finite element space for the coefficients
    assem2.push_mf(fracture->getMeshFEMPressure());
    assem2.push_mf(fracture->getMeshFEMPressure());

    // Assign the coefficients
    assem2.push_data(etaGammaUinvh);
    assem2.push_data(*v_diri);

    // Set the matrices to save the evaluations
    assem2.push_vec(Bvel_tot);

    // Assemble in each sub region
    const size_type shiftDirichlet =
            bcHandler->getFractureBC(fractureID)->getDirichlet().size();
    for ( size_type bndID = 0; bndID < shiftDirichlet; bndID++ )
    {
        const size_type val = bcHandler->getFractureBC(fractureID)->getDirichlet(bndID);
        assem2.assembly( fracture->getMeshFEMVelocity().linked_mesh().get_mpi_sub_region( val));
    }

    for ( size_type i = 0; i < shiftVelocity; ++i )
    {
        (*Bvel) [ i ] += Bvel_tot [ i ];
    }

    cout << "DARCY :: DATA (penal. bound.)    [OK]" << endl;

    // ----------------- Ext Stress ---------------

    const size_type shiftMapFactor1 =
            fracture->getMagnificationMapFactor1().size();
    const size_type shiftMapFactor2 =
            fracture->getMagnificationMapFactor2().size();

    scalarVector_Type Bs(shiftVelocity, 0.), coefx(shiftMapFactor1);
    scalarVector_Type coefy(shiftMapFactor2);

    getfem::generic_assembly assemb;

    if ( fracture->getMeshFEMVelocity().get_qdim() == 1 )
    {
        for ( size_type i = 0; i < shiftMapFactor1; ++i )
        {
            coefx [ i ] = 1 / fracture->getMagnificationMapFactor1(i);
        }
        assemb.set("p=data$1(#2);"
            "V(#1)+=-comp(vBase(#1).Base(#2))"
            "(:,1, h).p(h)");
    }
    else
    {

        //attenzione, quando integro sulla frattura devo tenere conto che quella vera non è flat quindi le lunghezze/aree devono essere convertite
        for ( size_type i = 0; i < shiftMapFactor1; ++i )
        {
            coefx [ i ] = 1 / fracture->getMagnificationMapFactor1(i);
            coefy [ i ] = 1 / fracture->getMagnificationMapFactor2(i);
        }
        assemb.set("p=data$1(#2);"
            "V(#1)+=-comp(vBase(#1).Normal().Base(#2))"
            "(:,k, k, h).p(h)");
    }

    // Assign the M_mediumMesh integration method
    assemb.push_mi(fracture->getIntegrationMethodVelocity());

    // Assign the M_mediumMesh finite element space
    assemb.push_mf(fracture->getMeshFEMVelocity());

    // Assign the M_mediumMesh finite element space for the coefficients
    assemb.push_mf(bcHandler->getFractureBC(fracture->getId())->getMeshFEM());

    // Assign the coefficients
    assemb.push_data(*Pneumann);

    // Set the vector to save the evaluations
    assemb.push_vec(Bs);

    // Assemble in each sub region
    const size_type shiftNeumann =
            bcHandler->getFractureBC(fractureID)->getNeumann().size();
    for ( size_type bndID = 0; bndID < shiftNeumann; bndID++ )
    {
        const size_type val = bcHandler->getFractureBC(fractureID)->getNeumann( bndID);
        assemb.assembly( fracture->getMeshFEMVelocity().linked_mesh().get_mpi_sub_region( val ));
    }

    for ( size_type i = 0; i < shiftVelocity; ++i )
    {
        (*Bstress) [ i ] = Bs [ i ];
    }

} // darcy_dataF


//termine sorgente per la frattura
void assembling_Source_BoundaryF ( scalarVectorPtr_Type& D,
                                   const scalarVectorPtr_Type& source,
                                   const FractureHandlerPtr_Type& fracture,
                                   const size_type& uncutRegionFlag )
{

    const size_type shiftPressure = fracture->getMeshFEMPressure().nb_dof();
    const size_type shiftMapFactor = fracture->getMagnificationMapFactor1().size();
    
    scalarVector_Type D_(shiftPressure, 0.0), invF(shiftMapFactor, 0);

    generic_assembly assem_Source, assem_Vx, assem_Vy;

    if ( fracture->getMeshFEMPressure().get_qdim() == 1 )
    {
        for ( size_type i = 0; i < shiftMapFactor; ++i )
        {
            invF [ i ] = 1.0 / (fracture->getMagnificationMapFactor1(i));
        }

        assem_Source.set("w=data$1(#2);" "q=data$2(#2);"
            "a=comp(Base(#1).Base(#2).Base(#2));"
            "V(#1)+=a(:, k,j).w(k).q(j)");

    }
    else
    {
        for ( size_type i = 0; i < shiftMapFactor; ++i )
        {
            invF [ i ] = 1.0 / (fracture->getMagnificationMapFactor1(i) * fracture->getMagnificationMapFactor2(i));
        }

        assem_Source.set("w=data$1(#2);" "q=data$2(#2);"
            "a=comp(Base(#1).Base(#2).Base(#2));"
            "V(#1)+=a(:, k,j).w(k).q(j)");

    }

    // Assign the M_mediumMesh integration method
    assem_Source.push_mi(fracture->getIntegrationMethodPressure());

    // Assign the M_mediumMesh finite element space
    assem_Source.push_mf(fracture->getMeshFEMPressure());

    // Assign the M_mediumMesh finite element space for the coefficients
    assem_Source.push_mf(fracture->getMeshFEMPressure());
    assem_Source.push_mf(fracture->getMeshFEMPressure());

    // Assign the coefficients
    assem_Source.push_data(*source);
    assem_Source.push_data(invF);

    // Set the vector to save the evaluations
    assem_Source.push_vec(D_);

    // Computes the matrices
    assem_Source.assembly ( uncutRegionFlag );

    gmm::add(D_, gmm::sub_vector(*D, gmm::sub_interval(0, shiftPressure)));

} // assembling_Source_BoundaryF


void coupleFractures ( sparseMatrixPtr_Type& M, const FracturesSetPtr_Type& fractures )
{
    const size_type numFractures = fractures->getNumberFractures ();
    
    const size_type numCross = fractures-> getIntersections ()-> getNumberCross ();
    const size_type numBifurcation = fractures-> getIntersections ()-> getNumberBifurcation ();

    for ( size_type i = 0; i < numFractures; ++i )
    {
        const FractureHandlerPtr_Type& fracture = fractures->getFracture(i);
        const pairSizeVectorContainer_Type& intersectElementsGlobalIndex = fracture->getFractureIntersectElementsGlobalIndex ();

        for ( size_type j = 0; j < numFractures; ++j )
        {
            const size_type numIntersections = intersectElementsGlobalIndex [j].size();
            for ( size_type k = 0; k < numIntersections; ++k )
            {
                const size_type first = intersectElementsGlobalIndex [j] [k].first;
                const size_type second = intersectElementsGlobalIndex [j] [k].second;
                
                if ( first < numCross )
                {
                	(*M)( first, first ) = 1;
                	(*M)( first, second ) = -1;
                }
                
                else if (first >= 2*numCross && first < 2*numCross + 5*numBifurcation )
                {
					(*M)( first, first ) = 1;
					(*M)( first, second ) = -1;
                	
                }
                
            }
        }
    }

} // coupleFractures



}// namespace getfem