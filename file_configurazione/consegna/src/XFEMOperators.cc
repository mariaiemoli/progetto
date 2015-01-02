
#include "../include/XFEMOperators.h"

/**************************************************************************/
/*  XFEMOperators.cc												  	  */
/*  Libreria che definisce le forme lineari e bilineari per il problema	  */
/*  di Darcy             												  */
/**************************************************************************/

namespace getfem
{

// Defining unit normal on a level set ------------------------------------

level_set_unit_normal::level_set_unit_normal ( const getfem::mesh_fem& mf_, const scalarVector_Type& U_ ) :
    mf(mf_), U(mf_.nb_basic_dof()), N(mf_.linked_mesh().dim()), gradU(1, N)
{
    sizes_.resize(1);
    sizes_ [ 0 ] = short_type(N);
    mf.extend_vector(U_, U);
}

void level_set_unit_normal::compute ( getfem::fem_interpolation_context& ctx, bgeot::base_tensor& t )
{
    size_type cv = ctx.convex_num();
    coeff.resize(mf.nb_basic_dof_of_element(cv));
    
    gmm::copy( gmm::sub_vector( U, gmm::sub_index( mf.ind_basic_dof_of_element( cv ) ) ), coeff );

    ctx.pf()->interpolation_grad( ctx, coeff, gradU, 1);
    
    scalar_type norm = gmm::vect_norm2( gmm::mat_row( gradU, 0 ));
    
    for ( size_type i = 0; i < N; ++i )
    {
        t [ i ] = gradU(0, i) / norm;
    }
    
    return;
}


//matrice A11 per la frattura non intersecata
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
    gmm::resize ( M_, shiftVelocity, shiftVelocity );
    gmm::clear ( M_ );

    scalarVector_Type etaGammaUinvh ( shiftData );

    for ( size_type i = 0; i < shiftData; ++i )
    {
        etaGammaUinvh [ i ] = invKTangentialInterpolated [ i ] * gammaU * fracture->getInverseMeshSize(i);
    }

    // Volume integration
    const size_type shiftMapFactor = fracture->getMagnificationMapFactor1().size();
    scalarVector_Type invF ( shiftMapFactor, 0. );

    generic_assembly assem;
    
    /**
     * getfem::generic_assembly assem;
     * 
     * assem.push_im(mim);
     * assem.push_mf(mf);
     * assem.push_mf(mfdata);
     * assem.push_data(F);
     * assem.push_vec(B);
     * 
     * assem.set("Z=data(#2);"
     * 			 "V(#1)+=comp(Base(#1).Base(#2))(:,j).Z(j);");
     * 
     * assem.assembly();
     * 
     * The first instructions declare the object, and set the data that it will use: a mesh_im object which holds the integration
     * methods, two mesh_fem objects, the input data F, and the destination vector B.
     * 
     * The input data is the vector F , defined on mfdata.
     * One wants to evaluate sum(j){ f_j* int_Ω (φ_i * ψ_j). The instruction must be seen as something that will be executed for each convex cv 
     * of the mesh.
     * The terms #1 and #2 refer to the first mesh_fem and the second one (i.e. mf and mfdata).
     * The instruction Z=data(#2); means that for each convex, the “tensor” Z will receive the values of the first data argument provided 
     * with push_data, at indexes corresponding to the degrees of freedom attached to the convex of the second (#2) mesh_fem 
     * (here, Z = F[mfdata.ind_dof_of_element(cv)].
     * The part V(#1)+=... means that the result of the next expression will be accumulated into the output vector (provided with push_vec).
     * Here again, #1 means that we will write the result at indexes corresponding to the degrees of freedom of the current convex with 
     * respect to the first (#1) mesh_fem.
     * 
     * The right hand side comp(Base(#1).Base(#2))(:,j).Z(j) contains two operations.
     * The first one is a computation of a tensor on the convex: comp(Base(#1).Base(#2)) is evaluated as a 2-dimensions tensor, int(φ_i*ψ_j) ,
     * for all degrees of freedom i of mf and j of mfdata attached to the current convex.
     * The next part is a reduction operation, C(:,j).Z(j): each named index (here j) is summed, i.e. the result is sum(j){ c_(i,j)*z_j }.
     * 
     * The integration method used inside comp(Base(#1).Base(#2)) is taken from mim.
     * 
     */
    if ( fracture->getMeshFEMVelocity().get_qdim() == 1 )
    {
        for ( size_type i = 0; i < shiftMapFactor; ++i )
        {
            invF [ i ] = 1 / fracture->getMagnificationMapFactor1(i);
        }
        /*
         *  definisce la forma bilineare:
         *  
         *  		a_i(u,w) = (eta_i * u, w)_L2
         *  		
         *  u velocità
         *  w funzione test per la velocità
         *  
         *  #1 velocità
         *  #2 pressione
         *  		
         */
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

    // Assegno il metodo di integrazione su M_mediumMesh 
    assem.push_mi(fracture->getIntegrationMethodVelocity());

    // Assegno lo spazio degli elementi finiti su M_mediumMesh 
    assem.push_mf(fracture->getMeshFEMVelocity());

    // Assegno lo spazio degli elementi finiti su M_mediumMesh per i coefficienti
    assem.push_mf(fracture->getMeshFEMPressure());
    assem.push_mf(fracture->getMeshFEMPressure());

    // Assegno i coefficienti
    assem.push_data(invKTangentialInterpolated);
    assem.push_data(invF);

    // Definisco la matrice dove salare i risultati
    assem.push_mat(M_);

    // Calcolo la matrice
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

    getfem::generic_assembly assem_surf;

    /*
     * tratta il termine di bordo non soggetto a condizione al contorno di dirichlet per la pressione
     * 
     */ 
    assem_surf.set("gamma=data$1(#2);"
    			   "t=comp(vBase(#1).Normal().vBase(#1).Normal().Base(#2));"		
    			   "M$1(#1,#1)+=(t(:,i, i, :,j, j, k).gamma(k));");

    // Assegno il metodo di integrazione su M_mediumMesh 
    assem_surf.push_mi(fracture->getIntegrationMethodVelocity());

    // Assegno lo spazio degli elementi finiti su M_mediumMesh 
    assem_surf.push_mf(fracture->getMeshFEMVelocity());

    // Assegno lo spazio degli elementi finiti su M_mediumMesh per i coefficienti
    assem_surf.push_mf(fracture->getMeshFEMPressure());

    // Assegno i coefficienti
    assem_surf.push_data(etaGammaUinvh);

    // Definisco la matrice dove salare i risultati
    assem_surf.push_mat(M_);

    // Assemblo la matrice su ogni sottoregione
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
    
    return;

} // darcy_A11F


//matrice tau tau per la frattura con un'intersezione di tipo Cross
void darcy_A11F_Cross ( sparseMatrixPtr_Type& M,
                  	  	const FractureHandlerPtr_Type& fracture,
                  	  	const scalarVector_Type& invKTangentialInterpolated,
                  	  	const FractureHandlerPtr_Type& otherFracture,
                  	  	const size_type& cutRegionFlag )
{
    const size_type shiftVelocity = fracture->getMeshFEMVelocity().nb_dof();
    
    sparseMatrix_Type MIn, MOut, Gamma;
    gmm::resize ( MOut, shiftVelocity, shiftVelocity );
    gmm::clear ( MOut );
    gmm::resize ( MIn, shiftVelocity, shiftVelocity );
    gmm::clear ( MIn );
    gmm::resize ( Gamma, shiftVelocity, shiftVelocity );
    gmm::clear ( Gamma );

    const size_type otherFractureId = otherFracture->getId();
    LevelSetHandlerPtr_Type& levelSetOtherFracture = otherFracture->getLevelSet();

    const size_type shiftMapFactor = fracture->getMagnificationMapFactor1().size();
    scalarVector_Type invF ( shiftMapFactor, 0. );

    for ( size_type i = 0; i < shiftMapFactor; ++i )
    {
        invF [ i ] = 1 / fracture->getMagnificationMapFactor1(i);
    }

    generic_assembly assemIn, assemOut, assemGam;

    assemIn.set("w=data(#2);" "q=data$2(#2);"
                "a=comp(Base(#1).Base(#1).Base(#2).Base(#2));"
                "M(#1,#1)+=a(:,:,j,k).w(j).q(k);");

    assemOut.set("w=data(#2);" "q=data$2(#2);"
                 "a=comp(Base(#1).Base(#1).Base(#2).Base(#2));"
                 "M(#1,#1)+=a(:,:,j,k).w(j).q(k);");

    const getfem::pintegration_method intTypeIM = getfem::int_method_descriptor ( "IM_GAUSS1D(3)" );
    getfem::mesh_im_level_set meshImLevelSetOut ( *fracture->getMeshLevelSetIntersect ( otherFractureId ),
                                                  getfem::mesh_im_level_set::INTEGRATE_OUTSIDE );

    getfem::mesh_im_level_set meshImLevelSetIn ( *fracture->getMeshLevelSetIntersect ( otherFractureId ),
                                                 getfem::mesh_im_level_set::INTEGRATE_INSIDE );

    meshImLevelSetOut.set_integration_method ( fracture->getMeshFlat().convex_index(), intTypeIM );
    meshImLevelSetIn.set_integration_method ( fracture->getMeshFlat().convex_index(), intTypeIM );

    meshImLevelSetOut.set_simplex_im ( intTypeIM );
    meshImLevelSetIn.set_simplex_im ( intTypeIM );


    // Assegno il metodo di integrazioned
    assemIn.push_mi ( meshImLevelSetIn );
    assemOut.push_mi ( meshImLevelSetOut );

    // Assegno lo spazio deglie elementi finiti per la velocità
    assemIn.push_mf ( fracture->getMeshFEMVelocity() );
    assemOut.push_mf ( fracture->getMeshFEMVelocity() );

    // Assegno lo spazio deglie elementi finiti per i coefficienti
    assemIn.push_mf ( fracture->getMeshFEMPressure() );
    assemOut.push_mf ( fracture->getMeshFEMPressure() );

    // Assegno i coefficienti
    assemIn.push_data ( invKTangentialInterpolated );
    assemIn.push_data ( invF );

    assemOut.push_data ( invKTangentialInterpolated );
    assemOut.push_data ( invF );

    assemOut.push_mat ( MOut );
    assemIn.push_mat ( MIn );

    assemOut.assembly ( cutRegionFlag );
    assemIn.assembly ( cutRegionFlag );

    // Aggiorno i gradi di libertà estesi
    const sizeVector_Type& extendedVelocity = fracture->getExtendedVelocity();
    const size_type extendedNumVelocity = fracture->getNumExtendedVelocity();

    for ( size_type i = 0; i < extendedNumVelocity; ++i )
    {
		const size_type ii = extendedVelocity [ i ];
		const base_node pointFlat = fracture->getMeshFEMVelocity().point_of_basic_dof(ii);
       	
		base_node pointMapped(0,0);
       	base_node pointMapped1(0,0);
       	scalar_type t = ii*1./(fracture->getData().getSpatialDiscretization () );
        pointMapped[0]= t;
        pointMapped1[0] = pointFlat[0];
        pointMapped1[1] = fracture->getLevelSet()->getData()->y_map( pointMapped );
	
        const scalar_type levelSetValue1 = levelSetOtherFracture->getData()->ylevelSetFunction ( pointMapped1 );

        for ( size_type j = 0; j < extendedNumVelocity; ++j )
        {
            const size_type jj = extendedVelocity [ j ];
            const base_node pointFlat = fracture->getMeshFEMVelocity().point_of_basic_dof(jj);
        
            base_node pointMapped(0,0);
            base_node pointMapped1(0,0);
            scalar_type t = jj*1./(fracture->getData().getSpatialDiscretization () );
            pointMapped[0] = t;
            pointMapped1[0] = pointFlat[0];
            pointMapped1[1] = fracture->getLevelSet()->getData()->y_map( pointMapped );

            const scalar_type levelSetValue2 = levelSetOtherFracture->getData()->ylevelSetFunction ( pointMapped1 );

            if ( levelSetValue1 < 0 && levelSetValue2 < 0 )
            {
                // i and j are both In
                (*M)(ii, jj) += MIn ( ii, jj );
                (*M)(i + shiftVelocity, j + shiftVelocity) += MOut ( ii, jj );
            }

            else if ( levelSetValue1 >= 0 && levelSetValue2 >= 0 )
            {
                // i and j are both Out
                (*M)(ii, jj) += MOut(ii, jj);
                (*M)(i + shiftVelocity, j + shiftVelocity) += MIn(ii, jj);
            }
            else if ( levelSetValue1 < 0 && levelSetValue2 >= 0 )
            {
                // i is In, j is Out
                (*M)(ii, j + shiftVelocity) += MIn(ii, jj);
                (*M)(i + shiftVelocity, jj) += MOut(ii, jj);
            }
            else
            {
                // i is Out, j is In
                (*M)(i + shiftVelocity, jj) += MIn(ii, jj);
                (*M)(ii, j + shiftVelocity) += MOut(ii, jj);
            }
        }
    }

    // Calcolo l'integrale {u dot n}{v dot n}

    scalarVector_Type otherFractureEtaNormal ( fracture->getMeshFEMPressure().nb_dof(), 0. );

    for ( size_type i = 0; i < fracture->getMeshFEMPressure().nb_dof(); ++i )
    {
        base_node nodo = fracture->getMeshFEMPressure().point_of_basic_dof( i );
        otherFractureEtaNormal [ i ] = otherFracture->getData().getEtaNormal () *
                    otherFracture->getData().etaNormalDistribution ( nodo ) / fracture->getData().getThickness();
    }

    assemGam.set("w=data(#2);"
        "a=comp(vBase(#1).NonLin(#3).vBase(#1).NonLin(#3).Base(#2));"
        "M(#1,#1)+=a(:,j,j,:,i,i,k).w(k);");

    level_set_unit_normal nterm ( fracture->getMeshLevelSetIntersect ( otherFractureId )->get_level_set(0)->get_mesh_fem(),
                                  fracture->getMeshLevelSetIntersect ( otherFractureId )->get_level_set(0)->values() );

    getfem::mesh_im_level_set meshImLevel ( *fracture->getMeshLevelSetIntersect ( otherFractureId ),
                                            getfem::mesh_im_level_set::INTEGRATE_BOUNDARY );

    meshImLevel.set_integration_method ( fracture->getMeshFlat().convex_index(), intTypeIM );

    meshImLevel.set_simplex_im ( intTypeIM );

    // Assegno il metodo di integrazione sulla mesh
    assemGam.push_mi ( meshImLevel );

    // Assegno lo spazio deglie elementi finiti
    assemGam.push_mf ( fracture->getMeshFEMVelocity() );
    assemGam.push_mf ( fracture->getMeshFEMPressure() );

    // Assegno il termine non lineare
    assemGam.push_nonlinear_term ( &nterm );

    // Assegno lo spazio deglie elementi finiti per i coefficienti
    assemGam.push_mf ( fracture->getMeshLevelSetIntersect ( otherFractureId )->get_level_set(0)->get_mesh_fem() );

    // Assegno i coefficienti
    assemGam.push_data ( otherFractureEtaNormal );

    // Set the matrices to save the evaluations
    assemGam.push_mat ( Gamma  );

    // Computes the matrices
    assemGam.assembly ( cutRegionFlag );

    // Add the extended degrees of freedom
    for ( size_type i = 0; i < extendedNumVelocity; ++i )
    {
        const size_type ii = extendedVelocity [ i ];
        for ( size_type j = 0; j < extendedNumVelocity; ++j )
        {
            const size_type jj = extendedVelocity [ j ];

            (*M)(ii, jj) += 0.25 * Gamma(ii, jj);
            (*M)(i + shiftVelocity, j + shiftVelocity) += 0.25 * Gamma(ii, jj);
            (*M)(i + shiftVelocity, jj) += 0.25 * Gamma(ii, jj);
            (*M)(ii, j + shiftVelocity) += 0.25 * Gamma(ii, jj);
        }
    }

    return;
    
} // darcy_A11F_Cross


//matrice A12 per la frattura non intersecata
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
    	/*
    	 * definisce la forma bilineare
    	 * 
    	 * 		b_i(q,w) = - (q, div(w))_L2
    	 * 	
    	 * 	w funzione test per la velocità
    	 * 	q funzione test per la pressione
    	 * 	
    	 */
        assem.set("M(#1,#2)+=-comp(vGrad(#1).Base(#2))" "(:, i,i,:);");	// ma non ci va un prodotto scalare in L2?
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
    
    return;
    
} // darcy_A12F


//lo stesso per la frattura con un'intersezione di tipo Cross 
void darcy_A12F_Cross ( sparseMatrixPtr_Type& M,
                  	    const FractureHandlerPtr_Type& fracture,
                  	    const FractureHandlerPtr_Type& otherFracture,
                  	    const size_type& cutRegionFlag )
{
    const size_type velocityShift = fracture->getMeshFEMVelocity().nb_dof();
    const size_type pressureShift = fracture->getMeshFEMPressure().nb_dof();

    sparseMatrix_Type MIn, MOut;

    gmm::resize(MIn, velocityShift, pressureShift);
    gmm::clear(MIn);
    gmm::resize(MOut, velocityShift, pressureShift);
    gmm::clear(MOut);

    const size_type otherFractureId = otherFracture->getId();
    LevelSetHandlerPtr_Type& levelSetOtherFracture = otherFracture->getLevelSet();

    // Volume integration
    getfem::generic_assembly assemIn, assemOut;

    assemIn.set("M(#1,#2)+=-comp(vGrad(#1).Base(#2))"
                "(:, i,i,:);");
    assemOut.set("M(#1,#2)+=-comp(vGrad(#1).Base(#2))"
                 "(:, i,i,:);");

    const getfem::pintegration_method intTypeIM = getfem::int_method_descriptor ( "IM_GAUSS1D(3)" );
    getfem::mesh_im_level_set meshImLevelSetOut ( *fracture->getMeshLevelSetIntersect ( otherFractureId ),
                                                  getfem::mesh_im_level_set::INTEGRATE_OUTSIDE );

    getfem::mesh_im_level_set meshImLevelSetIn ( *fracture->getMeshLevelSetIntersect ( otherFractureId ),
                                                 getfem::mesh_im_level_set::INTEGRATE_INSIDE );

    meshImLevelSetOut.set_integration_method ( fracture->getMeshFlat().convex_index(), intTypeIM );
    meshImLevelSetIn.set_integration_method ( fracture->getMeshFlat().convex_index(), intTypeIM );

    meshImLevelSetOut.set_simplex_im ( intTypeIM );
    meshImLevelSetIn.set_simplex_im ( intTypeIM );

    // Assign the mesh integration method
    assemIn.push_mi ( meshImLevelSetIn );
    assemOut.push_mi ( meshImLevelSetOut );

    // Assign the M_mediumMesh finite element space
    assemIn.push_mf ( fracture->getMeshFEMVelocity() );
    assemIn.push_mf ( fracture->getMeshFEMPressure() );

    assemOut.push_mf ( fracture->getMeshFEMVelocity() );
    assemOut.push_mf ( fracture->getMeshFEMPressure() );

    // Set the matrices to save the evaluations
    assemIn.push_mat ( MIn );
    assemOut.push_mat ( MOut );

    // Computes the matrices
    assemIn.assembly ( cutRegionFlag );
    assemOut.assembly ( cutRegionFlag );

    // Update the extended degrees of freedom
    const sizeVector_Type& extendedVelocity = fracture->getExtendedVelocity();
    const size_type extendedNumVelocity = fracture->getNumExtendedVelocity();
    const sizeVector_Type& extendedPressure = fracture->getExtendedPressure();
    const size_type extendedNumPressure = fracture->getNumExtendedPressure();

    for ( size_type i = 0; i < extendedNumVelocity; ++i )
    {
        const size_type ii = extendedVelocity [ i ];
        const base_node pointFlat = fracture->getMeshFEMVelocity().point_of_basic_dof(ii);
       
        base_node pointMapped(0,0);
        base_node pointMapped1(0,0);
        scalar_type t = ii*1./(fracture->getData().getSpatialDiscretization () );
        pointMapped[0] = t;
        pointMapped1[0] = pointFlat[0];
        pointMapped1[1] = fracture->getLevelSet()->getData()->y_map( pointMapped );

        const scalar_type levelSetValue1 = levelSetOtherFracture->getData()->ylevelSetFunction ( pointMapped1 );

        for ( size_type j = 0; j < extendedNumPressure; ++j )
        {
            const size_type jj = extendedPressure [ j ];
            const base_node pointFlat = fracture->getMeshFEMPressure().point_of_basic_dof(jj);
            
            base_node pointMapped(0,0);
            base_node pointMapped1(0,0);
            scalar_type t = jj*1./(fracture->getData().getSpatialDiscretization () );
            pointMapped[0] = t;
            pointMapped1[0] = pointFlat[0];
            pointMapped1[1] = fracture->getLevelSet()->getData()->y_map( pointMapped );
           
            const scalar_type levelSetValue2 = levelSetOtherFracture->getData()->ylevelSetFunction ( pointMapped1 );

            if ( levelSetValue1 < 0 && levelSetValue2 < 0 )
            {
                // i and j are both In
                (*M)(ii, jj) += MIn ( ii, jj );
                (*M)(i + velocityShift, j + pressureShift) += MOut ( ii, jj );
            }

            else if ( levelSetValue1 >= 0 && levelSetValue2 >= 0 )
            {
                // i and j are both Out
                (*M)(ii, jj) += MOut(ii, jj);
                (*M)(i + velocityShift, j + pressureShift) += MIn(ii, jj);
            }
            else if ( levelSetValue1 < 0 && levelSetValue2 >= 0 )
            {
                // i is In, j is Out
                (*M)(ii, j + pressureShift) += MIn(ii, jj);
                (*M)(i + velocityShift, jj) += MOut(ii, jj);
            }
            else
            {
                // i is Out, j is In
                (*M)(i + velocityShift, jj) += MIn(ii, jj);
                (*M)(ii, j + pressureShift) += MOut(ii, jj);
            }
        }
    }

    return;
    
} // darcy_A12F_Cross


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
    const size_type shiftVelocity = fracture->getMeshFEMVelocity().nb_dof();
    const size_type shiftCoefficinents = fracture->getMeshFEMPressure().nb_dof();
    
    scalarVector_Type etaGammaUinvh(shiftCoefficinents, 0.), Bvel_tot( shiftVelocity, 0.);

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
    const size_type shiftDirichlet = bcHandler->getFractureBC(fractureID)->getDirichlet().size();
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

    const size_type shiftMapFactor1 = fracture->getMagnificationMapFactor1().size();
    const size_type shiftMapFactor2 = fracture->getMagnificationMapFactor2().size();

    scalarVector_Type Bs( shiftVelocity, 0. ), coefx( shiftMapFactor1 );
    scalarVector_Type coefy( shiftMapFactor2 );

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
    const size_type shiftNeumann = bcHandler->getFractureBC(fractureID)->getNeumann().size();
    
    for ( size_type bndID = 0; bndID < shiftNeumann; bndID++ )
    {
        const size_type val = bcHandler->getFractureBC(fractureID)->getNeumann( bndID);
        assemb.assembly( fracture->getMeshFEMVelocity().linked_mesh().get_mpi_sub_region( val ));
    }

    for ( size_type i = 0; i < shiftVelocity; ++i )
    {
        (*Bstress) [ i ] = Bs [ i ];
    }

    return;
    
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

        /*
         * assembla il termine noto
         * 
         * 		sum(i) { - (f_i, q_i) }
         */
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
    
    return;

} // assembling_Source_BoundaryF


//termine sorgente per la frattura con intersezione
void assembling_SourceF ( scalarVectorPtr_Type& D,
                          const scalarVectorPtr_Type& source,
                          const FractureHandlerPtr_Type& fracture,
                          const FractureHandlerPtr_Type& otherFracture,
                          const size_type& cutRegionFlag )
{
    const size_type shiftPressure = fracture->getMeshFEMPressure().nb_dof();
    const size_type shiftMapFactor = fracture->getMagnificationMapFactor1().size();
    const size_type otherFractureId = otherFracture->getId();
    LevelSetHandlerPtr_Type& levelSetOtherFracture = otherFracture->getLevelSet();

    scalarVector_Type DIn ( shiftPressure, 0.0 ), DOut ( shiftPressure, 0.0 ),
                      invF ( shiftMapFactor, 0 );

    generic_assembly assemIn, assemOut;

    for ( size_type i = 0; i < shiftMapFactor; ++i )
    {
        invF [ i ] = 1.0 / (fracture->getMagnificationMapFactor1(i));
    }

    assemIn.set( "w=data$1(#2);" "q=data$2(#2);"
                 "a=comp(Base(#1).Base(#2).Base(#2));"
    			 "V(#1)+=a(:, k,j).w(k).q(j)" );

    assemOut.set( "w=data$1(#2);" "q=data$2(#2);"
                  "a=comp(Base(#1).Base(#2).Base(#2));"
                  "V(#1)+=a(:, k,j).w(k).q(j)" );

    const getfem::pintegration_method intTypeIM = getfem::int_method_descriptor ( "IM_GAUSS1D(3)" );
    getfem::mesh_im_level_set meshImLevelSetOut ( *fracture->getMeshLevelSetIntersect ( otherFractureId ),
                                                  getfem::mesh_im_level_set::INTEGRATE_OUTSIDE );

    getfem::mesh_im_level_set meshImLevelSetIn ( *fracture->getMeshLevelSetIntersect ( otherFractureId ),
                                                 getfem::mesh_im_level_set::INTEGRATE_INSIDE );

    meshImLevelSetOut.set_integration_method ( fracture->getMeshFlat().convex_index(), intTypeIM );
    meshImLevelSetIn.set_integration_method ( fracture->getMeshFlat().convex_index(), intTypeIM );

    meshImLevelSetOut.set_simplex_im ( intTypeIM );
    meshImLevelSetIn.set_simplex_im ( intTypeIM );

    // Assign the M_mediumMesh integration method
    assemIn.push_mi ( meshImLevelSetIn );
    assemOut.push_mi ( meshImLevelSetOut );

    // Assign the M_mediumMesh finite element space
    assemIn.push_mf ( fracture->getMeshFEMPressure() );
    assemOut.push_mf ( fracture->getMeshFEMPressure() );

    // Assign the M_mediumMesh finite element space for the coefficients
    assemIn.push_mf ( fracture->getMeshFEMPressure() );
    assemIn.push_mf ( fracture->getMeshFEMPressure() );
    assemOut.push_mf ( fracture->getMeshFEMPressure() );
    assemOut.push_mf ( fracture->getMeshFEMPressure() );

    // Assign the coefficients
    assemIn.push_data ( *source );
    assemIn.push_data ( invF );
    assemOut.push_data ( *source );
    assemOut.push_data ( invF );

    // Set the vector to save the evaluations
    assemIn.push_vec ( DIn );
    assemOut.push_vec ( DOut );

    // Computes the matrices
    assemIn.assembly ( cutRegionFlag );
    assemOut.assembly ( cutRegionFlag );

    const sizeVector_Type& extendedPressure = fracture->getExtendedPressure();
    const size_type extendedNumPressure = fracture->getNumExtendedPressure();

    for ( size_type i = 0; i < extendedNumPressure; ++i )
    {
        const size_type ii = extendedPressure [ i ];
        const base_node pointFlat = fracture->getMeshFEMPressure().point_of_basic_dof(ii);
        
        base_node pointMapped(0,0);
        base_node pointMapped1(0,0);
        scalar_type t = ii*1./(fracture->getData().getSpatialDiscretization () );
        pointMapped[0] = t;
        pointMapped1[0] = pointFlat[0];
        pointMapped1[1] = fracture->getLevelSet()->getData()->y_map( pointMapped );
        
        const scalar_type levelSetValue = levelSetOtherFracture->getData()->ylevelSetFunction ( pointMapped1 );

        if ( levelSetValue < 0 )
        {
            (*D) [ ii ] += DIn [ ii ];
            (*D) [ i + shiftPressure ] += DOut [ ii ];
        }
        else
        {
            (*D) [ ii ] += DOut [ ii ];
            (*D) [ i + shiftPressure ] += DIn [ ii ];
        }
    }

    return;
    
} // assembling_SourceF


// creiamo l'accoppiamento tra le fratture che hanno intersezione di tipo " Cross "
void coupleFractures ( sparseMatrixPtr_Type& M, const FracturesSetPtr_Type& fractures, const size_type i )
{
    const size_type numFractures = fractures->getNumberFractures ();
    
    const size_type numCross = fractures-> getIntersections ()-> getNumberCross ();
    
    const size_type numBifurcation = fractures-> getIntersections ()-> getNumberBifurcation ();
    
    const size_type numBifurcation2 = fractures-> getIntersections ()-> getNumberBifurcation2 ();
    
    if( i == 1 )
    {
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
				   
				}
			}

		}
    }
    else
    {
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
					
					if ( ( first > numCross*2 + numBifurcation ) && ( first < numCross*2 + numBifurcation + numBifurcation2 ) )
					{
						(*M)( first, first ) = 1;
						(*M)( first, second ) = -1;
					}
				   
				}
			}

		}

    }
    
    
    return;

} // coupleFractures


// Salto di velocità in presenza di intersezione di tipo " Cross "
void velocityJump_Cross ( sparseMatrixPtr_Type& M,
                    	  const FractureHandlerPtr_Type& fracture,
                    	  const FractureHandlerPtr_Type& otherFracture,
                    	  const size_type& convex )
{
    const size_type shiftVelocity = fracture->getMeshFEMVelocity().nb_dof();
    scalarVector_Type V ( shiftVelocity, 0. );

    const size_type otherFractureId = otherFracture->getId();
    LevelSetHandlerPtr_Type& levelSetOtherFracture = otherFracture->getLevelSet();

    const getfem::pintegration_method intTypeIM = getfem::int_method_descriptor ( "IM_GAUSS1D(3)" );
    getfem::mesh_im_level_set meshImLevel ( *fracture->getMeshLevelSetIntersect ( otherFractureId ),
                                            getfem::mesh_im_level_set::INTEGRATE_BOUNDARY );

    meshImLevel.set_integration_method ( fracture->getMeshFlat().convex_index(), intTypeIM );

    meshImLevel.set_simplex_im ( intTypeIM );

    getfem::mesh_region meshElement;
    meshElement.add ( convex );

    generic_assembly assem;

    assem.set ( "V(#1)+=comp(vBase(#1).NonLin(#2))"
                "(:,1,1)" );

    level_set_unit_normal nterm ( fracture->getMeshLevelSetIntersect ( otherFractureId )->get_level_set(0)->get_mesh_fem(),
                                  fracture->getMeshLevelSetIntersect ( otherFractureId )->get_level_set(0)->values() );

    assem.push_mi ( meshImLevel );

    // Assign the mesh finite element space
    assem.push_mf ( fracture->getMeshFEMVelocity() );

    // Assign the mesh finite element space for the coefficients
    assem.push_mf ( fracture->getMeshLevelSetIntersect ( otherFractureId )->get_level_set(0)->get_mesh_fem() );

    // Assign the non linear term
    assem.push_nonlinear_term ( &nterm );

    // Set the matrices to save the evaluations
    assem.push_vec ( V );

    assem.assembly ( meshElement );

    const sizeVector_Type& extendedVelocity = fracture->getExtendedVelocity();
    const size_type extendedNumVelocity = fracture->getNumExtendedVelocity();

    for ( size_type i = 0; i < extendedNumVelocity; ++i )
    {
        const size_type ii = extendedVelocity [ i ];
        const base_node pointFlat = fracture->getMeshFEMVelocity().point_of_basic_dof(ii);
        base_node pointMapped(0,0);
        base_node pointMapped1(0,0);
	scalar_type t = ii*1./(fracture->getData().getSpatialDiscretization () );
        pointMapped[0] = t;
        pointMapped1[0] = pointFlat[0];
        pointMapped1[1] = fracture->getLevelSet()->getData()->y_map( pointMapped );

	const scalar_type levelSetValue = levelSetOtherFracture->getData()->ylevelSetFunction ( pointMapped1 );

        if ( levelSetValue < 0 )
        {
            (*M) ( ii, 0 ) += V [ ii ];
            (*M) ( i + shiftVelocity, 0 ) -= V [ ii ];
        }
        else
        {
            (*M) ( ii, 0 ) -= V [ ii ];
            (*M) ( i + shiftVelocity, 0 ) += V [ ii ];

        }

    }
    
    return;

} // velocityJump_Cross


scalarVector_Type setDOF_v ( scalarVector_Type& DOF, FracturePtrContainer_Type& Fracture, 
							 MatrixBifurcationHandler_Type& Matrix )
{
	scalarVector_Type DOF_v(DOF.size());
	
		for( size_type i=0; i< DOF.size(); i++)
		{
			Matrix.SetDOFIntersecton( Fracture[ i ], DOF[ i ] );
		}
	
		for( size_type i=0;  i< DOF.size(); i++)
		{
			if( DOF[ i ] == 0 )
			{
				DOF_v[ i ] = DOF[ i ];
			}
			else
			{
				DOF_v[ i ] = DOF[ i ] + 1.;
			}
		}
	
	return DOF_v;
											   	
} //setDOF_v

void setAup_i ( sparseMatrixPtr_Type& Aup_i, 
				size_type id, size_type id_i, size_type id_j, size_type id_k, 
				scalarVector_Type& DOF, scalarVector_Type& DOF_v, 
				sizeVector_Type& shiftIntersect, sizeVector_Type& fractureNumberGlobalDOFVelocity,  
				Matrix4d T, const size_type Index,
				scalar_type s )
{						
	size_type id0, id1, id2;
	
	id0 = id_i;
	id1 = id_j;
	id2 = id_k;
	
	orderId( id0, id1, id2 );
	
	if( id != 3 )
	{
		(*Aup_i) ( 0 , shiftIntersect[ id_i ] + DOF_v[ id ] )  = 1.;
		(*Aup_i) ( 0 , shiftIntersect[ id0 ] + DOF[ 0 ] + fractureNumberGlobalDOFVelocity [ id0 ] ) = 1.*T( id , 0 );
		(*Aup_i) ( 0 , shiftIntersect[ id1 ] + DOF[ 1 ] + fractureNumberGlobalDOFVelocity [ id1 ] ) = 1.*T( id , 1 );
		(*Aup_i) ( 0 , shiftIntersect[ id2 ] + DOF[ 2 ] + fractureNumberGlobalDOFVelocity [ id2 ] ) = 1.*T( id , 2 );
		(*Aup_i) ( 0 , Index ) = -1.*( T( id , 0 ) + T( id , 1 ) + T( id , 2 ) );
	}
	else
	{
		// velocità: attenzione alla convenzione dei segni!!	
		if ( DOF_v[ 0 ] == 0 )
		{
			(*Aup_i) ( 0 , shiftIntersect[ id0 ] + DOF_v[ 0 ] )  = 1./( 3.0 * s );
		}
		else 
		{
			(*Aup_i) ( 0 , shiftIntersect[ id0 ] + DOF_v[ 0 ] )  = -1./( 3.0 * s );
		}
		if ( DOF_v[ 1 ] == 0 )
		{
			(*Aup_i) ( 0 , shiftIntersect[ id1 ] + DOF_v[ 1 ] )  = 1./( 3.0 * s );
		}
		else 
		{
			(*Aup_i) ( 0 , shiftIntersect[ id1 ] + DOF_v[ 1 ] )  = -1./( 3.0 * s );
		}
		if ( DOF_v[ 2 ] == 0 )
		{
			(*Aup_i) ( 0 , shiftIntersect[ id2 ] + DOF_v[ 2 ] )  = 1./( 3.0 * s );
		}
		else 
		{
			(*Aup_i) ( 0 , shiftIntersect[ id2 ] + DOF_v[ 2 ] )  = -1./( 3.0 * s );
		}

		// pressione 
		(*Aup_i) ( 0 , shiftIntersect[ id0 ] + DOF[ 0 ] + fractureNumberGlobalDOFVelocity [ id0 ] ) = -1./3.;
		(*Aup_i) ( 0 , shiftIntersect[ id1 ] + DOF[ 1 ] + fractureNumberGlobalDOFVelocity [ id1 ] ) = -1./3.;
		(*Aup_i) ( 0 , shiftIntersect[ id2 ] + DOF[ 2 ] + fractureNumberGlobalDOFVelocity [ id2 ] ) = -1./3.;

		// pressione media
		(*Aup_i) ( 0 , Index ) = 1.;
	}
	
	return;
											   	
} //setAup_i


}// namespace getfem

size_type GlobalIndex_Bifurcation( FracturesSetPtr_Type& M_fractures, size_type id0, size_type id1, size_type id2 )
{
    const pairSizeVectorContainer_Type& intersectElementsGlobalIndex0 = M_fractures->getFracture( id0 )->getFractureIntersectElementsGlobalIndex ();
    
	const size_type globalIndex01 =  intersectElementsGlobalIndex0[ id1 ][ 0 ].first;
	const size_type globalIndex02 =  intersectElementsGlobalIndex0[ id2 ][ 0 ].first;

    const pairSizeVectorContainer_Type& intersectElementsGlobalIndex1 = M_fractures->getFracture( id1 )->getFractureIntersectElementsGlobalIndex ();
    
	const size_type globalIndex10 =  intersectElementsGlobalIndex1[ id0 ][ 0 ].first;
	const size_type globalIndex12 =  intersectElementsGlobalIndex1[ id2 ][ 0 ].first;
	
	const pairSizeVectorContainer_Type& intersectElementsGlobalIndex2 = M_fractures->getFracture( id2 )->getFractureIntersectElementsGlobalIndex ();
	
	const size_type globalIndex20 =  intersectElementsGlobalIndex2[ id0 ][ 0 ].first;
	const size_type globalIndex21 =  intersectElementsGlobalIndex2[ id1 ][ 0 ].first;

	const size_type globalIndex = fmin ( globalIndex01, fmin( globalIndex02, fmin ( globalIndex10, fmin ( globalIndex12, fmin ( globalIndex20, globalIndex21 )))));
	
	return globalIndex;
	
} // GlobalIndex_Bifurcation