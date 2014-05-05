#include "../include/DarcyFractured.h"

DarcyFractured::DarcyFractured ( const MediumDataPtr_Type& medium,
                                 const MeshHandlerPtr_Type& mesh,
                                 const BCHandlerPtr_Type& bcHandler,
                                 const FracturesSetPtr_Type& fractures,
                                 const ExporterPtr_Type& exporter ) :
            M_mediumData(medium), M_mesh(mesh), M_bcHandler(bcHandler), M_fractures(
            fractures), M_exporter(exporter),M_fractureEtaNormalOnMedium(
            M_fractures->getNumberFractures()), M_fractureVelocity(M_fractures->getNumberFractures()),
            M_fracturePressure(M_fractures->getNumberFractures())
{}

/* Build the M_mediumMesh, set finite element
 * and integration methods and select the boundaries.
 */
void DarcyFractured::init ( )
{

    // Allocate data
    const size_type shifCoefficients =
            M_mesh->getMeshFEMCoefficients().nb_dof();

    const size_type numberFractures = M_fractures->getNumberFractures();

    M_mediumEtaInterpolated.resize(3);

    // Allocate the vector for M_mediumEtaInterpolated
    gmm::resize(M_mediumEtaInterpolated[0], shifCoefficients);
    gmm::clear(M_mediumEtaInterpolated[0]);

    gmm::resize(M_mediumEtaInterpolated[1], shifCoefficients);
    gmm::clear(M_mediumEtaInterpolated[1]);

    gmm::resize(M_mediumEtaInterpolated[2], shifCoefficients);
    gmm::clear(M_mediumEtaInterpolated[2]);

    for ( size_type f = 0; f < numberFractures; ++f )
    {
        // Allocate the vector for M_fracturesEtaNormalOnMedium
        gmm::resize(M_fractureEtaNormalOnMedium [ f ], shifCoefficients);
        gmm::clear(M_fractureEtaNormalOnMedium [ f ]);

    }

    // Fill the vectors M_mediumEtaInterpolated, M_mediumMuInterpolated, M_fracturesEtaNormalOnMedium and M_fracturesMuNormalOnMedium of the bulk
    for ( size_type i = 0; i < shifCoefficients; ++i )
    {
        const base_node node = M_mesh->getMeshFEMCoefficients().point_of_dof(i);
        M_mediumEtaInterpolated[0] [ i ] = M_mediumData->invKDistribution11(node)
                * M_mediumData->getInvK();

        M_mediumEtaInterpolated[1] [ i ] = M_mediumData->invKDistribution12(node)
                * M_mediumData->getInvK();

        M_mediumEtaInterpolated[2] [ i ] = M_mediumData->invKDistribution22(node)
                * M_mediumData->getInvK();

        for ( size_type f = 0; f < numberFractures; ++f )
        {
            M_fractureEtaNormalOnMedium [ f ] [ i ]
                    = M_fractures->getFracture( f )->getData().etaNormalDistribution(node)
                      * M_fractures->getFracture ( f )->getData().getEtaNormal();

        }
    }

} // init


// Assemble all the matrices and the right hand side for the Darcy and a piece of the matrix for the transport
void DarcyFractured::assembly ( )
{

    /*
     * Our system has the following structure (DARCY):
     *
     * [  A11   A12 ] [ V ]  = [ Bv ]
     * [ -A12'  0   ] [ P ]  = [ Bp ]
     *
     * (V,P = velocity, pressure)
     *
     * where
     *
     * BV = Mvd * Vdirichlet + Bstress
     * BP = Mpd * Vdirichlet
     */

    const scalar_type numberFractures = M_fractures->getNumberFractures();

  
    // Degrees of freedom of the fracture primal
    sizeVector_Type fractureNumberDOFPressure(numberFractures);

    sizeVector_Type fractureNumberGlobalDOFPressure(numberFractures);

    // Degrees of freedom of the fracture dual
    sizeVector_Type fractureNumberDOFVelocity(numberFractures);

    sizeVector_Type fractureNumberGlobalDOFVelocity(numberFractures);

    // Global number of Degrees of freedom of the fracture
    sizeVector_Type fractureNumberDOFVelocityPressure(numberFractures);

    // Global number of Degrees of freedom of the fracture
    sizeVector_Type fractureNumberBoundaryDOF(numberFractures);

    // Total number of degrees of freedom
    size_type fractureTotalNumberDOFVelocityPressure = 0;

    // Total number of degrees of freedom for the pressure/concentration in the fractures
    size_type fractureTotalNumberDOFPressure = 0;

    size_type fractureNumberIntersections = 0;

    for ( size_type f = 0; f < numberFractures; ++f )
    {
        fractureNumberDOFPressure [ f ]
                = M_fractures->getFracture( f )->getMeshFEMPressure().nb_dof();

        fractureNumberGlobalDOFPressure [ f ] =
            fractureNumberDOFPressure [ f ] + M_fractures->getFracture( f )->getNumExtendedPressure();

        fractureTotalNumberDOFPressure += fractureNumberGlobalDOFPressure [ f ];

        fractureNumberDOFVelocity [ f ]
                = M_fractures->getFracture( f )->getMeshFEMVelocity().nb_dof();

        fractureNumberGlobalDOFVelocity [ f ]
                = fractureNumberDOFVelocity [ f ] + M_fractures->getFracture( f )->getNumExtendedVelocity();

        fractureNumberDOFVelocityPressure [ f ]
                = fractureNumberGlobalDOFVelocity [ f ]
                        + fractureNumberGlobalDOFPressure [ f ];

        fractureTotalNumberDOFVelocityPressure
                += fractureNumberDOFVelocityPressure [ f ];

        fractureNumberBoundaryDOF [ f ]
                = M_bcHandler->getFractureBC(f)->getMeshFEM().nb_dof();

        fractureNumberIntersections += M_fractures->getFracture( f )->getNumIntersections();
    }

  

    // Initialize all the block matrices, global matrices and right hand side for the system and the solution vector

    // Allocate the global matrix M_darcyGlobalMatrix
    M_globalMatrix.reset(new sparseMatrix_Type(
            fractureTotalNumberDOFVelocityPressure + fractureNumberIntersections,
            fractureTotalNumberDOFVelocityPressure + fractureNumberIntersections ));
    gmm::clear(*M_globalMatrix);

    // Allocate the global vector M_darcyGlobalRightHandSide
    M_globalRightHandSide.reset(new scalarVector_Type(
             fractureTotalNumberDOFVelocityPressure + fractureNumberIntersections ));
    gmm::clear(*M_globalRightHandSide);

    // Allocate the global vector M_darcyVelocityAndPressure
    M_velocityAndPressure.reset(new scalarVector_Type(
             fractureTotalNumberDOFVelocityPressure + fractureNumberIntersections ));
    gmm::clear(*M_velocityAndPressure);

  
    // Blocks of matrix
    sparseMatrixPtr_Type  App;

    

    App.reset(new sparseMatrix_Type(fractureNumberIntersections/2,
            fractureNumberIntersections));
    gmm::clear(*App);

    // Block of fracture matrix
    sparseMatrixPtrContainer_Type A11F(numberFractures), A12F(numberFractures),
            E(numberFractures);

    for ( size_type f = 0; f < numberFractures; ++f )
    {
        // Allocate the fracture matrix A11F
        A11F [ f ].reset(new sparseMatrix_Type ( fractureNumberGlobalDOFVelocity [ f ],
                fractureNumberGlobalDOFVelocity [ f ]) );
        gmm::clear(*(A11F [ f ]));

        // Allocate the fracture matrix A12F
        A12F [ f ].reset(new sparseMatrix_Type ( fractureNumberGlobalDOFVelocity [ f ],
                fractureNumberGlobalDOFPressure [ f ]) );
        gmm::clear(*(A12F [ f ]));

     
    }

    // Assembling phase

    getfem::coupleFractures ( App, M_fractures );
    size_type shiftIntersect = 0;

    for ( size_type f = 0; f < numberFractures; ++f )
    {
        std::cout << "Fracture " << f << std::endl;

        
        // Computes the matrix \int_\gamma \tau_i \cdot \tau_j
        getfem::darcy_A11F ( A11F [ f ], M_fractures->getFracture( f ),
                             M_mediumData->getPenaltyVector(),
                             M_fractures->getFracture( f )->getEtaTangentialInterpolated(),
                             M_bcHandler->getFractureBC(f)->getDirichlet(), FractureHandler::FRACTURE_UNCUT * ( f + 1 ) );

        // Computes the matrix \int_\gamma \nabla v_i \cdot \tau_j
        getfem::darcy_A12F ( A12F [ f ], M_fractures->getFracture( f ), FractureHandler::FRACTURE_UNCUT * ( f + 1 ) );

        sparseMatrixPtr_Type Aup;
        Aup.reset ( new sparseMatrix_Type ( fractureNumberGlobalDOFVelocity [ f ], 1 ) );
        gmm::clear(*Aup);

        for ( size_type otherFracture = 0; otherFracture < numberFractures; ++otherFracture )
        {
            if ( M_fractures->getFracture( f )->getMeshLevelSetIntersect ( otherFracture ).get() )
            {
                getfem::darcy_A11F ( A11F [ f ], M_fractures->getFracture( f ),
                                     M_fractures->getFracture( f )->getEtaTangentialInterpolated(),
                                     M_fractures->getFracture ( otherFracture ),
                                     FractureHandler::FRACTURE_INTERSECT * ( f + 1 ) + otherFracture + 1 );

                getfem::darcy_A12F ( A12F [ f ], M_fractures->getFracture( f ), M_fractures->getFracture ( otherFracture ),
                                     FractureHandler::FRACTURE_INTERSECT * ( f + 1 ) + otherFracture + 1 );

                const pairSizeVectorContainer_Type& intersectElementsGlobalIndex =
                                         M_fractures->getFracture( f )->getFractureIntersectElementsGlobalIndex ();
                const size_type numIntersections = intersectElementsGlobalIndex [otherFracture].size();
                const sizeVectorContainer_Type& intersectElements = M_fractures->getFracture( f )->getFractureIntersectElements ();
                for ( size_type k = 0; k < numIntersections; ++k )
                {
                    gmm::clear (*Aup);

                    getfem::velocityJump ( Aup, M_fractures->getFracture( f ), M_fractures->getFracture ( otherFracture ),
                                           intersectElements[otherFracture][k] );

                    const size_type globalIndex = intersectElementsGlobalIndex [otherFracture] [k].first;
                    const size_type globalIndex2 = intersectElementsGlobalIndex [otherFracture] [k].second;

                    gmm::copy ( *Aup, gmm::sub_matrix (*M_globalMatrix,
                        gmm::sub_interval ( shiftIntersect, fractureNumberGlobalDOFVelocity [ f ] ),
                        gmm::sub_interval (  fractureTotalNumberDOFVelocityPressure + globalIndex, 1 ) ) );

                    gmm::copy ( gmm::transposed(*Aup), gmm::sub_matrix (*M_globalMatrix,
                        gmm::sub_interval (  fractureTotalNumberDOFVelocityPressure +
                        std::max(globalIndex, globalIndex2), 1 ),
                        gmm::sub_interval ( shiftIntersect, fractureNumberGlobalDOFVelocity [ f ] ) ) );

                }

            }
        }

        shiftIntersect += fractureNumberDOFVelocityPressure [ f ];
        //massLumping(*(A11F [ f ]));
//FINE MODIFICATO

    }



    // Copy blocks into the system matrix M_darcyGlobalMatrix
    //     [  A11  A12  0          ]
    // M = [ -A12  A22  0          ]
    //     [  0    0    A11F  A12F ]
    //     [  0    0   -A12F  0    ]


 
    // Shift for the fracture
    size_type fractureShift = 0;
    for ( size_type f = 0; f < numberFractures; ++f )
    {

        // Copy the matrix A11F in M_darcyGlobalMatrix in the correct position
        gmm::copy(*(A11F [ f ]),
                gmm::sub_matrix(*M_globalMatrix, gmm::sub_interval(
                        fractureShift, fractureNumberGlobalDOFVelocity [ f ]),
                        gmm::sub_interval(fractureShift,
                                fractureNumberGlobalDOFVelocity [ f ])));

        // Copy the matrix A12F in M_darcyGlobalMatrix in the correct position
        gmm::copy(*(A12F [ f ]),
                gmm::sub_matrix(*M_globalMatrix, gmm::sub_interval(
                        fractureShift, fractureNumberGlobalDOFVelocity [ f ]),
                        gmm::sub_interval(fractureShift
                                + fractureNumberGlobalDOFVelocity [ f ],
                                fractureNumberGlobalDOFPressure [ f ])));

        // Copy the matrix -A12F in M_darcyGlobalMatrix in the correct position
        gmm::copy(gmm::transposed(gmm::scaled(*(A12F [ f ]), -1.0)),
                gmm::sub_matrix(*M_globalMatrix, gmm::sub_interval(
                        fractureShift + fractureNumberGlobalDOFVelocity [ f ],
                        fractureNumberGlobalDOFPressure [ f ]), gmm::sub_interval(
                        fractureShift, fractureNumberGlobalDOFVelocity [ f ])));

   
        // Update the shift
        fractureShift += fractureNumberDOFVelocityPressure [ f ];

 
    }

   
    gmm::copy(*App, gmm::sub_matrix(*M_globalMatrix, gmm::sub_interval(
         fractureTotalNumberDOFVelocityPressure,
        fractureNumberIntersections/2 ), gmm::sub_interval(
         fractureTotalNumberDOFVelocityPressure,
        fractureNumberIntersections ) ));
 

    // Vector, Neumann conditions
    scalarVectorPtrContainer_Type BstressF(numberFractures);

    // Pressure term (Neumann)
    scalarVectorPtrContainer_Type PneumannF(numberFractures);

    // Velocity term (Dirichilet)
    scalarVectorPtrContainer_Type VdirichletF(numberFractures);

    // Right Hand Side: velocity and pressure
    scalarVectorPtrContainer_Type B_vF(numberFractures);

    // Right Hand Side: velocity and pressure
    scalarVectorPtrContainer_Type B_pF(numberFractures);

    for ( size_type f = 0; f < numberFractures; ++f )
    {
        BstressF [ f ].reset(new scalarVector_Type(
                fractureNumberGlobalDOFVelocity [ f ], 0.));

        PneumannF [ f ].reset(new scalarVector_Type(
                fractureNumberBoundaryDOF [ f ], 0.));

        VdirichletF [ f ].reset(new scalarVector_Type(
                fractureNumberDOFPressure [ f ], 0.));

        B_vF [ f ].reset(new scalarVector_Type(fractureNumberGlobalDOFVelocity [ f ],
                0.));

        B_pF [ f ].reset(new scalarVector_Type(fractureNumberGlobalDOFPressure [ f ],
                0.));
    }

  
    for ( size_type f = 0; f < numberFractures; ++f )
    {
        const BCPtr_Type& fractureBC = M_bcHandler->getFractureBC(f);

        for ( size_type i = 0; i < fractureNumberBoundaryDOF [ f ]; i++ )
        {
            const base_node node = fractureBC->getMeshFEM().point_of_basic_dof(
                    i);

            (*(PneumannF [ f ])) [ i ]
                    = M_fractures->getFracture( f )->getData().pressureExact(node);
        }

        (*(PneumannF [ f ])) [ 0 ] *= -1; //perché la normale ha il segno meno all'inizio di una roba 1D
    }

 

    fractureShift = 0;
    for ( size_type f = 0; f < numberFractures; ++f )
    {

        fractureShift += M_mesh->getExtendedDOFVector(f).size();

        // Computes the right hand side rhs which stores the boundary conditions for the fracture
        getfem::darcy_dataF ( BstressF [ f ], B_vF [ f ], M_bcHandler,
                              M_fractures->getFracture( f ), M_mediumData->getPenaltyVector(),
                              M_mediumData->getInvK(), PneumannF [ f ], VdirichletF [ f ] );
    }

    fractureShift = 0;
    for ( size_type f = 0; f < numberFractures; ++f )
    {
        for ( size_type i = 0; i < fractureNumberGlobalDOFVelocity [ f ]; ++i )
        {
            (*M_globalRightHandSide) [ fractureShift + i ]
                    += (*(BstressF [ f ])) [ i ];

            (*M_globalRightHandSide) [ fractureShift + i ]
                    += (*(B_vF [ f ])) [ i ];
        }

        // Update the shift
        fractureShift += fractureNumberDOFVelocityPressure [ f ];
    }

   

    scalarVectorPtrContainer_Type divF(numberFractures);

    fractureShift = 0;
    for ( size_type f = 0; f < numberFractures; ++f )
    {
        divF [ f ].reset(new scalarVector_Type(fractureNumberDOFPressure [ f ],
                0.));

        for ( size_type i = 0; i < fractureNumberDOFPressure [ f ]; i++ )
        {
            const base_node
                    node =
                            M_fractures->getFracture( f )->getMeshFEMPressure().point_of_basic_dof(
                                    i);

            (*(divF [ f ])) [ i ] = M_fractures->getFracture( f )->getData().darcySource(
                    node);
        }

        gmm::clear(*(B_vF [ f ]));
        gmm::clear(*(B_pF [ f ]));

        getfem::assembling_Source_BoundaryF ( B_pF [ f ], divF [ f ],
                M_fractures->getFracture( f ), FractureHandler::FRACTURE_UNCUT * ( f + 1 ) );

        for ( size_type otherFracture = 0; otherFracture < numberFractures; ++otherFracture )
        {
            if ( M_fractures->getFracture( f )->getMeshLevelSetIntersect ( otherFracture ).get() )
            {
                getfem::assembling_SourceF ( B_pF [ f ], divF [ f ], M_fractures->getFracture( f ),
                                             M_fractures->getFracture ( otherFracture ),
                                             FractureHandler::FRACTURE_INTERSECT * ( f + 1 ) + otherFracture + 1 );
            }
        }

        for ( size_type i = 0; i < fractureNumberGlobalDOFVelocity [ f ]; ++i )
        {
            (*M_globalRightHandSide) [ fractureShift + i ]
                    += (*(B_vF [ f ])) [ i ];
        }

        for ( size_type i = fractureNumberGlobalDOFVelocity [ f ]; i
                < fractureNumberDOFVelocityPressure [ f ]; ++i )
        {
                (*M_globalRightHandSide) [ fractureShift + i ]
                    += (*(B_pF [ f ])) [ i - fractureNumberGlobalDOFVelocity [ f ] ];
        }

        // Update the shift
        fractureShift += fractureNumberDOFVelocityPressure [ f ];
    }

  
  
    // Curo la matrice
    for ( size_type i = 0; i < M_globalRightHandSide->size(); ++i )
    {
        scalar_type somma = 0;

        for ( size_type j = 0; j < M_globalRightHandSide->size(); ++j )
        {
            somma += std::fabs( (*M_globalMatrix) (i, j) );
        }

        if ( somma == 0 )
        {
            (*M_globalRightHandSide) [ i ] = 0.;
            (*M_globalMatrix) (i, i) = 1.;
        }
    }

    M_exporter->spy(M_globalMatrix, "./matrice.mm");
    M_exporter->spy(M_globalRightHandSide, "./rhs.mm");

} // assembly

// Solve the Darcy for the governing flux and do the time loop for the evolution problem
void DarcyFractured::solve ( )
{  
    const scalar_type numberFractures = M_fractures->getNumberFractures();
    const scalar_type numberIntersections = M_fractures->getIntersections()
                                            ->getNumberIntersections();
   
    sizeVector_Type fractureNumberDOFVelocity(numberFractures);
    sizeVector_Type fractureNumberGlobalDOFVelocity(numberFractures);

    sizeVector_Type fractureNumberDOFPressure(numberFractures);
    sizeVector_Type fractureNumberGlobalDOFPressure(numberFractures);

    sizeVector_Type fractureNumberDOFVelocityPressure(numberFractures);

    size_type fractureTotalNumberDOFVelocityPressure(0);

    size_type fractureNumberIntersections = 0;

    for ( size_type f = 0; f < numberFractures; ++f )
    {
        fractureNumberDOFVelocity [ f ]
                = M_fractures->getFracture( f )->getMeshFEMVelocity().nb_dof();

        fractureNumberGlobalDOFVelocity [ f ]
                = fractureNumberDOFVelocity [ f ] + M_fractures->getFracture( f )->getNumExtendedVelocity();

        fractureNumberDOFPressure [ f ]
                = M_fractures->getFracture( f )->getMeshFEMPressure().nb_dof();

        fractureNumberGlobalDOFPressure [ f ] =
                fractureNumberDOFPressure [ f ] + M_fractures->getFracture( f )->getNumExtendedPressure();

        fractureNumberDOFVelocityPressure [ f ]
                = fractureNumberGlobalDOFVelocity [ f ]
                        + fractureNumberGlobalDOFPressure [ f ];

        fractureTotalNumberDOFVelocityPressure
                += fractureNumberDOFVelocityPressure [ f ];

        M_fracturePressure [ f ].reset(new scalarVector_Type(
                fractureNumberGlobalDOFPressure [ f ], 0));

        M_fractureVelocity [ f ].reset(new scalarVector_Type(
                fractureNumberGlobalDOFVelocity [ f ], 0));

        fractureNumberIntersections += M_fractures->getFracture( f )->getNumIntersections();

    }

   

    // Solve the Darcy problem
    std::cout << std::endl << "Solving problem in Omega..." << std::flush;
    scalar_type roundConditionNumber;
    int ris;
    ris=SuperLU_solve(*M_globalMatrix, *M_velocityAndPressure,
                  *M_globalRightHandSide, roundConditionNumber);
    std::cout << " completed!" << std::endl;

 
    size_type fractureShift = 0;
    for ( size_type f = 0; f < numberFractures; ++f )
    { 
        // Extract the dual in the fracture
        gmm::copy(gmm::sub_vector(*M_velocityAndPressure, gmm::sub_interval(
                fractureShift, fractureNumberGlobalDOFVelocity [ f ])),
                *(M_fractureVelocity [ f ]));

        // Extract the primal in the fracture
        gmm::copy(gmm::sub_vector(*M_velocityAndPressure, gmm::sub_interval(
                fractureShift + fractureNumberGlobalDOFVelocity [ f ],
                fractureNumberGlobalDOFPressure [ f ])), *(M_fracturePressure [ f ]));

        // Update the shift
        fractureShift += fractureNumberDOFVelocityPressure [ f ];
    }

    // the extra dof start at the end of the fracture dofs
    size_type intersectShift = 0;
  
    getfem::pfem fractureFETypePressure;
    if ( numberFractures > 0 )
    {
        fractureFETypePressure = getfem::fem_descriptor(
                M_fractures->getFracture( 0 )->getData().getFEMTypePressure());
    }


    //--------export solution in the fracture-------------------
    for ( size_type f = 0; f < numberFractures; ++f )
    {
        FractureHandlerPtr_Type& fracture = M_fractures->getFracture(f);
        std::ostringstream osFileName;

        getfem::mesh_level_set meshFLevelSetCutFlat ( fracture->getMeshFlat() );
        for ( size_type otherF = 0; otherF < numberFractures; ++otherF )
        {
            GFLevelSetPtr_Type levelSetPtr = fracture->getLevelSetIntersect ( otherF );

            if ( levelSetPtr.get() != NULL )
            {
                meshFLevelSetCutFlat.add_level_set ( *levelSetPtr );
            }
        }

        meshFLevelSetCutFlat.adapt();

        getfem::mesh meshFcutFlat;
        getfem::mesh meshFcutMapped;
        meshFLevelSetCutFlat.global_cut_mesh ( meshFcutFlat );

        getfem::mesh_fem meshFEMcutFlat ( meshFcutFlat, fracture->getMeshFEMPressure().get_qdim());
        meshFEMcutFlat.set_finite_element ( fracture->getMeshFEMVelocity().fem_of_element(0) );

        scalarVector_Type ordinataUncut ( fractureNumberDOFVelocity [ f ], 0. );
        scalarVector_Type ordinataCut ( meshFEMcutFlat.nb_dof(), 0. );
        for ( size_type i = 0; i < fractureNumberDOFVelocity [ f ]; ++i )
        {
            const bgeot::base_node P = fracture->getMeshFEMVelocity().point_of_basic_dof(i);
            ordinataUncut [ i ] = fracture->getLevelSet()->getData()->z_map(P);
        }

        getfem::interpolation ( fracture->getMeshFEMVelocity(), meshFEMcutFlat, ordinataUncut, ordinataCut );

        for ( size_type i = 0; i < meshFEMcutFlat.nb_dof(); ++i )
        {
            bgeot::base_node P ( fracture->getData().getSpaceDimension() + 2 );
            P [ 0 ] = meshFEMcutFlat.point_of_basic_dof(i)[0];
            P [ 1 ] = ordinataCut [ i ];
            meshFcutMapped.add_point(P);
        }

        const size_type numConvex = meshFcutFlat.convex_index().size();
        for ( size_type i = 0; i < numConvex; ++i )
        {
            std::vector<bgeot::size_type> point(3);
            point [ 0 ] = meshFcutFlat.ind_points_of_convex(i) [ 0 ];
            point [ 1 ] = meshFcutFlat.ind_points_of_convex(i) [ 1 ];
            point [ 2 ] = meshFcutFlat.ind_points_of_convex(i) [ 2 ];
            meshFcutMapped.add_convex( fracture->getGeometricTransformation(), point.begin() );
        }

        osFileName << "cmeshF2-" << f << ".vtk";
        exportMesh(M_exporter->getFolder() + osFileName.str().c_str(), meshFcutMapped );

        getfem::mesh_fem mfproj ( meshFcutMapped, fracture->getMeshFEMPressure().get_qdim() );

        mfproj.set_classical_discontinuous_finite_element ( 0, 0.01 );

        scalarVector_Type fracturePressureMeanUNCUT ( fractureNumberDOFPressure [ f ], 0. );
        scalarVector_Type fracturePressureInCut ( fractureNumberDOFPressure [ f ], 0. );
        scalarVector_Type fracturePressureOutCut ( fractureNumberDOFPressure [ f ], 0. );

        for ( size_type i = 0; i < fractureNumberDOFPressure [ f ]; ++i )
        {
            const size_type el =  fracture->getMeshFEMPressure().first_convex_of_basic_dof(i);

            if ( fracture->getMeshFlat().region ( FractureHandler::FRACTURE_UNCUT * ( f + 1 ) ).is_in(el) )
            {
                fracturePressureMeanUNCUT [ i ] = (*(M_fracturePressure [ f ])) [ i ];
            }
        }

        scalarVector_Type fracturePressureMeanUNCUTInterpolated ( mfproj.nb_dof(), 0. );
        scalarVector_Type fracturePressureMeanInCutInterpolated ( mfproj.nb_dof(), 0. );
        scalarVector_Type fracturePressureMeanOutCutInterpolated ( mfproj.nb_dof(), 0. );

        getfem::mesh_fem mfprojUncut ( fracture->getMeshMapped(),
                                       fracture->getMeshFEMPressure().get_qdim());
        mfprojUncut.set_finite_element ( fractureFETypePressure );

        getfem::interpolation ( mfprojUncut, mfproj,
                                fracturePressureMeanUNCUT, fracturePressureMeanUNCUTInterpolated );

        const sizeVector_Type& extendedPressure = fracture->getExtendedPressure();

        for ( size_type otherFracture = 0; otherFracture < numberFractures; ++otherFracture )
        {
            if ( fracture->getLevelSetIntersect( otherFracture ).get() )
            {
                getfem::mesh_region regionMesh = fracture->getMeshFlat().region ( FractureHandler::FRACTURE_INTERSECT * ( f + 1 ) + otherFracture + 1 );
                size_type i_cv = 0;
                dal::bit_vector bc_cv = regionMesh.index();
                gmm::clear ( fracturePressureInCut );
                gmm::clear ( fracturePressureOutCut );

                for ( i_cv << bc_cv; i_cv != size_type(-1); i_cv << bc_cv )
                {
                    const size_type ibase = fracture->getMeshFEMPressure().ind_basic_dof_of_element ( i_cv )[0];
                    const size_type position = size_type ( std::find ( extendedPressure.begin(), extendedPressure.end(), ibase )
                                                - extendedPressure.begin());

                    const base_node point = mfprojUncut.point_of_basic_dof ( ibase );
                    const scalar_type levelSetValue = M_fractures->getFracture( otherFracture )->getLevelSet()->getData()->levelSetFunction( point );

                    if ( levelSetValue < 0 )
                    {
                        fracturePressureInCut [ ibase ] = (*(M_fracturePressure [ f ])) [ ibase ];
                        fracturePressureOutCut [ ibase ] = (*(M_fracturePressure [ f ])) [ fractureNumberDOFPressure [ f ] + position ];
                    }
                    else
                    {
                        fracturePressureOutCut [ ibase ] = (*(M_fracturePressure [ f ])) [ ibase ];
                        fracturePressureInCut [ ibase ] = (*(M_fracturePressure [ f ])) [ fractureNumberDOFPressure [ f ] + position ];
                    }

                }

                gmm::clear ( fracturePressureMeanInCutInterpolated );

                getfem::interpolation ( mfprojUncut, mfproj,
                                        fracturePressureInCut, fracturePressureMeanInCutInterpolated );

                gmm::clear ( fracturePressureMeanOutCutInterpolated );

                getfem::interpolation ( mfprojUncut, mfproj,
                                        fracturePressureOutCut, fracturePressureMeanOutCutInterpolated );

                i_cv = 0;
                bc_cv = meshFcutMapped.convex_index();
                for ( i_cv << bc_cv; i_cv != size_type(-1); i_cv << bc_cv )
                {
                    getfem::mesh_fem::ind_dof_ct idofs = mfproj.ind_basic_dof_of_element ( i_cv );
                    for ( size_type k = 0; k < idofs.size(); ++k )
                    {
                        const base_node node = mfproj.point_of_basic_dof(idofs [ k ]);
                        const scalar_type levelSetValue = M_fractures->getFracture( otherFracture )->getLevelSet()->getData()->levelSetFunction( node );
                        if ( levelSetValue < 0 )
                        {
                            fracturePressureMeanUNCUTInterpolated [ idofs[k] ] += fracturePressureMeanInCutInterpolated [ idofs[k] ];
                        }
                        else
                        {
                            fracturePressureMeanUNCUTInterpolated [ idofs[k] ] += fracturePressureMeanOutCutInterpolated [ idofs[k] ];
                        }

                    }

                }

            }
        }

        osFileName.str("");
        osFileName << "fracturePressure" << f << ".vtk";

        exportSolution ( M_exporter->getFolder() + osFileName.str(), "Pressure",
                     mfproj, fracturePressureMeanUNCUTInterpolated );
    }

   
   
}
