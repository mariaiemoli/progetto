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

    const size_type mediumNumberDOFVelocity =
            M_mesh->getMeshFEMVector().nb_dof();

    const size_type mediumNumberDOFPressure =
            M_mesh->getMeshFEMScalar().nb_dof();

    const size_type mediumNumberDOFCoefficients =
            M_mesh->getMeshFEMCoefficients().nb_dof();

    const size_type mediumNumberBoundaryDOF =
            M_bcHandler->getMediumBC()->getMeshFEM().nb_dof();

    const size_type intersectDOFPressure = M_mesh->getCountExtendedIntersectDOFScalar();
    const size_type intersectDOFVelocity = M_mesh->getCountExtendedIntersectDOFVector();

    // Degrees of freedom of the bulk primal: total degree plus extended elements
    const size_type mediumNumberDOFPressureGlobal = mediumNumberDOFPressure
            + M_mesh->getCountExtendedDOFScalar(numberFractures - 1) + intersectDOFPressure;

    // Degrees of freedom of the bulk dual: total degree plus extended elements
    const size_type mediumNumberDOFVelocityGlobal = mediumNumberDOFVelocity
            + M_mesh->getCountExtendedDOFVector(numberFractures - 1) + intersectDOFVelocity;

    // Degrees of freedom of the bulk: primal plus dual degrees for freedom
    const size_type mediumNumberDOFVelocityPressureGlobal =
            mediumNumberDOFVelocityGlobal + mediumNumberDOFPressureGlobal;

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

    // Initialize the dofs of intersect regions
    FractureIntersect::mapIntersection_Type& mapIntersect = M_fractures->getIntersections()->
                                                            getIntersections();

    FractureIntersect::mapIntersection_Type::iterator begin = mapIntersect.begin();
    FractureIntersect::mapIntersection_Type::iterator end = mapIntersect.end();
    FractureIntersect::mapIntersection_Type::iterator it;
    size_type mapIDIntersectType = 0;
    for ( it = begin; it != end; ++it, ++mapIDIntersectType )
    {
        IntersectDataContainer_Type& intersection = it->second;
        const FractureIntersect::IntersectionType intersectionType = it->first;
        const size_type num = intersection.size();

        for ( size_type i = 0; i < num; ++i )
        {

            sizeVector_Type dofPressure ( intersection[i].getRegionActive().size(), 0);
            sizeVectorContainer_Type dofVelocity ( 3 );
            for ( size_type i = 0; i < 3; ++i )
            {
                dofVelocity [i].resize( dofPressure.size() );
            }


            const size_type elementID = intersection[i].getElementID();
            const size_type numFractures = intersection[i].getNumFractures();

            // evaluation of level set in the dofs
            scalarVectorContainer_Type levelSetSignsVelocity ( 3 );

            scalarVector_Type levelSetSignsPressure;

            for ( size_type f = 0; f < numFractures; ++f )
            {
                size_type ind_dof = M_mesh->getMeshFEMScalar().ind_basic_dof_of_element ( elementID )[0];
                const LevelSetHandlerPtr_Type& levelSet = intersection[i].getFracture(f)->getLevelSet();
                levelSetSignsPressure.push_back ( levelSet->getBaricenterValue( ind_dof ) );
            }

            for ( size_type k = 0; k < 3; ++k )
            {
                for ( size_type f = 0; f < numFractures; ++f )
                {
                    size_type ind_dof = M_mesh->getMeshFEMVector().ind_basic_dof_of_element ( elementID )[k];
                    const LevelSetHandlerPtr_Type& levelSet = intersection[i].getFracture(f)->getLevelSet();
                    levelSetSignsVelocity[k].push_back ( levelSet->getDOFValue( ind_dof ) );
                }
            }

            for ( size_type k = 0; k < intersection[i].getRegionActive().size(); ++k )
            {
                const size_type dofPressureBase = M_mesh->getMeshFEMScalar().ind_basic_dof_of_element ( elementID )[0];

                const std::pair < std::string, size_type > comparaPressure = comparaSegni ( intersection[i].getRegionActive()[k], levelSetSignsPressure );

                if ( comparaPressure.first == "Base" )
                {
                    dofPressure [ k ] = dofPressureBase;
                }
                else if ( comparaPressure.first == "Extended" )
                {
                    const scalar_type fractureID = intersection[i].getFracture(comparaPressure.second)->getId();

                    const sizeVector_Type& extendedDOFScalar = M_mesh->getExtendedDOFScalar ( fractureID );

                    const size_type position = std::find ( extendedDOFScalar.begin(), extendedDOFScalar.end(), dofPressureBase ) -
                                                           extendedDOFScalar.begin();

                    dofPressure [ k ] = position + mediumNumberDOFPressure +
                                        M_mesh->getCountExtendedDOFScalar( fractureID - 1. );


                }
                if ( comparaPressure.first == "Extra" )
                {
                    const sizeVector_Type& extendedInterDOFScalar = M_mesh->getExtendedIntersectDOFScalar();

                    const size_type position = std::find ( extendedInterDOFScalar.begin(), extendedInterDOFScalar.end(), dofPressureBase ) -
                                                           extendedInterDOFScalar.begin();

                    dofPressure [ k ] = position + mediumNumberDOFPressureGlobal - intersectDOFPressure;
                }

                for ( size_type j = 0; j < 3; ++j )
                {
                    const size_type dofVelocityBase = M_mesh->getMeshFEMVector().ind_basic_dof_of_element ( elementID )[j];

                    const std::pair < std::string, size_type > compara = comparaSegni ( intersection[i].getRegionActive()[k], levelSetSignsVelocity[j] );

                    if ( compara.first == "Base" )
                    {
                        dofVelocity [ j ][ k ] = dofVelocityBase;
                    }
                    else if ( compara.first == "Extended" )
                    {
                        const scalar_type fractureID = intersection[i].getFracture(compara.second)->getId();

                        const sizeVector_Type& extendedDOFVector = M_mesh->getExtendedDOFVector ( fractureID );
                        const size_type positionVelocity = std::find ( extendedDOFVector.begin(), extendedDOFVector.end(),
                                                           dofVelocityBase ) - extendedDOFVector.begin();

                        dofVelocity [ j ][ k ] = positionVelocity + mediumNumberDOFVelocity + M_mesh->getCountExtendedDOFVector( fractureID - 1. );

                    }
                    if ( compara.first == "Extra" )
                    {
                        const sizeVector_Type& extendedInterDOFVector = M_mesh->getExtendedIntersectDOFVector();

                        const size_type positionVelocity = std::find ( extendedInterDOFVector.begin(), extendedInterDOFVector.end(),
                                                                       dofVelocityBase ) - extendedInterDOFVector.begin();

                        dofVelocity [ j ] [ k ] = positionVelocity + mediumNumberDOFVelocityGlobal - intersectDOFVelocity;

                    }

                }

            }

        intersection[i].setDOFPosition ( dofPressure, dofVelocity );
        }

    }

    // Initialize all the block matrices, global matrices and right hand side for the system and the solution vector

    // Allocate the global matrix M_darcyGlobalMatrix
    M_globalMatrix.reset(new sparseMatrix_Type(
            mediumNumberDOFVelocityPressureGlobal + fractureTotalNumberDOFVelocityPressure + fractureNumberIntersections,
            mediumNumberDOFVelocityPressureGlobal + fractureTotalNumberDOFVelocityPressure + fractureNumberIntersections ));
    gmm::clear(*M_globalMatrix);

    // Allocate the global vector M_darcyGlobalRightHandSide
    M_globalRightHandSide.reset(new scalarVector_Type(
            mediumNumberDOFVelocityPressureGlobal + fractureTotalNumberDOFVelocityPressure + fractureNumberIntersections ));
    gmm::clear(*M_globalRightHandSide);

    // Allocate the global vector M_darcyVelocityAndPressure
    M_velocityAndPressure.reset(new scalarVector_Type(
            mediumNumberDOFVelocityPressureGlobal + fractureTotalNumberDOFVelocityPressure + fractureNumberIntersections ));
    gmm::clear(*M_velocityAndPressure);

    // Allocate the vector M_darcyMediumVelocity
    M_mediumVelocity.reset(new scalarVector_Type(mediumNumberDOFVelocityGlobal));
    gmm::clear(*M_mediumVelocity);

    // Allocate the vector M_darcyMediumPressure
    M_mediumPressure.reset(new scalarVector_Type(mediumNumberDOFPressureGlobal));
    gmm::clear(*M_mediumPressure);

    // Blocks of matrix
    sparseMatrixPtr_Type A11, A12, A22, F, App;

    // Allocate the matrix A11
    A11.reset(new sparseMatrix_Type(mediumNumberDOFVelocityGlobal,
            mediumNumberDOFVelocityGlobal));
    gmm::clear(*A11);

    // Allocate the matrix A12
    A12.reset(new sparseMatrix_Type(mediumNumberDOFVelocityGlobal,
            mediumNumberDOFPressureGlobal));
    gmm::clear(*A12);

    // Allocate the matrix A22
    A22.reset(new sparseMatrix_Type(mediumNumberDOFPressureGlobal,
            mediumNumberDOFPressureGlobal));
    gmm::clear(*A22);

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

        // Allocate the matrix E
        // DA SISTEMARE!!!
        E [ f ].reset(new sparseMatrix_Type(mediumNumberDOFVelocityGlobal,
                mediumNumberDOFPressure));
        gmm::clear(*(E [ f ]));
    }

    // Assembling phase

    // Compute the edge measure of the fracture

    scalarVectorContainer_Type fractureEdgeMeasure(numberFractures);
    sparseMatrixPtrContainer_Type fractureMediumInterpolationMatrix(
            numberFractures);

    for ( size_type f = 0; f < numberFractures; ++f )
    {
        gmm::resize(fractureEdgeMeasure [ f ], fractureNumberDOFPressure [ f ]);

        getfem::stimaLati(fractureEdgeMeasure [ f ], M_fractures->getFracture( f ));

        // Compute the matrix fractureMediumInterpolationMatrix for interpolation between the fracture and the bulk
        fractureMediumInterpolationMatrix [ f ].reset(new sparseMatrix_Type(
                mediumNumberDOFCoefficients,
                fractureNumberGlobalDOFPressure[f] ) );

        gmm::clear(*fractureMediumInterpolationMatrix [ f ]);

        getfem::interpolationMatrix ( fractureMediumInterpolationMatrix [ f ],
                M_mesh, M_fractures->getFracture ( f ), f + FractureData::FRACTURE );

    }
    for ( it = mapIntersect.begin(); it != mapIntersect.end(); ++it, ++mapIDIntersectType )
    {
        IntersectDataContainer_Type& intersection = it->second;
        const FractureIntersect::IntersectionType intersectionType = it->first;
        const size_type num = intersection.size();

        for ( size_type i = 0; i < num; ++i )
        {
            FracturePtrContainer_Type fractures;
            sparseMatrixPtrContainer_Type interpolationMatrices;

            for ( size_type f = 0; f < intersection[i].getNumFractures(); ++f )
            {
                fractures.push_back( intersection[i].getFracture(f) );
                interpolationMatrices.push_back ( fractureMediumInterpolationMatrix [ fractures[f]->getId() ] );
            }

            getfem::interpolationMatrix ( interpolationMatrices, M_mesh, fractures, intersection[i] );

        }

    }

    std::cout << "Assembling medium" << std::endl;

    // Computes the matrix \int_K \tau_i \cdot \tau_j
    getfem::massHdiv ( A11, M_mesh, M_mediumEtaInterpolated,
                       MeshHandler::UNCUT_REGION );

    getfem::nitsche ( A11, M_mesh, M_mediumEtaInterpolated,
                      M_mediumData->getPenaltyVector(), M_bcHandler->getDirichletUncut(),
                      MeshHandler::UNCUT_REGION );

    // Computes the matrix \int_K \nabla v_i \cdot \tau_j
    getfem::divHdiv ( A12, M_mesh, MeshHandler::UNCUT_REGION );

    std::cout << "App dim " << fractureNumberIntersections/2 << " e " << fractureNumberIntersections << std::endl;

    getfem::coupleFractures ( App, M_fractures );
    size_type shiftIntersect = mediumNumberDOFVelocityPressureGlobal;

    for ( size_type f = 0; f < numberFractures; ++f )
    {
        std::cout << "Fracture " << f << std::endl;

        // Computes the matrix \int_K \tau_i \cdot \tau_j +  \int_K {\tau_i \cdot n}{\tau_j \cdot n} + \int_K [\tau_i \cdot n][\tau_j \cdot n]
        getfem::massHdiv ( A11, M_mesh, M_fractures->getFracture ( f ),
                           M_fractureEtaNormalOnMedium [ f ], M_mediumEtaInterpolated, f
                           + FractureData::FRACTURE  );

        getfem::darcy_E ( E [ f ], M_mesh, M_fractures->getFracture( f ), f
                          + FractureData::FRACTURE );

        getfem::nitsche ( A11, M_mesh, M_fractures->getFracture( f ), M_mediumEtaInterpolated,
                          M_mediumData->getPenaltyVector(), M_bcHandler->getDirichletCut(
                          FractureData::FRACTURE + f) );

        getfem::divHdiv ( A12, M_mesh, M_fractures->getFracture( f ), f
                          + FractureData::FRACTURE );

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
                        gmm::sub_interval ( mediumNumberDOFVelocityPressureGlobal + fractureTotalNumberDOFVelocityPressure + globalIndex, 1 ) ) );

                    gmm::copy ( gmm::transposed(*Aup), gmm::sub_matrix (*M_globalMatrix,
                        gmm::sub_interval ( mediumNumberDOFVelocityPressureGlobal + fractureTotalNumberDOFVelocityPressure +
                        std::max(globalIndex, globalIndex2), 1 ),
                        gmm::sub_interval ( shiftIntersect, fractureNumberGlobalDOFVelocity [ f ] ) ) );

                }

            }
        }

        shiftIntersect += fractureNumberDOFVelocityPressure [ f ];
        //massLumping(*(A11F [ f ]));
//FINE MODIFICATO

    }

    mapIDIntersectType = 0;
    for ( it = mapIntersect.begin(); it != mapIntersect.end(); ++it, ++mapIDIntersectType )
    {
        IntersectDataContainer_Type& intersection = it->second;
        const FractureIntersect::IntersectionType intersectionType = it->first;
        const size_type num = intersection.size();

        for ( size_type i = 0; i < num; ++i )
        {
            std::vector < scalarVector_Type* > fracturesEta;
            for ( size_type f = 0; f < intersection[i].getNumFractures(); ++f )
            {
                const FractureHandlerPtr_Type& fracture = intersection[i].getFracture(f);
                const size_type fractureID = fracture->getId();
                getfem::darcy_E ( E[fractureID], fracture, intersection[i], M_mesh, f );
                fracturesEta.push_back ( &(M_fractureEtaNormalOnMedium [ fractureID ]) );
            }

            getfem::massHdiv ( A11, intersection[i], M_mesh, M_mediumEtaInterpolated, fracturesEta );
            getfem::divHdiv ( A12, intersection[i], M_mesh );

        }

    }

    // Copy blocks into the system matrix M_darcyGlobalMatrix
    //     [  A11  A12  0          ]
    // M = [ -A12  A22  0          ]
    //     [  0    0    A11F  A12F ]
    //     [  0    0   -A12F  0    ]


    // Copy the matrix A11 in M_darcyGlobalMatrix in the correct position
    gmm::copy(*A11, gmm::sub_matrix(*M_globalMatrix, gmm::sub_interval(0,
            mediumNumberDOFVelocityGlobal), gmm::sub_interval(0,
            mediumNumberDOFVelocityGlobal)));

    // Copy the matrix A12 in M_darcyGlobalMatrix in the correct position
    gmm::copy(*A12, gmm::sub_matrix(*M_globalMatrix, gmm::sub_interval(0,
            mediumNumberDOFVelocityGlobal), gmm::sub_interval(
            mediumNumberDOFVelocityGlobal, mediumNumberDOFPressureGlobal)));

    // Copy the matrix -A12 in M_darcyGlobalMatrix in the correct position
    gmm::copy(gmm::transposed(gmm::scaled(*A12, -1.0)), gmm::sub_matrix(
            *M_globalMatrix, gmm::sub_interval(mediumNumberDOFVelocityGlobal,
                    mediumNumberDOFPressureGlobal), gmm::sub_interval(0,
                    mediumNumberDOFVelocityGlobal)));

    // Shift for the fracture
    size_type fractureShift = mediumNumberDOFVelocityPressureGlobal;
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

        //--------------------------------qui E-----------------------------
        // E è il termine che accoppia la soluzione nella frattura a quella esterna.
        // uso E e -E', dovrebbe essere una garanzia per l'energia perché uso la
        //stessa interpolazione e approssimazione "nei due sensi"

        scalarVector_Type provv ( fractureNumberGlobalDOFPressure [ f ] ),
                          provv1( mediumNumberDOFPressure );

        for ( size_type i = 0; i < mediumNumberDOFVelocityGlobal; ++i )
        {
            // n columns of the block
            gmm::clear(provv);
            gmm::clear(provv1);

            for ( size_type j = 0; j < mediumNumberDOFPressure; ++j )
            {
                provv1 [ j ] = -(*(E [ f ]))(i, j);
            }

            gmm::mult_add(gmm::transposed(
                    *fractureMediumInterpolationMatrix [ f ]), provv1, provv);
            for ( size_type j = 0; j < fractureNumberGlobalDOFPressure [ f ]; ++j )
            {
                (*M_globalMatrix)(fractureShift + fractureNumberGlobalDOFVelocity [ f ] + j, i) = provv [ j ];
            }
        }

        for ( size_type i = 0; i < mediumNumberDOFVelocityGlobal; ++i )
        {
            // n rows of the block
            gmm::clear(provv);
            gmm::clear(provv1);
            for ( size_type j = 0; j < mediumNumberDOFPressure; ++j )
            {
                provv1 [ j ] = (*(E [ f ]))(i, j);
            }
            gmm::mult_add(gmm::transposed(
                    *fractureMediumInterpolationMatrix [ f ]), provv1, provv);

            for ( size_type j = 0; j < fractureNumberGlobalDOFPressure [ f ]; ++j )
            {
               (*M_globalMatrix)(i, fractureShift + fractureNumberGlobalDOFVelocity [ f ] + j) = provv [ j ];
            }
        }

        // Update the shift
        fractureShift += fractureNumberDOFVelocityPressure [ f ];

        const sizeVector_Type& extendedDOF = M_mesh->getExtendedDOFScalar(f);
        for ( size_type i = 0; i < extendedDOF.size(); ++i )
        {
            size_type j = extendedDOF [ i ];
            sizeVector_Type::iterator it;
            sizeVector_Type ele(1);
            sizeVector_Type& nonCut = M_mesh->getNonCut();
            ele [ 0 ] = M_mesh->getMeshFEMScalar().first_convex_of_basic_dof(j);
            it = std::search(nonCut.begin(), nonCut.end(), ele.begin(),
                    ele.end());
            if ( it != nonCut.end())
            {
                (*A22)(i + mediumNumberDOFPressure, i + mediumNumberDOFPressure)
                        = M_mesh->getMesh().convex_area_estimate(
                                M_mesh->getMeshFEMScalar().first_convex_of_basic_dof(
                                        j));
            }
        }

    }

    // Copy the matrix A22 in M_darcyGlobalMatrix in the correct position
    gmm::copy(*A22, gmm::sub_matrix(*M_globalMatrix, gmm::sub_interval(
            mediumNumberDOFVelocityGlobal, mediumNumberDOFPressureGlobal),
            gmm::sub_interval(mediumNumberDOFVelocityGlobal,
                    mediumNumberDOFPressureGlobal)));

    gmm::copy(*App, gmm::sub_matrix(*M_globalMatrix, gmm::sub_interval(
        mediumNumberDOFVelocityPressureGlobal + fractureTotalNumberDOFVelocityPressure,
        fractureNumberIntersections/2 ), gmm::sub_interval(
        mediumNumberDOFVelocityPressureGlobal + fractureTotalNumberDOFVelocityPressure,
        fractureNumberIntersections ) ));
    // Data for the boundary conditions
    // Vector, Neumann conditions
    scalarVectorPtr_Type Bstress;
    Bstress.reset(new scalarVector_Type(mediumNumberDOFVelocityGlobal, 0.));

    // Pressure term (Neumann)
    scalarVectorPtr_Type Pneumann;
    Pneumann.reset(new scalarVector_Type(mediumNumberBoundaryDOF, 0.));

    // Pressure term (Neumann)
    scalarVectorPtr_Type PneumannIn;
    PneumannIn.reset(new scalarVector_Type(mediumNumberBoundaryDOF, 0.));

    // Pressure term (Neumann)
    scalarVectorPtr_Type PneumannOut;
    PneumannOut.reset(new scalarVector_Type(mediumNumberBoundaryDOF, 0.));

    // Velocity term (Dirichilet)
    scalarVectorPtr_Type Vdirichlet;
    Vdirichlet.reset(new scalarVector_Type(mediumNumberDOFCoefficients, 0.));

    // Right Hand Side velocity
    scalarVectorPtr_Type B_v;
    B_v.reset(new scalarVector_Type(mediumNumberDOFVelocityGlobal, 0.));

    // Right Hand Side: velocity and pressure
    scalarVectorPtr_Type B_p;
    B_p.reset(new scalarVector_Type(mediumNumberDOFPressureGlobal, 0.));

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

    // Fill External Stress values for pressure
    const BCPtr_Type& mediumBC = M_bcHandler->getMediumBC();

    for ( size_type i = 0; i < mediumNumberBoundaryDOF; i++ )
    {
        const base_node node = mediumBC->getMeshFEM().point_of_basic_dof(i);

        (*Pneumann) [ i ] = M_mediumData->exact(node);

        (*PneumannIn) [ i ] = M_mediumData->exactInlet(node);

        (*PneumannOut) [ i ] = M_mediumData->exactOutlet(node);
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

    // Computes the right hand side rhs which stores the boundary conditions
    getfem::stressRHS ( Bstress, M_mesh, M_bcHandler, Pneumann );

    getfem::nitscheRHS ( B_v, M_mesh, M_bcHandler, M_mediumData->getInvK(),
                         M_mediumData->getPenaltyVector(), Vdirichlet );

    fractureShift = 0;
    for ( size_type f = 0; f < numberFractures; ++f )
    {

        getfem::stressRHS ( Bstress, M_mesh, M_bcHandler, M_fractures->getFracture( f ),
                            fractureShift, PneumannIn, PneumannOut );

        getfem::nitscheRHS ( B_v, M_mesh, M_bcHandler, M_fractures->getFracture( f ),
                             fractureShift, M_mediumData->getInvK(),
                             M_mediumData->getPenaltyVector(), Vdirichlet );

//        fractureShift += fractureNumberGlobalDOFVelocity [ f ];
        fractureShift += M_mesh->getExtendedDOFVector(f).size();

        // Computes the right hand side rhs which stores the boundary conditions for the fracture
        getfem::darcy_dataF ( BstressF [ f ], B_vF [ f ], M_bcHandler,
                              M_fractures->getFracture( f ), M_mediumData->getPenaltyVector(),
                              M_mediumData->getInvK(), PneumannF [ f ], VdirichletF [ f ] );
    }

    for ( size_type i = 0; i < mediumNumberDOFVelocityGlobal; ++i )
    {
        (*M_globalRightHandSide) [ i ] += (*Bstress) [ i ];
    }

    for ( size_type i = 0; i < mediumNumberDOFVelocityGlobal; ++i )
    {
        (*M_globalRightHandSide) [ i ] += (*B_v) [ i ];
    }

    fractureShift = mediumNumberDOFVelocityPressureGlobal;
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

    // Assemble the source term
    scalarVectorPtr_Type div;
    div.reset(new scalarVector_Type(mediumNumberDOFCoefficients, 0.));

    scalarVectorPtrContainer_Type divF(numberFractures);

    fractureShift = mediumNumberDOFVelocityPressureGlobal;
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

    for ( size_type i = 0; i < mediumNumberDOFCoefficients; i++ )
    {
        const base_node node =
                M_mesh->getMeshFEMCoefficients().point_of_basic_dof(i);
        (*div) [ i ] = M_mediumData->source(node);
    }

    gmm::clear(*B_p);

    //questo tiene conto del termine sorgenre
    getfem::sourceL2(B_p, div, M_mesh, MeshHandler::UNCUT_REGION);

    fractureShift = 0;
    for ( size_type f = 0; f < numberFractures; ++f )
    {
        getfem::sourceL2 ( B_p, div, M_mesh, M_fractures->getFracture( f ), f + FractureData::FRACTURE );
    }

    mapIDIntersectType = 0;
    for ( it = mapIntersect.begin(); it != mapIntersect.end(); ++it, ++mapIDIntersectType )
    {
        IntersectDataContainer_Type& intersection = it->second;
        const FractureIntersect::IntersectionType intersectionType = it->first;
        const size_type num = intersection.size();

        for ( size_type i = 0; i < num; ++i )
        {
                getfem::sourceL2 ( B_p, div, intersection[i], M_mesh );
        }

    }

    for ( size_type i = mediumNumberDOFVelocityGlobal; i
            < mediumNumberDOFVelocityPressureGlobal; ++i )
    {
        (*M_globalRightHandSide) [ i ] += (*B_p) [ i
                - mediumNumberDOFVelocityGlobal ];
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

//(*M_globalRightHandSide) [M_globalRightHandSide->size()-1 ] = 1.;
//(*M_globalRightHandSide) [M_globalRightHandSide->size()-2 ] = 1.;

    M_exporter->spy(M_globalMatrix, "./matrice.mm");
    M_exporter->spy(M_globalRightHandSide, "./rhs.mm");

} // assembly

// Solve the Darcy for the governing flux and do the time loop for the evolution problem
void DarcyFractured::solve ( )
{
    const scalar_type numberFractures = M_fractures->getNumberFractures();
    const scalar_type numberIntersections = M_fractures->getIntersections()
                                            ->getNumberIntersections();

    const size_type dofVelocity = M_mesh->getMeshFEMVector().nb_basic_dof_of_element(0);
    const size_type dofPressure = M_mesh->getMeshFEMScalar().nb_basic_dof_of_element(0);

    const size_type mediumNumberDOFVelocity =
            M_mesh->getMeshFEMVector().nb_dof();

    const size_type mediumNumberDOFPressure =
            M_mesh->getMeshFEMScalar().nb_dof();

    const size_type intersectDOFPressure = M_mesh->getCountExtendedIntersectDOFScalar();
    const size_type intersectDOFVelocity = M_mesh->getCountExtendedIntersectDOFVector();

    const size_type mediumNumberDOFVelocityGlobal = mediumNumberDOFVelocity
            + M_mesh->getCountExtendedDOFVector(numberFractures - 1) + intersectDOFVelocity;

    const size_type mediumNumberDOFPressureGlobal = mediumNumberDOFPressure
            + M_mesh->getCountExtendedDOFScalar(numberFractures - 1) + intersectDOFPressure;

    const size_type mediumNumberDOFVelocityPressureGlobal =
            mediumNumberDOFPressureGlobal + mediumNumberDOFVelocityGlobal;

    gmm::resize(*M_mediumPressure, mediumNumberDOFPressureGlobal);
    gmm::resize(*M_mediumVelocity, mediumNumberDOFVelocityGlobal);

    M_intersectionVelocity.reset ( new scalarVector_Type ( intersectDOFVelocity, 0 ) );
    M_intersectionPressure.reset ( new scalarVector_Type ( intersectDOFPressure, 0 ) );

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

    /*for ( size_type f = 0; f < numberFractures; ++f )
    {
        std::ostringstream osFileNameMesh;
        osFileNameMesh << "cmeshF2-" << f << ".vtk";
        exportMesh(M_exporter->getFolder() + osFileNameMesh.str(),
                M_fractures->getFracture ( f )->getMeshMapped());
    }
    */
//    scalarVector_Type Vx(mediumNumberDOFPressure, 0.), Vy(
//            mediumNumberDOFPressure, 0.);

    M_mediumVelocityInlet.reset(new scalarVector_Type(mediumNumberDOFVelocity,
            0.));
    M_mediumVelocityOutlet.reset(new scalarVector_Type(mediumNumberDOFVelocity,
            0.));

    // Solve the Darcy problem
    std::cout << std::endl << "Solving problem in Omega..." << std::flush;
    scalar_type roundConditionNumber;
    int ris;
    ris=SuperLU_solve(*M_globalMatrix, *M_velocityAndPressure,
                  *M_globalRightHandSide, roundConditionNumber);
    std::cout << " completed!" << std::endl;

    // Extract the dual in the bulk
    gmm::copy(gmm::sub_vector(*M_velocityAndPressure, gmm::sub_interval(0,
            mediumNumberDOFVelocityGlobal)), *M_mediumVelocity);

    // Extract the primal in the bulk
    gmm::copy(gmm::sub_vector(*M_velocityAndPressure, gmm::sub_interval(
            mediumNumberDOFVelocityGlobal, mediumNumberDOFPressureGlobal)),
            *M_mediumPressure);

    size_type fractureShift = mediumNumberDOFVelocityPressureGlobal;
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
    size_type intersectShift = mediumNumberDOFVelocityPressureGlobal
                               - intersectDOFVelocity - intersectDOFPressure;

    // extra dof for the dual
    gmm::copy ( gmm::sub_vector(*M_velocityAndPressure,
                gmm::sub_interval( intersectShift, intersectDOFVelocity ) ),
                *M_intersectionVelocity );

    // extra dof for the primal
    gmm::copy ( gmm::sub_vector(*M_velocityAndPressure,
                gmm::sub_interval( intersectShift + intersectDOFVelocity, intersectDOFPressure ) ),
                *M_intersectionPressure );

    getfem::pfem fractureFETypePressure;
    if ( numberFractures > 0 )
    {
        fractureFETypePressure = getfem::fem_descriptor(
                M_fractures->getFracture( 0 )->getData().getFEMTypePressure());
    }

    getfem::mesh_level_set mLevelSetCut(M_mesh->getMesh());
    for ( size_type f = 0; f < numberFractures; ++f )
    {
        mLevelSetCut.add_level_set(
                M_fractures->getFracture( f )->getLevelSet()->getLevelSet());
    }

    mLevelSetCut.adapt();

    // The following finite element space is useful for visualization
    getfem::mesh mcut;
    mLevelSetCut.global_cut_mesh(mcut);

    // Discontinuous finite element space
    getfem::mesh_fem mfcut(mcut, M_mesh->getMeshFEMScalar().get_qdim());
    mfcut.set_classical_discontinuous_finite_element(3, 0.01);

    exportMesh(M_exporter->getFolder() + "cmesh.vtk", mcut);

    scalarVector_Type darcyPressureInlet(mediumNumberDOFPressure, 0.0);
    scalarVector_Type darcyPressureOutlet(mediumNumberDOFPressure, 0.0);
    scalarVector_Type darcyPressureMean(mfcut.nb_dof(), 0.0);

    scalarVector_Type darcyVelocityMean(mediumNumberDOFVelocity, 0.0);
    scalarVector_Type darcyVelocityInlet1, darcyVelocityInlet2,
            darcyVelocityOutlet1, darcyVelocityOutlet2;

  
    //valore del level set
    for ( size_type f = 0; f < numberFractures; ++f )
    {
        std::ostringstream osFileName;
        osFileName << "MediumLevelSet" << f << ".vtk";
        exportSolution(M_exporter->getFolder() + osFileName.str(), "LevelSet",
                M_mesh->getMeshFEMScalar(),
                M_fractures->getFracture( f )->getLevelSet()->getBaricenterValue());
    }
   
    size_type shiftVelocity(mediumNumberDOFVelocityGlobal);
    fractureShift = 0;
    for ( size_type f = 0; f < numberFractures; ++f )
    {
        scalarVector_Type darcyPressureInletFracture(mfcut.nb_dof(), 0.);
        scalarVector_Type darcyPressureOutletFracture(mfcut.nb_dof(), 0.);

        gmm::clear(darcyPressureInlet);
        gmm::clear(darcyPressureOutlet);
        const sizeVector_Type& extendedDOF = M_mesh->getExtendedDOFScalar(f);
        const size_type shiftExtended = extendedDOF.size();

        for ( size_type i = 0; i < shiftExtended; ++i )
        {
            const size_type j = extendedDOF [ i ];
            const size_type el =
                    M_mesh->getMeshFEMScalar().first_convex_of_basic_dof(j);

            if ( M_mesh->getRegion(f + FractureData::FRACTURE).is_in(el) )
            {
                if ( M_fractures->getFracture( f )->getLevelSet()->getBaricenterValue(j) < 0 )
                {

                    darcyPressureInlet [ j ]
                            = (*M_velocityAndPressure) [ shiftVelocity + j ];

                    darcyPressureOutlet [ j ]
                            = (*M_velocityAndPressure) [ shiftVelocity
                                    + mediumNumberDOFPressure + i
                                    + fractureShift ];
                }
                else
                {
                    darcyPressureOutlet [ j ]
                            = (*M_velocityAndPressure) [ shiftVelocity + j ];

                    darcyPressureInlet [ j ]
                            = (*M_velocityAndPressure) [ shiftVelocity
                                    + mediumNumberDOFPressure + i
                                    + fractureShift ];
                }
            }

        }

        fractureShift += shiftExtended;

        getfem::interpolation(M_mesh->getMeshFEMScalar(), mfcut,
                darcyPressureInlet, darcyPressureInletFracture);

        getfem::interpolation(M_mesh->getMeshFEMScalar(), mfcut,
                darcyPressureOutlet, darcyPressureOutletFracture);

        for ( dal::bv_visitor cv(mcut.convex_index()); !cv.finished(); ++cv )
        {
            getfem::mesh_fem::ind_dof_ct idofs =
                    mfcut.ind_basic_dof_of_element(cv);
            scalar_type dmean = 0;

            for ( size_type i = 0; i < idofs.size(); ++i )
            {
                const base_node node = mfcut.point_of_basic_dof(idofs [ i ]);
                dmean
                        += M_fractures->getFracture( f )->getLevelSet()->getData()->levelSetFunction(
                                node);
            }

            if ( dmean < 0 )
            {
                for ( size_type i = 0; i < idofs.size(); ++i )
                {
                    size_type j = idofs [ i ];
                    darcyPressureOutletFracture [ j ]
                            = darcyPressureInletFracture [ j ];

                }
            }

        }
        for ( size_type i = 0; i < darcyPressureOutletFracture.size(); ++i )
        {
            darcyPressureMean [ i ] += darcyPressureOutletFracture [ i ];
        }

    }

    // Take the solution for the intersected dof
    FractureIntersect::mapIntersection_Type& mapIntersect = M_fractures->getIntersections()->
                                                            getIntersections();

    FractureIntersect::mapIntersection_Type::const_iterator begin = mapIntersect.begin();
    FractureIntersect::mapIntersection_Type::const_iterator end = mapIntersect.end();
    FractureIntersect::mapIntersection_Type::const_iterator it;
    size_type mapIDIntersectType = 0;

    scalarVectorContainer_Type darcyPressureIntersect;
    scalarVectorContainer_Type darcyPressureIntersectCut;

    for ( it = begin; it != end; ++it, ++mapIDIntersectType )
    {
        const IntersectDataContainer_Type& intersection = it->second;
        const FractureIntersect::IntersectionType intersectionType = it->first;
        const size_type num = intersection.size();

        const size_type basisFunction = M_fractures->getIntersections()->
                                        getBasisFunctionOfType ( intersectionType );

        const sizeVector_Type& extendDOFScalar =
                               M_mesh->getExtendedIntersectDOFScalar ();

        const sizeVector_Type& extendDOFVector =
                               M_mesh->getExtendedIntersectDOFVector ();

        for ( size_type i = 0; i < num; ++i )
        {
            const stringContainer_Type& subRegion = intersection[i].getRegionActive();
            const size_type numSubRegion = subRegion.size();
            darcyPressureIntersect.resize ( numSubRegion );
            darcyPressureIntersectCut.resize ( numSubRegion );
            for ( size_type u = 0; u < numSubRegion; ++u )
            {
                darcyPressureIntersect[u].resize( mediumNumberDOFPressure );
                darcyPressureIntersectCut[u].resize( mfcut.nb_dof() );

                gmm::clear ( darcyPressureIntersect[u] );
                gmm::clear ( darcyPressureIntersectCut[u] );
            }

            const size_type elementID = intersection[i].getElementID();
            const size_type numFractures = intersection[i].getNumFractures();
            const size_type basicDOF = M_mesh->getMeshFEMScalar().ind_basic_dof_of_element ( elementID )[0];

            for ( size_type u = 0; u < numSubRegion; ++u )
            {
                const size_type dof = intersection[i].getDOFPressure (u);
                darcyPressureIntersect [ u ][ basicDOF ] = (*M_velocityAndPressure) [ dof + mediumNumberDOFVelocityGlobal ];
            }

            // interpolo
            for ( size_type u = 0; u < numSubRegion; ++u )
            {
                getfem::interpolation ( M_mesh->getMeshFEMScalar(), mfcut,
                        darcyPressureIntersect[u], darcyPressureIntersectCut[u] );

            }

            for ( dal::bv_visitor cv(mcut.convex_index()); !cv.finished(); ++cv )
            {

                getfem::mesh_fem::ind_dof_ct idofs = mfcut.ind_basic_dof_of_element ( cv );
                scalarVector_Type levelSetValue ( numFractures );

                for ( size_type k = 0; k < idofs.size(); ++k )
                {
                    const base_node node = mfcut.point_of_basic_dof(idofs [ k ]);
                    for ( size_type f = 0; f < numFractures; ++f )
                    {
                        levelSetValue[f] = intersection[i].getFractures()[f]->
                                           getLevelSet()->getData()->levelSetFunction(node);
                    }
                    std::string regionSign = regionSigns( levelSetValue );
                    size_type position = intersection[i].getIndexRegion ( regionSign );
    
                    if ( isInTriangle ( M_mesh->getMesh(), elementID, node) )
                    {
                        darcyPressureIntersectCut[0][ idofs[k] ] = darcyPressureIntersectCut[position][ idofs[k] ];
                    }
                }

            }
            for ( size_type k = 0; k < darcyPressureIntersectCut[0].size(); ++k )
            {
                    darcyPressureMean[k] += darcyPressureIntersectCut[0][k];
            }

        }
    }

    scalarVector_Type darcyPressureMeanUNCUT(mediumNumberDOFPressure, 0.);

    for ( size_type i = 0; i < mediumNumberDOFPressure; ++i )
    {
        const size_type el =
                M_mesh->getMeshFEMScalar().first_convex_of_basic_dof(i);

        if ( M_mesh->getRegion(MeshHandler::UNCUT_REGION).is_in(el) )
        {
            darcyPressureMeanUNCUT [ i ]
                    = (*M_velocityAndPressure) [ shiftVelocity + i ];
        }
    }

    scalarVector_Type darcyPressureMeanUNCUTInterpolated(mfcut.nb_dof(), 0.);

    getfem::interpolation(M_mesh->getMeshFEMScalar(), mfcut,
            darcyPressureMeanUNCUT, darcyPressureMeanUNCUTInterpolated);

    for ( size_type i = 0; i < darcyPressureMean.size(); ++i )
    {
        darcyPressureMean [ i ] += darcyPressureMeanUNCUTInterpolated [ i ];
    }

  
    // Export discontinuous solution
    exportSolution ( M_exporter->getFolder() + "mediumPressure.vtk", "Pressure",
            mfcut, darcyPressureMean );

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

   
    for ( size_type f = 0; f < numberFractures; ++f )
    {
        // Split the velocitin in "in" and "out", useful for the transport problem
        for ( size_type i = 0; i < mediumNumberDOFVelocity; ++i )
        {
            if ( M_fractures->getFracture ( f )->getLevelSet()->getDOFValue(i) < 0 )
            {
                (*M_mediumVelocityInlet) [ i ] = (*M_velocityAndPressure) [ i ];
                (*M_mediumVelocityOutlet) [ i ] = (*M_velocityAndPressure) [ i ];
            }
            else
            {
                (*M_mediumVelocityOutlet) [ i ] = (*M_velocityAndPressure) [ i ];
                (*M_mediumVelocityInlet) [ i ] = (*M_velocityAndPressure) [ i ];
            }
        }

        // Do the same for the extended degrees of freedom
        const sizeVector_Type& extendedDOFVelocity =
                M_mesh->getExtendedDOFVector(f);
        const size_type shiftExtended = extendedDOFVelocity.size();
        for ( size_type i = 0; i < shiftExtended; ++i )
        {
            size_type j = extendedDOFVelocity [ i ];
            if ( M_fractures->getFracture( f )->getLevelSet()->getDOFValue(j) < 0 )
            {
                (*M_mediumVelocityInlet) [ j ] = (*M_velocityAndPressure) [ j ];
                (*M_mediumVelocityOutlet) [ j ]
                        = (*M_velocityAndPressure) [ mediumNumberDOFVelocity
                                + i ];
            }
            if ( M_fractures->getFracture( f )->getLevelSet()->getDOFValue(j) >= 0 )
            {
                (*M_mediumVelocityOutlet) [ j ] = (*M_velocityAndPressure) [ j ];
                (*M_mediumVelocityInlet) [ j ]
                        = (*M_velocityAndPressure) [ mediumNumberDOFVelocity
                                + i ];
            }

        }
    }

    // Compute the reconstruction of the Darcy velocity in each elements
    M_mediumVelocityInterpolatedAbscissa.reset(new scalarVector_Type(
            mediumNumberDOFPressure, 0));
    M_mediumVelocityInterpolatedOrdinate.reset(new scalarVector_Type(
            mediumNumberDOFPressure, 0));

    for ( size_type i = 0; i < mediumNumberDOFPressure; ++i )
    {
        const size_type elem =
                M_mesh->getMeshFEMScalar().ind_basic_dof_of_element(i) [ 0 ];

        const base_node point = M_mesh->getMeshFEMScalar().point_of_basic_dof(
                elem);

        scalarVector_Type V(2, 0.);
        getfem::interpolateVelocity(V, point, i, *M_mediumVelocity, M_mesh,
                M_globalMatrix);

        (*M_mediumVelocityInterpolatedAbscissa) [ elem ] = V [ 0 ];
        (*M_mediumVelocityInterpolatedOrdinate) [ elem ] = V [ 1 ];

    }

}
