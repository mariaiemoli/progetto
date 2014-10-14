/** DarcyFractured.cc
 *
 * structure for the Darcy fractured problem
 *
 */

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
{} // costruttore


void DarcyFractured::init ( )
{

    // Allocate data
    const size_type shifCoefficients = M_mesh->getMeshFEMCoefficients().nb_dof();

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
        
        M_mediumEtaInterpolated[0] [ i ] = M_mediumData->invKDistribution11(node) * M_mediumData->getInvK();

        M_mediumEtaInterpolated[1] [ i ] = M_mediumData->invKDistribution12(node) * M_mediumData->getInvK();

        M_mediumEtaInterpolated[2] [ i ] = M_mediumData->invKDistribution22(node) * M_mediumData->getInvK();

        for ( size_type f = 0; f < numberFractures; ++f )
        {
            M_fractureEtaNormalOnMedium [ f ] [ i ] = M_fractures->getFracture( f )->getData().etaNormalDistribution(node)
                    								  * M_fractures->getFracture ( f )->getData().getEtaNormal();
        }
    }

} // init



void DarcyFractured::assembly ( )
{
	// numero totale delle fratture
    const scalar_type numberFractures = M_fractures->getNumberFractures();
	
    // vettore dei gradi di libertà per la pressione per le fratture
    sizeVector_Type fractureNumberDOFPressure(numberFractures);

    // numero di gradi di libertà totali per la pressione (con quelli estesi) per ogni frattura 
    sizeVector_Type fractureNumberGlobalDOFPressure(numberFractures);

    // vettore dei gradi di libertà per la velocità per le fratture
    sizeVector_Type fractureNumberDOFVelocity(numberFractures);

    // numero di gradi di libertà totali per la velocità (con quelli estesi) per ogni frattura 
    sizeVector_Type fractureNumberGlobalDOFVelocity(numberFractures);

    // numero totale dei gradi di libertà ( velocità + pressione ) per ogni frattura
    sizeVector_Type fractureNumberDOFVelocityPressure(numberFractures);

    // numero di gradi di libertà vincolati alla condizione al bordo per ogni frattura
    sizeVector_Type fractureNumberBoundaryDOF(numberFractures);

    // numero complessivo dei gradi di libertà ( pressione + velocità ) 
    size_type fractureTotalNumberDOFVelocityPressure = 0;

    // numero complessivo dei gradi di libertà di pressione 
    size_type fractureTotalNumberDOFPressure = 0;

    // numero totale di intersezioni
    size_type fractureNumberCross = 0;
    size_type fractureNumberBifurcation = 0;
    size_type fractureNumberIntersections = 0;
    size_type globalFractureNumber =0;
	
    for ( size_type f = 0; f < numberFractures; ++f )
    {
    	// numero di gradi di libertà per la pressione senza quelli estesi per la frattura f
        fractureNumberDOFPressure [ f ] = M_fractures->getFracture( f )->getMeshFEMPressure().nb_dof();

    	// numero di gradi di libertà per la pressione con quelli estesi per la frattura f
        fractureNumberGlobalDOFPressure [ f ] = fractureNumberDOFPressure [ f ] + M_fractures->getFracture( f )->getNumExtendedPressure();

    	// numero totale di gradi di libertà per la pressione
        fractureTotalNumberDOFPressure += fractureNumberGlobalDOFPressure [ f ];

    	// numero di gradi di libertà per la velocità senza quelli estesi per la frattura f
        fractureNumberDOFVelocity [ f ] = M_fractures->getFracture( f )->getMeshFEMVelocity().nb_dof();

    	// numero di gradi di libertà per la velocità con quelli estesi per la frattura f
        fractureNumberGlobalDOFVelocity [ f ] = fractureNumberDOFVelocity [ f ] + M_fractures->getFracture( f )->getNumExtendedVelocity();

    	// numero totale di gradi di libertà per la velocità e la pressione per la frattura f
        fractureNumberDOFVelocityPressure [ f ] = fractureNumberGlobalDOFVelocity [ f ] + fractureNumberGlobalDOFPressure [ f ];

    	// numero totale di gradi di libertà per la pressione e la velocità
        fractureTotalNumberDOFVelocityPressure += fractureNumberDOFVelocityPressure [ f ];

    	// numero di gradi di libertà per i bordi per la frattura f
        fractureNumberBoundaryDOF [ f ] = M_bcHandler->getFractureBC(f)->getMeshFEM().nb_dof();

    }
	
    // numero intersezioni
    fractureNumberCross = M_fractures->getIntersections ()->getNumberCross ();
   
    fractureNumberBifurcation = M_fractures->getIntersections ()->getNumberBifurcation ();
	
    fractureNumberIntersections = fractureNumberCross + 5*fractureNumberBifurcation;
    
    globalFractureNumber = fractureNumberCross*2 + fractureNumberBifurcation*6;
        
    // inizializziamo tutte le matrici a blocchi, la matrice globale e il termine di destra per il sistema e il vettore delle soluzioni 
    
    // Allochiamo la matrice globale del sistema:  M_darcyGlobalMatrix
    M_globalMatrix.reset( new sparseMatrix_Type( fractureTotalNumberDOFVelocityPressure + fractureNumberIntersections,
             	 	 	 	 	 	 	 	 	 fractureTotalNumberDOFVelocityPressure + globalFractureNumber ));
    gmm::clear( *M_globalMatrix );
    
    // Allochiamo il vettore del termine noto di destra del sistema: M_darcyGlobalRightHandSide
    M_globalRightHandSide.reset( new scalarVector_Type( fractureTotalNumberDOFVelocityPressure + globalFractureNumber ));
    gmm::clear( *M_globalRightHandSide );

    
    // Allochiamo il vettore globale del termine incognito del sistema: M_darcyVelocityAndPressure
    M_velocityAndPressure.reset(new scalarVector_Type( fractureTotalNumberDOFVelocityPressure + globalFractureNumber ));
    gmm::clear(*M_velocityAndPressure);

    
    // Blocks of matrix
    sparseMatrixPtr_Type  App;

    App.reset(new sparseMatrix_Type( fractureNumberIntersections , globalFractureNumber ));
    gmm::clear(*App);

    
    // Block of fracture matrix
    sparseMatrixPtrContainer_Type A11F(numberFractures), A12F(numberFractures), E(numberFractures);

    for ( size_type f = 0; f < numberFractures; ++f )
    {
        // Allocate the fracture matrix A11F
        A11F [ f ].reset(new sparseMatrix_Type ( fractureNumberGlobalDOFVelocity [ f ], fractureNumberGlobalDOFVelocity [ f ]) );
        gmm::clear(*(A11F [ f ]));

        // Allocate the fracture matrix A12F
        A12F [ f ].reset(new sparseMatrix_Type ( fractureNumberGlobalDOFVelocity [ f ], fractureNumberGlobalDOFPressure [ f ]) );
        gmm::clear(*(A12F [ f ]));
    }

 
    
    // Accoppiamo le fratture
    getfem::coupleFractures ( App, M_fractures );
    size_type shiftIntersect = 0;
    
    for ( size_type f = 0; f < numberFractures; ++f )
    {
        std::cout <<std::endl<< "Fracture " << f << std::endl;

        
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
            	// cosa fare?!
            	
            	std::cout << "gestire intersezione" << std::endl;
            }

        }
        
        shiftIntersect += fractureNumberDOFVelocityPressure [ f ];
      }
    
    // Shift for the fracture
    size_type fractureShift = 0;
    for ( size_type f = 0; f < numberFractures; ++f )
    {

        // Copy the matrix A11F in M_darcyGlobalMatrix in the correct position
        gmm::copy(*(A11F [ f ]), gmm::sub_matrix(*M_globalMatrix, gmm::sub_interval( fractureShift, fractureNumberGlobalDOFVelocity [ f ]),
                        										  gmm::sub_interval(fractureShift, fractureNumberGlobalDOFVelocity [ f ])));

        // Copy the matrix A12F in M_darcyGlobalMatrix in the correct position
        gmm::copy(*(A12F [ f ]), gmm::sub_matrix(*M_globalMatrix, gmm::sub_interval( fractureShift, fractureNumberGlobalDOFVelocity [ f ]),
        														  gmm::sub_interval(fractureShift + fractureNumberGlobalDOFVelocity [ f ], 
        																  	        fractureNumberGlobalDOFPressure [ f ])));

        // Copy the matrix -A12F in M_darcyGlobalMatrix in the correct position
        gmm::copy(gmm::transposed(gmm::scaled(*(A12F [ f ]), -1.0)), 
        						  gmm::sub_matrix( *M_globalMatrix, 
        								  gmm::sub_interval( fractureShift + fractureNumberGlobalDOFVelocity [ f ], fractureNumberGlobalDOFPressure [ f ]), 
        								  gmm::sub_interval( fractureShift, fractureNumberGlobalDOFVelocity [ f ])));
   
        // Update the shift
        fractureShift += fractureNumberDOFVelocityPressure [ f ];
 
    }

    gmm::copy(*App, gmm::sub_matrix(*M_globalMatrix, 
    		gmm::sub_interval( fractureTotalNumberDOFVelocityPressure, fractureNumberIntersections ), 
    		gmm::sub_interval( fractureTotalNumberDOFVelocityPressure, globalFractureNumber ) ));
    
    //Costruiamo il termine noto
    
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
        BstressF [ f ].reset(new scalarVector_Type( fractureNumberGlobalDOFVelocity [ f ], 0.));

        PneumannF [ f ].reset(new scalarVector_Type( fractureNumberBoundaryDOF [ f ], 0.));

        VdirichletF [ f ].reset(new scalarVector_Type( fractureNumberDOFPressure [ f ], 0.));

        B_vF [ f ].reset(new scalarVector_Type(fractureNumberGlobalDOFVelocity [ f ], 0.));

        B_pF [ f ].reset(new scalarVector_Type(fractureNumberGlobalDOFPressure [ f ], 0.));
        
    }

    for ( size_type f = 0; f < numberFractures; ++f )
    {
        const BCPtr_Type& fractureBC = M_bcHandler->getFractureBC(f);

        for ( size_type i = 0; i < fractureNumberBoundaryDOF [ f ]; i++ )
        {
            const base_node node = fractureBC->getMeshFEM().point_of_basic_dof(
                    i);

            (*(PneumannF [ f ])) [ i ] = M_fractures->getFracture( f )->getData().pressureExact(node);
        }

        (*(PneumannF [ f ])) [ 0 ] *= -1; //perché la normale ha il segno meno all'inizio di una roba 1D
    }

    for ( size_type f = 0; f < numberFractures; ++f )
    {
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
            (*M_globalRightHandSide) [ fractureShift + i ] += (*(BstressF [ f ])) [ i ];

            (*M_globalRightHandSide) [ fractureShift + i ] += (*(B_vF [ f ])) [ i ];
        }

        // Update the shift
        fractureShift += fractureNumberDOFVelocityPressure [ f ];
    }


    scalarVectorPtrContainer_Type divF(numberFractures);

    fractureShift = 0;
    
    for ( size_type f = 0; f < numberFractures; ++f )
    {
        divF [ f ].reset(new scalarVector_Type(fractureNumberDOFPressure [ f ], 0.));

        for ( size_type i = 0; i < fractureNumberDOFPressure [ f ]; i++ )
        {
            const base_node node = M_fractures->getFracture( f )->getMeshFEMPressure().point_of_basic_dof( i );

            (*(divF [ f ])) [ i ] = M_fractures->getFracture( f )->getData().darcySource( node );
        }

        gmm::clear(*(B_vF [ f ]));
        gmm::clear(*(B_pF [ f ]));

        getfem::assembling_Source_BoundaryF ( B_pF [ f ], divF [ f ], M_fractures->getFracture( f ), FractureHandler::FRACTURE_UNCUT * ( f + 1 ) );

        for ( size_type otherFracture = 0; otherFracture < numberFractures; ++otherFracture )
        {
            if ( M_fractures->getFracture( f )->getMeshLevelSetIntersect ( otherFracture ).get() )
            {
               /* getfem::assembling_SourceF ( B_pF [ f ], divF [ f ], M_fractures->getFracture( f ),
                                             M_fractures->getFracture ( otherFracture ),
                                             FractureHandler::FRACTURE_INTERSECT * ( f + 1 ) + otherFracture + 1 );
                                             */
            	std::cout << "sistemare intersezione" << std::endl;
            }
        }

    /*    for ( size_type i = 0; i < fractureNumberGlobalDOFVelocity [ f ]; ++i )
        {
            (*M_globalRightHandSide) [ fractureShift + i ] += (*(B_vF [ f ])) [ i ];
        }

        for ( size_type i = fractureNumberGlobalDOFVelocity [ f ]; i < fractureNumberDOFVelocityPressure [ f ]; ++i )
        {
                (*M_globalRightHandSide) [ fractureShift + i ] += (*(B_pF [ f ])) [ i - fractureNumberGlobalDOFVelocity [ f ] ];
        }
	   */
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

    
}// assembly