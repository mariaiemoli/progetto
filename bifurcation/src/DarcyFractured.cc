
#include "../include/DarcyFractured.h"

/**************************************************************************/
/*  DarcyFractured.cc													  */
/*  Classe che assembla e risolve il sistema per il problema di Darcy	  */ 
/*  sulle fratture											              */
/**************************************************************************/


DarcyFractured::DarcyFractured ( const MediumDataPtr_Type& medium,
                                 const MeshHandlerPtr_Type& mesh,
                                 const BCHandlerPtr_Type& bcHandler,
                                 const FracturesSetPtr_Type& fractures,
                                 const ExporterPtr_Type& exporter ) :
            M_mediumData ( medium ), 
            M_mesh ( mesh ), 
            M_bcHandler ( bcHandler ), 
            M_fractures ( fractures ),
            M_exporter ( exporter ),
            M_fractureEtaNormalOnMedium ( M_fractures->getNumberFractures () ), 
            M_fractureVelocity ( M_fractures->getNumberFractures () ),
            M_fracturePressure ( M_fractures->getNumberFractures () )
{} // costruttore


void DarcyFractured::init ( )
{

    // Alloco i dati
    const size_type shifCoefficients = M_mesh->getMeshFEMCoefficients().nb_dof();

    const size_type numberFractures = M_fractures->getNumberFractures();

    for ( size_type f = 0; f < numberFractures; ++f )
    {
        // Alloco il vettore per M_fracturesEtaNormalOnMedium
        gmm::resize(M_fractureEtaNormalOnMedium [ f ], shifCoefficients);
        gmm::clear(M_fractureEtaNormalOnMedium [ f ]);

    }

    return;

} // init



void DarcyFractured::assembly ( const GetPot& dataFile )
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
	size_type fractureNumberBifurcation2 = 0;
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

    // Numero intersezioni
    fractureNumberCross = M_fractures->getIntersections ()->getNumberCross ();
           
    fractureNumberBifurcation = M_fractures->getIntersections ()->getNumberBifurcation ();
	
	fractureNumberBifurcation2 = M_fractures->getIntersections ()->getNumberBifurcation2 ();
    
    globalFractureNumber = fractureNumberCross*2 + fractureNumberBifurcation + 2*fractureNumberBifurcation2;
    
    // Inizializziamo tutte le matrici a blocchi, la matrice globale e il termine di destra per il sistema e il vettore delle soluzioni 
    
    // Allochiamo la matrice globale del sistema:  M_darcyGlobalMatrix
    M_globalMatrix.reset( new sparseMatrix_Type( fractureTotalNumberDOFVelocityPressure + globalFractureNumber,
             	 	 	 	 	 	 	 	 	 fractureTotalNumberDOFVelocityPressure + globalFractureNumber ));
												 
	
    gmm::clear( *M_globalMatrix );
 
    // Allochiamo il vettore del termine noto di destra del sistema: M_darcyGlobalRightHandSide
    M_globalRightHandSide.reset( new scalarVector_Type( fractureTotalNumberDOFVelocityPressure + globalFractureNumber ));
    gmm::clear( *M_globalRightHandSide );
    
    // Allochiamo il vettore globale del termine incognito del sistema: M_darcyVelocityAndPressure
    M_velocityAndPressure.reset(new scalarVector_Type( fractureTotalNumberDOFVelocityPressure + globalFractureNumber ));
    gmm::clear(*M_velocityAndPressure);
    
    // Allochiamo la matrice che usereme per accoppiare le fratture che si intersecano
    sparseMatrixPtr_Type  App;
	
    App.reset(new sparseMatrix_Type( fractureNumberCross , 2*fractureNumberCross ));
       gmm::clear(*App);

    // Matrici a blocchi per le fratture
    sparseMatrixPtrContainer_Type A11F(numberFractures), A12F(numberFractures);

    for ( size_type f = 0; f < numberFractures; ++f )
    {
        // Alloco la matrice A11F per la fratture f-esima
        A11F [ f ].reset(new sparseMatrix_Type ( fractureNumberGlobalDOFVelocity [ f ], fractureNumberGlobalDOFVelocity [ f ]) );
        gmm::clear(*(A11F [ f ]));

        // Alloco la matrice A12F per la fratture f-esima
        A12F [ f ].reset(new sparseMatrix_Type ( fractureNumberGlobalDOFVelocity [ f ], fractureNumberGlobalDOFPressure [ f ]) );
        gmm::clear(*(A12F [ f ]));
    }

    
    // Accoppio le fratture
    getfem::coupleFractures ( App, M_fractures, 1 );
    sizeVector_Type shiftIntersect ( numberFractures );
    
    shiftIntersect [ 0 ] = 0;
    
    // Calcolo i blocchi di matrici per ogni frattura
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
        
        if( f != 0)
        {
        	shiftIntersect [ f ] = shiftIntersect [ f-1 ] + fractureNumberDOFVelocityPressure [ f-1 ];	
        }
        
    }

    std::cout << std::endl;
    
    // Introduco il termine legato all'intersezione
    // A seconda del tipo di intersezione, Cross o Biforcazione, devo introdurre dei termini diversi
  
    IntersectDataContainer_Type IntCross = M_fractures->getIntersections ()-> getCrossIntersections ();
    
    IntersectDataContainer_Type IntBifurcation = M_fractures->getIntersections ()-> getBifurcationIntersections ();
        
    IntersectDataContainer_Type IntBifurcation2 = M_fractures->getIntersections ()-> getBifurcation2Intersections ();
	
    // Aggiorno prima le matrici di tutte le fratture che si intersecano formando un " Cross "
    for ( size_type i = 0; i < IntCross.size(); i++ )
    {
    	sparseMatrixPtr_Type Aup0, Aup1;
    	   	
    	FractureHandlerPtr_Type f0 = IntCross [ i ].getFracture (0);
    	FractureHandlerPtr_Type f1 = IntCross [ i ].getFracture (1);
    	
    	size_type id0 = f0->getId();
    	size_type id1 = f1->getId();
    	
        std::cout << "Cross: " << std::endl;
		std::cout << " Fractures " << id0 << ", " << id1 << std::endl;

    	
        Aup0.reset ( new sparseMatrix_Type ( fractureNumberGlobalDOFVelocity [ id0 ], 1 ) );
        gmm::clear(*Aup0);

        Aup1.reset ( new sparseMatrix_Type ( fractureNumberGlobalDOFVelocity [ id1 ], 1 ) );
        gmm::clear(*Aup1);
       
        // aggiorno per la frattura 0
        getfem::darcy_A11F_Cross ( A11F [ id0 ],f0,
								   f0->getEtaTangentialInterpolated(), 
								   f1,
								   FractureHandler::FRACTURE_INTERSECT * ( id0 + 1 ) + id1 + 1);

        getfem::darcy_A12F_Cross ( A12F [ id0 ], f0, f1,
                             	   FractureHandler::FRACTURE_INTERSECT * ( id0 + 1 ) + id1 + 1 );
        std::cout << std::endl;
        
        // aggiorno per la frattura 1
        getfem::darcy_A11F_Cross ( A11F [ id1 ],f1,
							 	   f1->getEtaTangentialInterpolated(),
							 	   f0,
							 	   FractureHandler::FRACTURE_INTERSECT * ( id1 + 1 ) + id0 + 1 );

		getfem::darcy_A12F_Cross ( A12F [ id1 ], f1, f0,
								   FractureHandler::FRACTURE_INTERSECT * ( id1 + 1 ) + id0 + 1 );

		std::cout << std::endl;
		
        const pairSizeVectorContainer_Type& intersectElementsGlobalIndex0 = f0->getFractureIntersectElementsGlobalIndex ();
        
        const size_type numIntersections = intersectElementsGlobalIndex0 [id1].size();
        
        const sizeVectorContainer_Type& intersectElements0 = f0->getFractureIntersectElements ();
        const sizeVectorContainer_Type& intersectElements1 = f1->getFractureIntersectElements ();
		
		for ( size_type k = 0; k < numIntersections; ++k )
		{
			gmm::clear (*Aup0);
			gmm::clear (*Aup1);

			getfem::velocityJump_Cross ( Aup0, f0, f1, intersectElements0 [ id1 ][ k ] );
			getfem::velocityJump_Cross ( Aup1, f1, f0, intersectElements1 [ id0 ][ k ] );

			const size_type globalIndex = intersectElementsGlobalIndex0 [ id1 ] [ k ].first;
			const size_type globalIndex2 = intersectElementsGlobalIndex0 [ id1 ] [ k ].second;
			
			gmm::copy ( *Aup0, 
					    gmm::sub_matrix (*M_globalMatrix,
									    gmm::sub_interval ( shiftIntersect [ id0 ], fractureNumberGlobalDOFVelocity [ id0 ] ),
									    gmm::sub_interval (  fractureTotalNumberDOFVelocityPressure + globalIndex, 1 ) ) );

			gmm::copy ( gmm::transposed(*Aup0), 
					    gmm::sub_matrix (*M_globalMatrix,
									    gmm::sub_interval (  fractureTotalNumberDOFVelocityPressure + fractureNumberCross + std::min(globalIndex, globalIndex2), 1 ),
									    gmm::sub_interval ( shiftIntersect [ id0 ], fractureNumberGlobalDOFVelocity [ id0 ] ) ) );

			gmm::copy ( *Aup1, 
					    gmm::sub_matrix (*M_globalMatrix,
									    gmm::sub_interval ( shiftIntersect [ id1 ], fractureNumberGlobalDOFVelocity [ id1 ] ),
									    gmm::sub_interval (  fractureTotalNumberDOFVelocityPressure + globalIndex2, 1 ) ) );

			gmm::copy ( gmm::transposed(*Aup1), 
					    gmm::sub_matrix (*M_globalMatrix,
									    gmm::sub_interval (  fractureTotalNumberDOFVelocityPressure + fractureNumberCross + std::min(globalIndex, globalIndex2), 1 ),
									    gmm::sub_interval ( shiftIntersect [ id1 ], fractureNumberGlobalDOFVelocity [ id1 ] ) ) );
						
         }
 					   	  
    }
    
    // Aggiorno prima le matrici di tutte le fratture che si intersecano formando un " IntBifurcation2 "
	for ( size_type i = 0; i < IntBifurcation2.size(); i++ )
	{
		sparseMatrixPtr_Type Aup0;
		
		// frattura lunga
		FractureHandlerPtr_Type f0, f1; 	
			
		if( IntBifurcation2 [ i ].getFracture (0)->getDOFBifurcation().size() != 0 )
		{
			f0 = IntBifurcation2 [ i ].getFracture (0);
			f1 = IntBifurcation2 [ i ].getFracture (1);
		}
		else
		{
			f0 = IntBifurcation2 [ i ].getFracture (1);
			f1 = IntBifurcation2 [ i ].getFracture (0);
		}

		
		size_type id0 = f0->getId();
		size_type id1 = f1->getId();
		
		Aup0.reset ( new sparseMatrix_Type ( fractureNumberGlobalDOFVelocity [ id0 ], 1 ) );
		gmm::clear(*Aup0);

		// aggiorno per la frattura 0
		getfem::darcy_A11F_Cross ( A11F [ id0 ],f0,
								   f0->getEtaTangentialInterpolated(), 
								   f1,
								   FractureHandler::FRACTURE_INTERSECT * ( id0 + 1 ) + id1 + 1);

		getfem::darcy_A12F_Cross ( A12F [ id0 ], f0, f1,
								   FractureHandler::FRACTURE_INTERSECT * ( id0 + 1 ) + id1 + 1 );

		
		const pairSizeVectorContainer_Type& intersectElementsGlobalIndex0 = f0->getFractureIntersectElementsGlobalIndex ();
		
		const size_type numIntersections = intersectElementsGlobalIndex0 [id1].size();
		
		const sizeVectorContainer_Type& intersectElements0 = f0->getFractureIntersectElements ();
		
		
		for ( size_type k = 0; k < numIntersections; ++k )
		{
			gmm::clear (*Aup0);

			getfem::velocityJump_Cross ( Aup0, f0, f1, intersectElements0 [ id1 ][ k ] );
			
			gmm::copy ( *Aup0, 
						gmm::sub_matrix (*M_globalMatrix,
										gmm::sub_interval ( shiftIntersect [ id0 ], fractureNumberGlobalDOFVelocity [ id0 ] ),
										gmm::sub_interval (  fractureTotalNumberDOFVelocityPressure , 1 ) ) );

			gmm::copy ( gmm::transposed(*Aup0), 
						gmm::sub_matrix (*M_globalMatrix,
										gmm::sub_interval (  fractureTotalNumberDOFVelocityPressure, 1 ),
										gmm::sub_interval ( shiftIntersect [ id0 ], fractureNumberGlobalDOFVelocity [ id0 ] ) ) );

		}
						  
	}
        
        ///////////////////////////////////////////////////////////////////////
    
	/*
	* 		Copy blocks into the system matrix M_darcyGlobalMatrix
	* 		
	* 				[  A11  A12  0          ]
	* 		    M = [ -A12  A22  0          ]
	* 		    	[  0    0    A11F  A12F ]
	* 		        [  0    0   -A12F  0    ]    	
	*/ 

	// Shift for the fracture
	size_type fractureShift = 0;
	for ( size_type f = 0; f < numberFractures; ++f )
	{
	   // Copy the matrix A11F in M_darcyGlobalMatrix in the correct position
	   gmm::copy( *( A11F [ f ] ), gmm::sub_matrix( *M_globalMatrix, gmm::sub_interval( fractureShift, fractureNumberGlobalDOFVelocity [ f ] ),
	                   										     gmm::sub_interval(fractureShift, fractureNumberGlobalDOFVelocity [ f ] )));

	   // Copy the matrix A12F in M_darcyGlobalMatrix in the correct position
	   gmm::copy( *( A12F [ f ] ), gmm::sub_matrix(*M_globalMatrix, gmm::sub_interval( fractureShift, fractureNumberGlobalDOFVelocity [ f ]),
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
	
	sizeVector_Type shiftVelocity ( numberFractures );
	shiftVelocity [ 0 ] = fractureNumberDOFVelocity [ 0 ];
	                 
	for ( size_type f = 0; f < numberFractures; ++f )
	{
		if( f != 0)
		{
			shiftVelocity [ f ] = shiftVelocity [ f-1 ] + M_fractures->getFracture( f-1 )->getNumExtendedVelocity() + fractureNumberGlobalDOFPressure [ f-1 ] + fractureNumberGlobalDOFVelocity [ f ];	
		}
		  
	}

    // Aggiorno ora la matrice globale imponendo le condizioni di interfaccia per la biforcazione
	
    for ( size_type i = 0; i < IntBifurcation.size(); i++ )
    {
    	sparseMatrixPtr_Type Aup0, Aup1, Aup2, Aup3;
    	   	
    	FractureHandlerPtr_Type f0 = IntBifurcation [ i ].getFracture (0);
    	FractureHandlerPtr_Type f1 = IntBifurcation [ i ].getFracture (1);
    	FractureHandlerPtr_Type f2 = IntBifurcation [ i ].getFracture (2);
    	
    	size_type id0 = f0->getId();
    	size_type id1 = f1->getId();
    	size_type id2 = f2->getId();

        std::cout << "Bifurcation: " << std::endl;
		std::cout << " Fractures " << id0 << ", " << id1 << ", " << id2 << std::endl;
		
		const size_type globalIndex = GlobalIndex_Bifurcation( M_fractures , id0,id1, id2 );
		const size_type Index = fractureTotalNumberDOFVelocityPressure + globalIndex;
		
        Aup0.reset ( new sparseMatrix_Type ( 1, fractureTotalNumberDOFVelocityPressure + globalFractureNumber) );
        gmm::clear(*Aup0);

        Aup1.reset ( new sparseMatrix_Type ( 1, fractureTotalNumberDOFVelocityPressure + globalFractureNumber) );
        gmm::clear(*Aup1);

        Aup2.reset ( new sparseMatrix_Type ( 1, fractureTotalNumberDOFVelocityPressure + globalFractureNumber) );
        gmm::clear(*Aup2);
		
        Aup3.reset ( new sparseMatrix_Type ( 1, fractureTotalNumberDOFVelocityPressure + globalFractureNumber) );
	    gmm::clear(*Aup3);
       		

        MatrixBifurcationHandler_Type Matrix( dataFile );
		FracturePtrContainer_Type Fracture( 3 );
		Fracture[ 0 ] =f0;
		Fracture[ 1 ] =f1;
		Fracture[ 2 ] =f2;
		

		FracturePtrContainer_Type Fracture_copy( 3 );
		Fracture_copy = Fracture;
		
		Matrix.setMatrices( Fracture_copy );
		
		Matrix4d T = Matrix.T();
		
		scalar_type s = 0.;
		Matrix.computeScap ( s );
				
		scalarVector_Type DOF( 3 );
		scalarVector_Type DOF_v( 3 );
		
		DOF_v = getfem::setDOF_v( DOF, Fracture, Matrix );
		
		getfem::setAup_i( Aup0, 0 , id0, id1, id2, DOF, DOF_v, shiftIntersect, fractureNumberGlobalDOFVelocity,  T, Index );
		getfem::setAup_i( Aup1, 1 , id1, id0, id2, DOF, DOF_v, shiftIntersect, fractureNumberGlobalDOFVelocity,  T, Index );
		getfem::setAup_i( Aup2, 2 , id2, id0, id1, DOF, DOF_v, shiftIntersect, fractureNumberGlobalDOFVelocity,  T, Index );
		getfem::setAup_i( Aup3, 3 , id2, id0, id1, DOF, DOF_v, shiftIntersect, fractureNumberGlobalDOFVelocity,  T, Index, s );
        
		gmm::copy(*Aup0, gmm::sub_matrix(*M_globalMatrix, 
	    		gmm::sub_interval( shiftVelocity[ id0 ] + DOF[ 0 ], 1), 
	    		gmm::sub_interval( 0, fractureTotalNumberDOFVelocityPressure + globalFractureNumber ) ));
				
		gmm::copy(*Aup1, gmm::sub_matrix(*M_globalMatrix, 
	    		gmm::sub_interval( shiftVelocity[ id1 ] + DOF[ 1 ], 1), 
	    		gmm::sub_interval( 0, fractureTotalNumberDOFVelocityPressure + globalFractureNumber ) ));
		
		gmm::copy(*Aup2, gmm::sub_matrix(*M_globalMatrix, 
	    		gmm::sub_interval( shiftVelocity[ id2 ] + DOF[ 2 ], 1), 
	    		gmm::sub_interval( 0, fractureTotalNumberDOFVelocityPressure + globalFractureNumber ) ));
		
		gmm::copy(*Aup3, gmm::sub_matrix(*M_globalMatrix, 
	    		gmm::sub_interval( fractureTotalNumberDOFVelocityPressure  + 2*fractureNumberCross + i, 1), 
	    		gmm::sub_interval( 0, fractureTotalNumberDOFVelocityPressure + globalFractureNumber ) ));	
				
			
    }
    
    // Aggiorno ora la matrice globale imponendo le condizioni di interfaccia per la biforcazione di tipo 2
    for ( size_type i = 0; i < IntBifurcation2.size(); i++ )
    {
    	
		sparseMatrixPtr_Type Aup0, Aup1, Aup2, Aup3, Aup4;

    	FractureHandlerPtr_Type f0, f1;
		
		if( IntBifurcation2 [ i ].getFracture (0)->getDOFBifurcation().size() != 0 )
		{
			f0 = IntBifurcation2 [ i ].getFracture (0);
			f1 = IntBifurcation2 [ i ].getFracture (1);
		}
		else
		{
			f1 = IntBifurcation2 [ i ].getFracture (0);
			f0 = IntBifurcation2 [ i ].getFracture (1);
		}
			
    	size_type id0 = f0->getId();
    	size_type id1 = f1->getId();
		

        std::cout << "Bifurcation: " << std::endl;
		std::cout << " Fractures " << id0 << ", " << id1 << std::endl;


    	const pairSizeVectorContainer_Type& intersectElementsGlobalIndex0 = M_fractures->getFracture( id0 )->getFractureIntersectElementsGlobalIndex ();

		const size_type globalIndex01 =  intersectElementsGlobalIndex0[ id1 ][ 0 ].first;

        const pairSizeVectorContainer_Type& intersectElementsGlobalIndex1 = M_fractures->getFracture( id1 )->getFractureIntersectElementsGlobalIndex ();

		const size_type globalIndex10 =  intersectElementsGlobalIndex1[ id0 ][ 0 ].first;


		const size_type globalIndex = fmin ( globalIndex01, globalIndex10 );
		
		const size_type globalShift = fractureTotalNumberDOFVelocityPressure + globalIndex;

        Aup0.reset ( new sparseMatrix_Type ( 1, fractureTotalNumberDOFVelocityPressure + globalFractureNumber) );
        gmm::clear(*Aup0);

        Aup1.reset ( new sparseMatrix_Type ( 1, fractureTotalNumberDOFVelocityPressure + globalFractureNumber) );
        gmm::clear(*Aup1);

        Aup2.reset ( new sparseMatrix_Type ( 1, fractureTotalNumberDOFVelocityPressure + globalFractureNumber) );
        gmm::clear(*Aup2);

        Aup3.reset ( new sparseMatrix_Type ( 1, fractureTotalNumberDOFVelocityPressure + globalFractureNumber) );
	    gmm::clear(*Aup3);
		
        Aup4.reset ( new sparseMatrix_Type ( 1, fractureTotalNumberDOFVelocityPressure + globalFractureNumber) );
	    gmm::clear(*Aup4);


        MatrixBifurcationHandler_Type Matrix( dataFile, "Bifurcation2" );
		FracturePtrContainer_Type Fracture( 2 );
		Fracture[ 0 ] =f0;
		Fracture[ 1 ] =f1;


		FracturePtrContainer_Type Fracture_copy( 2 );
		Fracture_copy = Fracture;
		Matrix.setMatrices( Fracture_copy );

		
		Matrix4d T = Matrix.T();
		
		
		scalarVector_Type DOF_p0( 2 );
		scalarVector_Type DOF_v0( 4 );
		scalar_type DOF_p1;
		scalar_type DOF_v1;

		Matrix.SetDOFIntersecton( Fracture[ 0 ], DOF_p0[ 0 ] );
		Matrix.SetDOFIntersecton( Fracture[ 1 ], DOF_p1 );

		const size_type nbDOF =  Fracture[ 1 ]-> getMeshFEMPressure().nb_basic_dof();
			
		if( DOF_p1 == nbDOF - 1 )
		{
			DOF_v1 = DOF_p1 + 1.;
		}
		else
		{
			DOF_v1 = DOF_p1;   
		}
		
		DOF_p0[ 1 ] = shiftIntersect[ id0 ] + fractureNumberDOFVelocity [ 0 ] + fractureNumberDOFPressure  [ 0 ] + Fracture[ 0 ]->getNumExtendedVelocity();
		DOF_v0[ 0 ] = DOF_p0[ 0 ];
		DOF_v0[ 1 ] = DOF_p0[ 0 ]+1;
		DOF_v0[ 2 ] = shiftIntersect[ id0 ] + fractureNumberDOFVelocity [ 0 ];
		DOF_v0[ 3 ] = shiftIntersect[ id0 ] + fractureNumberDOFVelocity [ 0 ] + 1;
	
		
		(*Aup0) ( 0 , shiftIntersect[ id0 ] + DOF_v0[ 0 ] )  = 0.5;
		(*Aup0) ( 0 , shiftIntersect[ id0 ] + DOF_v0[ 3 ] )  = 0.5;
		(*Aup0) ( 0 , shiftIntersect[ id0 ] + DOF_p0[ 0 ] + fractureNumberGlobalDOFVelocity [ id0 ] ) = 1.*T( 0 , 0 );
		(*Aup0) ( 0 , shiftIntersect[ id1 ] + DOF_p1 + fractureNumberGlobalDOFVelocity [ id1 ] ) = 1.*T( 0 , 1 );
		(*Aup0) ( 0 , shiftIntersect[ id0 ] + DOF_p0[ 1 ] ) = 1.*T( 0 , 2 );
		(*Aup0) ( 0 , globalShift + 1 ) = 1.*T( 0 , 3 );
		(*Aup0) ( 0 , globalShift ) = -1.*( T( 0 , 0 ) + T( 0 , 1 ) + T( 0 , 2 ) + T( 0 , 3 ) );

		gmm::copy(*Aup0, gmm::sub_matrix(*M_globalMatrix,
	    		gmm::sub_interval( shiftVelocity[ id0 ]  + DOF_p0[ 0 ], 1),
	    		gmm::sub_interval( 0, fractureTotalNumberDOFVelocityPressure + globalFractureNumber ) ));

		
		(*Aup1) ( 0 , shiftIntersect[ id1 ] + DOF_v1 )  = 1.;
		(*Aup1) ( 0 , shiftIntersect[ id0 ] + DOF_p0[ 0 ] + fractureNumberGlobalDOFVelocity [ id0 ] ) = 1.*T( 1 , 0 );
		(*Aup1) ( 0 , shiftIntersect[ id1 ] + DOF_p1 + fractureNumberGlobalDOFVelocity [ id1 ] ) = 1.*T( 1 , 1 );
		(*Aup1) ( 0 , shiftIntersect[ id0 ] + DOF_p0[ 1 ] ) = 1.*T( 1 , 2 );
		(*Aup1) ( 0 , globalShift + 1 ) = 1.*T( 1 , 3 );
		(*Aup1) ( 0 , globalShift ) = -1.*( T( 1 , 0 ) + T( 1 , 1 ) + T( 1 , 2 ) + T( 1 , 3 ) );


		gmm::copy(*Aup1, gmm::sub_matrix(*M_globalMatrix,
	    		gmm::sub_interval( shiftVelocity[ id1 ] + DOF_p1, 1),
	    		gmm::sub_interval( 0, fractureTotalNumberDOFVelocityPressure + globalFractureNumber ) ));

		
		(*Aup2) ( 0 , shiftIntersect[ id0 ] + DOF_v0[ 1 ] ) = 0.5;
		(*Aup2) ( 0 , shiftIntersect[ id0 ] + DOF_v0[ 2 ] ) = 0.5;
		(*Aup2) ( 0 , shiftIntersect[ id0 ] + DOF_p0[ 0 ] + fractureNumberGlobalDOFVelocity [ id0 ] ) = 1.*T( 2 , 0 );
		(*Aup2) ( 0 , shiftIntersect[ id1 ] + DOF_p1 + fractureNumberGlobalDOFVelocity [ id1 ] ) = 1.*T( 2 , 1 );
		(*Aup2) ( 0 , shiftIntersect[ id0 ] + DOF_p0[ 1 ] ) = 1.*T( 2 , 2 );
		(*Aup2) ( 0 , globalShift + 1 ) = 1.*T( 2 , 3 );
		(*Aup2) ( 0 , globalShift ) = -1.*( T( 2 , 0 ) + T( 2 , 1 ) + T( 2 , 2 ) + T( 2 , 3 ) );

		gmm::copy(*Aup2, gmm::sub_matrix(*M_globalMatrix,
	    		gmm::sub_interval( DOF_p0[ 1 ], 1),
	    		gmm::sub_interval( 0, fractureTotalNumberDOFVelocityPressure + globalFractureNumber ) ));

		
		// Imponiamo condizione di non slip alla parete
		(*Aup3) ( 0 , shiftIntersect[ id0 ] + DOF_p0[ 0 ] + fractureNumberGlobalDOFVelocity [ id0 ] ) = 1.*T( 3 , 0 );
		(*Aup3) ( 0 , shiftIntersect[ id1 ] + DOF_p1 + fractureNumberGlobalDOFVelocity [ id1 ] ) = 1.*T( 3 , 1 );
		(*Aup3) ( 0 , shiftIntersect[ id0 ] + DOF_p0[ 1 ] ) = 1.*T( 3 , 2 );
		(*Aup3) ( 0 , globalShift + 1 ) = 1.*T( 3 , 3 );
		(*Aup3) ( 0 , globalShift ) = -1.*( T( 3 , 0 ) + T( 3 , 1 ) + T( 3 , 2 ) + T( 3 , 3 ) );

		gmm::copy(*Aup3, gmm::sub_matrix(*M_globalMatrix,
	    		gmm::sub_interval( globalShift + 1 , 1),
	    		gmm::sub_interval( 0, fractureTotalNumberDOFVelocityPressure + globalFractureNumber ) ));
				
		
		scalar_type s = 0.;

		Matrix.computeScap ( s );
		
		// velocità: attenzione alla convenzione dei segni!!
		
		(*Aup4) ( 0 , shiftIntersect[ id0 ] + DOF_v0[ 0 ] )  = -1./( 8.0 * s );
		(*Aup4) ( 0 , shiftIntersect[ id0 ] + DOF_v0[ 3 ] )  = 1./( 8.0 * s );
		(*Aup4) ( 0 , shiftIntersect[ id0 ] + DOF_v0[ 1 ]  )  = 1./( 8.0 * s );
		(*Aup4) ( 0 , shiftIntersect[ id0 ] + DOF_v0[ 2 ]  )  = -1./( 8.0 * s );
		
		if ( DOF_v1 == 0 )
		{
			(*Aup4) ( 0 , shiftIntersect[ id1 ] + DOF_v1 )  = 1./( 4.0 * s );
		}
		else
		{
			(*Aup4) ( 0 , shiftIntersect[ id1 ] + DOF_v1 )  = -1./( 4.0 * s );
		}


		// pressione
		(*Aup4) ( 0 , shiftIntersect[ id0 ] + DOF_p0[ 0 ] + fractureNumberGlobalDOFVelocity [ id0 ] ) = -1./4.;
		(*Aup4) ( 0 , shiftIntersect[ id1 ] + DOF_p1 + fractureNumberGlobalDOFVelocity [ id1 ] ) = -1./4.;
		(*Aup4) ( 0 , shiftIntersect[ id0 ] + DOF_p0[ 1 ] ) = -1./4.;
		(*Aup4) ( 0 , globalShift + 1 ) = -1./4.;

		// pressione media
		(*Aup4) ( 0 , globalShift ) = 1.;
		
		
		gmm::copy(*Aup4, gmm::sub_matrix(*M_globalMatrix,
	    		gmm::sub_interval( globalShift, 1),
	    		gmm::sub_interval( 0, fractureTotalNumberDOFVelocityPressure + globalFractureNumber ) ));
		
    }

     gmm::copy(*App, gmm::sub_matrix(*M_globalMatrix, 
    		gmm::sub_interval( fractureTotalNumberDOFVelocityPressure, fractureNumberCross ), 
    		gmm::sub_interval( fractureTotalNumberDOFVelocityPressure, 2*fractureNumberCross ) ));
    
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
        const BCPtr_Type& fractureBC = M_bcHandler->getFractureBC( f );

        for ( size_type i = 0; i < fractureNumberBoundaryDOF [ f ]; i++ )
        {
            const base_node node = fractureBC->getMeshFEM().point_of_basic_dof( i );

            (*(PneumannF [ f ])) [ i ] = M_fractures->getFracture( f )->getData().pressureExact( node );
        }

        (*(PneumannF [ f ])) [ 0 ] *= -1; 
    }
	
    for ( size_type f = 0; f < numberFractures; ++f )
    {
    	std::cout << "Fracture " << f << std::endl;
    	
        // Computes the right hand side rhs which stores the boundary conditions for the fracture
        getfem::darcy_dataF ( BstressF [ f ], B_vF [ f ], M_bcHandler,
                              M_fractures->getFracture( f ), M_mediumData->getPenaltyVector(),
                              M_mediumData->getInvK(), PneumannF [ f ], VdirichletF [ f ] );
        
        std::cout<< std::endl;
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

        gmm::clear(*(B_pF [ f ]));

        getfem::assembling_Source_BoundaryF ( B_pF [ f ], divF [ f ], M_fractures->getFracture( f ), FractureHandler::FRACTURE_UNCUT * ( f + 1 ) );

    }
    
    // Aggiorno prima le matrici di tutte le fratture che si intersecano formando un " Cross "
    for ( size_type i = 0; i < IntCross.size(); i++ )
    {
    	FractureHandlerPtr_Type f0 = IntCross [ i ].getFracture (0);
    	FractureHandlerPtr_Type f1 = IntCross [ i ].getFracture (1);
    	
    	size_type id0 = f0->getId();
    	size_type id1 = f1->getId();

    	getfem::assembling_SourceF ( B_pF [ id0 ], divF [ id0 ],
    								f0, f1, FractureHandler::FRACTURE_INTERSECT * ( id0 + 1 ) + id1 + 1 );
    	                
    	getfem::assembling_SourceF ( B_pF [ id1 ], divF [ id1 ],
    								f1, f0, FractureHandler::FRACTURE_INTERSECT * ( id1 + 1 ) + id0 + 1 );
    	
    }

    
    // Aggiorno le matrici di tutte le fratture che si intersecano formando un " IntBifurcation2 "
    
    for ( size_type i = 0; i < IntBifurcation2.size(); i++ )
    {
    	FractureHandlerPtr_Type f0 = IntBifurcation2 [ i ].getFracture (0);
    	FractureHandlerPtr_Type f1 = IntBifurcation2 [ i ].getFracture (1);
    	
    	size_type id0 = f0->getId();
    	size_type id1 = f1->getId();

    	getfem::assembling_SourceF ( B_pF [ id0 ], divF [ id0 ],
    								f0, f1, FractureHandler::FRACTURE_INTERSECT * ( id0 + 1 ) + id1 + 1 );
    	                
    }
    
    
    fractureShift = 0;
    
    for ( size_type f = 0; f < numberFractures; ++f )
    {

        for ( size_type i = fractureNumberGlobalDOFVelocity [ f ]; i < fractureNumberDOFVelocityPressure [ f ]; ++i )
        {
                (*M_globalRightHandSide) [ fractureShift + i ] += (*(B_pF [ f ])) [ i - fractureNumberGlobalDOFVelocity [ f ] ];
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

    M_exporter->spy(M_globalMatrix, "./Matlab/matrice.mm");
    M_exporter->spy(M_globalRightHandSide, "./Matlab/rhs.mm");

    return;
    
}// assembly



// Solve the Darcy for the governing flux and do the time loop for the evolution problem
void DarcyFractured::solve ( )
{  
	// numero totale delle fratture
	const scalar_type numberFractures = M_fractures->getNumberFractures();
   
	// vettore dei gradi di libertà per la velocità per le fratture
    sizeVector_Type fractureNumberDOFVelocity(numberFractures);
    
    // numero di gradi di libertà totali per la velocità (con quelli estesi) per ogni frattura 
    sizeVector_Type fractureNumberGlobalDOFVelocity(numberFractures);

    // vettore dei gradi di libertà per la pressione per le fratture
    sizeVector_Type fractureNumberDOFPressure(numberFractures);
    
    // numero di gradi di libertà totali per la pressione (con quelli estesi) per ogni frattura 
    sizeVector_Type fractureNumberGlobalDOFPressure(numberFractures);

    // numero totale dei gradi di libertà ( velocità + pressione ) per ogni frattura
    sizeVector_Type fractureNumberDOFVelocityPressure(numberFractures);

    // numero complessivo dei gradi di libertà ( pressione + velocità ) 
    size_type fractureTotalNumberDOFVelocityPressure(0);

    
    for ( size_type f = 0; f < numberFractures; ++f )
    {
    	// numero di gradi di libertà per la velocità senza quelli estesi per la frattura f
        fractureNumberDOFVelocity [ f ] = M_fractures->getFracture( f )->getMeshFEMVelocity().nb_dof();

        // numero di gradi di libertà per la velocità con quelli estesi per la frattura f
        fractureNumberGlobalDOFVelocity [ f ] = fractureNumberDOFVelocity [ f ] + M_fractures->getFracture( f )->getNumExtendedVelocity();

        // numero di gradi di libertà per la pressione senza quelli estesi per la frattura f
        fractureNumberDOFPressure [ f ] = M_fractures->getFracture( f )->getMeshFEMPressure().nb_dof();

        // numero di gradi di libertà per la pressione con quelli estesi per la frattura f
        fractureNumberGlobalDOFPressure [ f ] = fractureNumberDOFPressure [ f ] + M_fractures->getFracture( f )->getNumExtendedPressure();

        // numero totale di gradi di libertà per la velocità e la pressione per la frattura f
        fractureNumberDOFVelocityPressure [ f ] = fractureNumberGlobalDOFVelocity [ f ] + fractureNumberGlobalDOFPressure [ f ];

        // numero totale di gradi di libertà per la pressione e la velocità
        fractureTotalNumberDOFVelocityPressure += fractureNumberDOFVelocityPressure [ f ];

        M_fracturePressure [ f ].reset(new scalarVector_Type( fractureNumberGlobalDOFPressure [ f ], 0));

        M_fractureVelocity [ f ].reset(new scalarVector_Type( fractureNumberGlobalDOFVelocity [ f ], 0));

    }
    
    // Solve the Darcy problem
    std::cout << std::endl << "Solving problem in Omega..." << std::flush;
    scalar_type roundConditionNumber;

    SuperLU_solve(*M_globalMatrix, *M_velocityAndPressure,
                  *M_globalRightHandSide, roundConditionNumber);
  
    size_type fractureShift = 0;
        
    for ( size_type f = 0; f < numberFractures; ++f )
    { 
        // Extract the dual in the fracture
        gmm::copy( gmm::sub_vector( *M_velocityAndPressure, gmm::sub_interval( fractureShift, fractureNumberGlobalDOFVelocity [ f ])),
                  *( M_fractureVelocity [ f ] ) );

        // Extract the primal in the fracture
        gmm::copy( gmm::sub_vector( *M_velocityAndPressure, 
        							gmm::sub_interval( fractureShift + fractureNumberGlobalDOFVelocity [ f ], fractureNumberGlobalDOFPressure [ f ])), 
        		   *(M_fracturePressure [ f ]));

        // Update the shift
        fractureShift += fractureNumberDOFVelocityPressure [ f ];
    }
    

    // the extra dof start at the end of the fracture dofs
  
    getfem::pfem fractureFETypePressure;
    getfem::pfem fractureFETypeVelocity;
    
    if ( numberFractures > 0 )
    {
        fractureFETypePressure = getfem::fem_descriptor( M_fractures->getFracture( 0 )->getData().getFEMTypePressure());
        fractureFETypeVelocity = getfem::fem_descriptor( M_fractures->getFracture( 0 )->getData().getFEMTypeVelocity());
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

            bgeot::base_node P1 = fracture->getMeshMapped().points( ) [i];


            scalar_type c= 1./fracture->getData().getSpatialDiscretization ();

            P1 [0] =  i*c;

            ordinataUncut [ i ] = fracture->getLevelSet()->getData()->y_map(P1);
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

        // per ogni frattura esporto la mesh reale
        osFileName << "cmeshF2-" << f << ".vtk";
        exportMesh(M_exporter->getFolder() + osFileName.str().c_str(), meshFcutMapped );

        getfem::mesh_fem mfproj ( meshFcutMapped, fracture->getMeshFEMPressure().get_qdim() );
        getfem::mesh_fem mfproj_v ( meshFcutMapped, fracture->getMeshFEMVelocity().get_qdim() );

        mfproj.set_classical_discontinuous_finite_element ( 0, 0.01 );
        mfproj_v.set_classical_discontinuous_finite_element ( 0, 0.01 );

        scalarVector_Type fracturePressureMeanUNCUT ( fractureNumberDOFPressure [ f ], 0. );
        scalarVector_Type fracturePressureInCut ( fractureNumberDOFPressure [ f ], 0. );
        scalarVector_Type fracturePressureOutCut ( fractureNumberDOFPressure [ f ], 0. );

        scalarVector_Type fractureVelocityMeanUNCUT ( fractureNumberDOFVelocity [ f ], 0. );
        scalarVector_Type fractureVelocityInCut ( fractureNumberDOFVelocity [ f ], 0. );
        scalarVector_Type fractureVelocityOutCut ( fractureNumberDOFVelocity [ f ], 0. );

        for ( size_type i = 0; i < fractureNumberDOFPressure [ f ]; ++i )
        {
            const size_type el =  fracture->getMeshFEMPressure().first_convex_of_basic_dof(i);
		
            if ( fracture->getMeshFlat().region ( FractureHandler::FRACTURE_UNCUT * ( f + 1 ) ).is_in(el) )
            {	
                fracturePressureMeanUNCUT [ i ] = (*(M_fracturePressure [ f ])) [ i ];
            }
         }

        for ( size_type i = 0; i < fractureNumberDOFVelocity [ f ]; ++i )
        {
            const size_type el =  fracture->getMeshFEMVelocity().first_convex_of_basic_dof(i);
		
            if ( fracture->getMeshFlat().region ( FractureHandler::FRACTURE_UNCUT * ( f + 1 ) ).is_in(el) )
            {	
                fractureVelocityMeanUNCUT [ i ] = (*(M_fractureVelocity [ f ])) [ i ];
            }
        }
        
        /*
         * getfem risolve il sistema per la mesh " piatta ", per avere i corretti valori di pressione devo interpolare
         * sulla mesh mappata i valori che ho ottenuto
         * 
         */
        scalarVector_Type fracturePressureMeanUNCUTInterpolated ( mfproj.nb_dof(), 0. );
        scalarVector_Type fracturePressureMeanInCutInterpolated ( mfproj.nb_dof(), 0. );
        scalarVector_Type fracturePressureMeanOutCutInterpolated ( mfproj.nb_dof(), 0. );

        scalarVector_Type fractureVelocityMeanUNCUTInterpolated ( mfproj_v.nb_dof(), 0. );
        scalarVector_Type fractureVelocityMeanInCutInterpolated ( mfproj_v.nb_dof(), 0. );
        scalarVector_Type fractureVelocityMeanOutCutInterpolated ( mfproj_v.nb_dof(), 0. );

        
        getfem::mesh_fem mfprojUncut ( fracture->getMeshMapped(), fracture->getMeshFEMPressure().get_qdim());
        mfprojUncut.set_finite_element ( fractureFETypePressure );

        getfem::mesh_fem mfprojUncut_v ( fracture->getMeshMapped(), fracture->getMeshFEMVelocity().get_qdim());
        mfprojUncut_v.set_finite_element ( fractureFETypeVelocity );

        getfem::interpolation ( mfprojUncut, mfproj, fracturePressureMeanUNCUT, fracturePressureMeanUNCUTInterpolated );
        
        getfem::interpolation ( mfprojUncut_v, mfproj_v, fractureVelocityMeanUNCUT, fractureVelocityMeanUNCUTInterpolated );

        const sizeVector_Type& extendedPressure = fracture->getExtendedPressure();
        
        const sizeVector_Type& extendedVelocity = fracture->getExtendedVelocity();

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
        
        for ( size_type otherFracture = 0; otherFracture < numberFractures; ++otherFracture )
        {
            if ( fracture->getLevelSetIntersect( otherFracture ).get() )
            {
                getfem::mesh_region regionMesh = fracture->getMeshFlat().region ( FractureHandler::FRACTURE_INTERSECT * ( f + 1 ) + otherFracture + 1 );
                size_type i_cv = 0;
                dal::bit_vector bc_cv = regionMesh.index();
                gmm::clear ( fractureVelocityInCut );
                gmm::clear ( fractureVelocityOutCut );

                for ( i_cv << bc_cv; i_cv != size_type(-1); i_cv << bc_cv )
                {
                    const size_type ibase = fracture->getMeshFEMVelocity().ind_basic_dof_of_element ( i_cv )[0];
                    const size_type position = size_type ( std::find ( extendedVelocity.begin(), extendedVelocity.end(), ibase )
                                                - extendedVelocity.begin());

                    const base_node point = mfprojUncut_v.point_of_basic_dof ( ibase );
                    const scalar_type levelSetValue = M_fractures->getFracture( otherFracture )->getLevelSet()->getData()->levelSetFunction( point );

                    if ( levelSetValue < 0 )
                    {
                        fractureVelocityInCut [ ibase ] = (*(M_fractureVelocity [ f ])) [ ibase ];
                        fractureVelocityOutCut [ ibase ] = (*(M_fractureVelocity [ f ])) [ fractureNumberDOFVelocity [ f ] + position ];
                    }
                    else
                    {
                        fractureVelocityOutCut [ ibase ] = (*(M_fractureVelocity [ f ])) [ ibase ];
                        fractureVelocityInCut [ ibase ] = (*(M_fractureVelocity [ f ])) [ fractureNumberDOFVelocity [ f ] + position ];
                    }

                }

                gmm::clear ( fractureVelocityMeanInCutInterpolated );

                getfem::interpolation ( mfprojUncut_v, mfproj_v,
                                        fractureVelocityInCut, fractureVelocityMeanInCutInterpolated );

                gmm::clear ( fractureVelocityMeanOutCutInterpolated );

                getfem::interpolation ( mfprojUncut_v, mfproj_v,
                                        fractureVelocityOutCut, fractureVelocityMeanOutCutInterpolated );

                i_cv = 0;
                bc_cv = meshFcutMapped.convex_index();
                for ( i_cv << bc_cv; i_cv != size_type(-1); i_cv << bc_cv )
                {
                    getfem::mesh_fem::ind_dof_ct idofs = mfproj_v.ind_basic_dof_of_element ( i_cv );
                    for ( size_type k = 0; k < idofs.size(); ++k )
                    {
                        const base_node node = mfproj_v.point_of_basic_dof(idofs [ k ]);
                        const scalar_type levelSetValue = M_fractures->getFracture( otherFracture )->getLevelSet()->getData()->levelSetFunction( node );
                        if ( levelSetValue < 0 )
                        {
                            fractureVelocityMeanUNCUTInterpolated [ idofs[k] ] += fractureVelocityMeanInCutInterpolated [ idofs[k] ];
                        }
                        else
                        {
                            fractureVelocityMeanUNCUTInterpolated [ idofs[k] ] += fractureVelocityMeanOutCutInterpolated [ idofs[k] ];
                        }

                    }

                }

            }
        }
        
		std::cout<<std::endl;
		std::cout<<std::endl;
		std::cout << "Frattura " << f <<" Pressione " << fracturePressureMeanUNCUTInterpolated << std::endl;
		std::cout<<std::endl;
		
		if( M_fractures->getFracture( f )->getDofIntersection().size() != 0 )
		{
			std::cout << "Frattura " << f <<" Pressione nel punto di intersezione: " << fracturePressureMeanUNCUTInterpolated[ M_fractures->getFracture( f )->getDofIntersection()[ 0 ] ] << std::endl;
			std::cout << " DOF =  " << M_fractures->getFracture( f )->getDofIntersection()[ 0 ] <<std::endl;
			if( M_fractures->getFracture( f )->getDofIntersection().size()==2 )
			{
				std::cout << "Frattura " << f <<" Pressione nel punto di intersezione: " << fracturePressureMeanUNCUTInterpolated[ M_fractures->getFracture( f )->getDofIntersection()[ 1 ] ] << std::endl;
				std::cout << " DOF =  " << M_fractures->getFracture( f )->getDofIntersection()[ 1 ] <<std::endl;
			}
		}
	
        osFileName.str("");
        osFileName << "fracturePressure" << f << ".vtk";

        exportSolution ( M_exporter->getFolder() + osFileName.str(), "Pressure",
                     mfproj, fracturePressureMeanUNCUTInterpolated );
        
        osFileName.str("");
        osFileName << "fractureVelocity" << f << ".vtk";
        
        exportSolution ( M_exporter->getFolder() + osFileName.str(), "Velocity",
                     mfproj_v, fractureVelocityMeanUNCUTInterpolated );
        
                
    }
        
    return;   
   
}
