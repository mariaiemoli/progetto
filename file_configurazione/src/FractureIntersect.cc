/** FractureIntersect.cc
 *
 */

#include "../include/FractureIntersect.h"


FractureIntersect::FractureIntersect ()
{
    // fill the number of basis functions for each extended element
    M_basisFunctionOfType [ Parallel ] = 1; 
    M_basisFunctionOfType [ Cross ] = 1;
    M_basisFunctionOfType [ Bifurcation ] = 1;

    regionLevelSetPair_Type coppia;
    
    // Assegno dei flags ad ogni sottoregione per distinguere il tipo di intersezione
    coppia.second = 3;
    coppia.first = 0;
    M_subRegionIntersection [ coppia ] = Bifurcation;
    coppia.first = 1;
    M_subRegionIntersection [ coppia ] = Cross;
    coppia.first = 2;
    M_subRegionIntersection [ coppia ] = Parallel;
    coppia.second = 2;
    coppia.first = 0;
    M_subRegionIntersection [ coppia ] = Cross;
    coppia.first = 1;
    M_subRegionIntersection [ coppia ] = Parallel;

} // costruttore



void FractureIntersect::
constructIntesection ( const getfem::mesh& mesh, getfem::mesh_level_set& meshLevelSet, const FracturePtrContainer_Type& fractures )
{
    // Assegno una numerazione globale ai Cross e alle Biforcazioni
    size_type globalIndexBifurcation = 0;
    size_type globalIndexCross = 0;

    // Cerco le intersezioni
    sizeVector_Type listOfConvex;
    std::vector<dal::bit_vector> listOfLevelSet_bitVector;

    /*
	 * Andiamo a pulire il vettore dei DOF_Free in modo che solo nel caso ci siano biforcazioni venga toccato
	 */	
	for ( size_type f=0; f < fractures.size(); f++ )
	{
		fractures[ f ] -> clearDofFree( );
	}
	
	/*
     * La chiamata di funzione 
     * 				meshLevelSet.find_level_set_potential_intersections ( listOfConvex, listOfLevelSet_bitVector )
     * restituisce un vettore con gli indici degli elementi della mesh di supporto tagliati da più di una frattura e per ciascun elemento un
     * vettore con gli indici delle fratture che lo attraversano. 
     * Da questa lista di possibili intersezioni devo estrarre quelle che effettivamente sono intersezioni e di che tipo sono.
     */
    meshLevelSet.find_level_set_potential_intersections ( listOfConvex, listOfLevelSet_bitVector );

    sizeVectorContainer_Type listOfLevelSet ( listOfConvex.size() );

    if( listOfConvex.size() > 0 )
	{
    	for ( size_type i = 0 ; i < listOfConvex.size(); ++i )
		{
    		// Conversione da un vettore di bit a un vettore di interi
	        fromBitVectorToStdVector ( listOfLevelSet_bitVector [ i ], listOfLevelSet [ i ] );

	        // Per ogni elemento della mesh di supporto in cui passano almeno due fratture verifico il tipo di intersezione
	        IntersectionType type = intersectionType ( meshLevelSet, listOfConvex [ i ], listOfLevelSet [ i ] );
	        
			// Prendo i puntatori alle fratture coinvolte
			FracturePtrContainer_Type fracturesInvolved ( listOfLevelSet[i].size() );

			for ( size_type f = 0; f < fracturesInvolved.size(); ++f )
			{
				 fracturesInvolved [ f ] = fractures [ listOfLevelSet [ i ] [ f ] ];

			}

			// Costruisco la classe IntersectData per la nuova intersezione e la aggiungo in base al tipo
			IntersectData intersection;

			intersection.setIntersection ( listOfConvex[i], fracturesInvolved );
			
			
			M_intersections [ type ].push_back ( intersection );


		}

        const size_type crossNum = M_intersections [ Cross ].size();
        std::cout << "Num cross " << crossNum << std::endl;

        const size_type bifuNum = M_intersections [ Bifurcation ].size();
        std::cout << "Num bifurcation " << bifuNum << std::endl;


        const size_type parNum = M_intersections [ Parallel ].size();
        std::cout << "Num parallel " << parNum << std::endl;

        std::map < sizeVector_Type, IntersectionType> intersezioni;

        // Per ogni elemento della mesh di supporto in cui passano almeno due fratture verifico il tipo di intersezione
        for ( size_type i = 0; i < listOfConvex.size(); ++i )
        {
           fromBitVectorToStdVector ( listOfLevelSet_bitVector [ i ], listOfLevelSet [ i ] );
           IntersectionType type = intersectionType ( meshLevelSet, listOfConvex [ i ], listOfLevelSet [ i ] );

           // Prendo i puntatori alle fratture coinvolte
           FracturePtrContainer_Type fracturesInvolved ( listOfLevelSet[i].size() );
           
           for ( size_type f = 0; f < fracturesInvolved.size(); ++f )
           {
                fracturesInvolved [ f ] = fractures [ listOfLevelSet [ i ] [ f ] ];
                
           }

           // Intersezione di tipo Cross
           if ( type == Cross )
           {
			    if ( fracturesInvolved.size() == 2 )
                {
        		    // Do una numerazione globale alle intersezioni di tipo Cross
                    std::cout << " globalIndexCross " << globalIndexCross << std::endl;

                    /*
                     * indexTmp è un indice che mi serve per tenere traccia degli indici dei gradi di libertà aggiuntivi associati ad ogni frattura
                     * l'idea è quella di considerare prima tutte le intersezioni di tipo Cross e poi quelle di tipo Bifurcation 
                     */
                    size_type indexTmp = globalIndexCross;
                    
                    std::string C ( "Cross" );
                    fracturesInvolved [ 0 ]->setMeshLevelSetFracture ( *fracturesInvolved [ 1 ], indexTmp, C );
                    
                    indexTmp = globalIndexCross + crossNum;
                    globalIndexCross += fracturesInvolved [ 1 ]->setMeshLevelSetFracture ( *fracturesInvolved [ 0 ], indexTmp, C );
                    std::cout << " fracturesInvolved [ 0 ]: " << fracturesInvolved[0]->getId() << std::endl;
                                        std::cout << " fracturesInvolved [ 1 ]: " << fracturesInvolved[1]->getId() << std::endl;
                }

        	   else if ( fracturesInvolved.size() == 3 )
        	   {
        		   std::cout << " caso cross da sistemare, tre fratture in un triangolo ma due sole si intersecano" << std::endl;
        	   }
        	   
           }
           
           // Intersezione di tipo Bifurcation
           if ( type == Bifurcation )
           {
			   
			   // Anche per le intersezioni di tipo Bifurcation assegno una numerazione globale
               std::cout << " globalIndexBifurcation " << globalIndexBifurcation << std::endl;
               
               std::string B ( "Bifurcation" );
               
               // Attenzione: nel caso della biforcazione ogni frattura interseca due altre fratture!
               size_type indexTmp = 2*crossNum + globalIndexBifurcation;
               fracturesInvolved [ 0 ]->setMeshLevelSetFracture ( *fracturesInvolved [ 1 ], indexTmp, B );
  
               indexTmp = 2*crossNum + globalIndexBifurcation + bifuNum;
               fracturesInvolved [ 0 ]->setMeshLevelSetFracture ( *fracturesInvolved [ 2 ], indexTmp, B );

               indexTmp = 2*crossNum + globalIndexBifurcation + 2*bifuNum;
               fracturesInvolved [ 1 ]->setMeshLevelSetFracture ( *fracturesInvolved [ 0 ], indexTmp, B );

               indexTmp = 2*crossNum + globalIndexBifurcation + 3*bifuNum ;
               fracturesInvolved [ 1 ]->setMeshLevelSetFracture ( *fracturesInvolved [ 2 ], indexTmp, B );

               indexTmp = 2*crossNum + globalIndexBifurcation + 4*bifuNum;
               fracturesInvolved [ 2 ]->setMeshLevelSetFracture ( *fracturesInvolved [ 0 ], indexTmp, B );

               indexTmp = 2*crossNum + globalIndexBifurcation + 5*bifuNum;
               fracturesInvolved [ 2 ]->setMeshLevelSetFracture ( *fracturesInvolved [ 1 ], indexTmp, B );

               globalIndexBifurcation++;
			   
			   //Andiamo a settare il vettore di DOF_FREE
			   for( size_type f=0; f < fracturesInvolved.size(); f++ )
			   {
				   setDofFree( mesh, fracturesInvolved[ f ], listOfConvex[ i ] );
			   }

           }

        	
		   // Stampe dei GlobalIndex
           if ( type == Cross )
           {
        	   pairSizeVectorContainer_Type& fracture0 = fracturesInvolved[0]->getFractureIntersectElementsGlobalIndex();
        	   pairSizeVectorContainer_Type& fracture1 = fracturesInvolved[1]->getFractureIntersectElementsGlobalIndex();
				
				fracture1 [ fracturesInvolved[0]->getId()][ 0 ].second = fracture0 [ fracturesInvolved[1]->getId()][ 0 ].first;
	
				fracture0 [ fracturesInvolved[1]->getId()][ 0 ].second = fracture1 [ fracturesInvolved[0]->getId()][ 0 ].first;
				
				/*
				 * funziona solo se due fratture si intersecano una volta sola
				 * per ora noi lavoriamo solo con rette quindi il problema non si pone
				 */
	
				std::cout << "Fratture: " << fracturesInvolved[0]->getId() << "		" << fracturesInvolved[1]->getId() << std::endl;
	
				std::cout << fracturesInvolved[0]->getId()
								<< " :	 first " << fracture0 [ fracturesInvolved[1]->getId() ][ 0 ].first
								<< " second " << fracture0 [ fracturesInvolved[1]->getId() ][ 0 ].second << std::endl;
	
				std::cout << fracturesInvolved[1]->getId()
								<< " :	 first " << fracture1 [ fracturesInvolved[0]->getId() ][ 0 ].first
								<< " second " << fracture1 [ fracturesInvolved[0]->getId() ][ 0 ].second << std::endl;
	
			}

			/*
			else if ( type == Bifurcation )
			{    
				
				pairSizeVectorContainer_Type& fracture0 = fracturesInvolved[0]->getFractureIntersectElementsGlobalIndex();
				pairSizeVectorContainer_Type& fracture1 = fracturesInvolved[1]->getFractureIntersectElementsGlobalIndex();
				pairSizeVectorContainer_Type& fracture2 = fracturesInvolved[2]->getFractureIntersectElementsGlobalIndex();
				
				
				fracture0 [ fracturesInvolved[2]->getId() ][ 0 ].second = fracture0 [ fracturesInvolved[1]->getId() ][ 0 ].first;
	 		    
				fracture1 [ fracturesInvolved[0]->getId() ][ 0 ].second = fracture0 [ fracturesInvolved[2]->getId() ][ 0 ].first;
	 		    
				fracture1 [ fracturesInvolved[2]->getId() ][ 0 ].second = fracture1 [ fracturesInvolved[0]->getId() ][ 0 ].first;
	 		    
				fracture2 [ fracturesInvolved[0]->getId() ][ 0 ].second = fracture1 [ fracturesInvolved[2]->getId() ][ 0 ].first;
	      	   
				fracture2 [ fracturesInvolved[1]->getId() ][ 0 ].second = fracture2 [ fracturesInvolved[0]->getId() ][ 0 ].first;
	 		    
				fracture0 [ fracturesInvolved[1]->getId() ][ 0 ].second = fracture2 [ fracturesInvolved[1]->getId() ][ 0 ].first;
	 
				std::cout << "Fratture: " << fracturesInvolved[0]->getId() << "		" << fracturesInvolved[1]->getId() << "		" << fracturesInvolved[2]->getId() << std::endl;

				std::cout << fracturesInvolved[0]->getId()
								<< " :	 first " << fracture0 [ fracturesInvolved[1]->getId() ][ 0 ].first
								<< " second " << fracture0 [ fracturesInvolved[1]->getId() ][ 0 ].second << std::endl;
	
				std::cout << fracturesInvolved[0]->getId()
								<< " :	 first " << fracture0 [ fracturesInvolved[2]->getId() ][ 0 ].first
								<< " second " << fracture0 [ fracturesInvolved[2]->getId() ][ 0 ].second << std::endl;
	
				std::cout << fracturesInvolved[1]->getId()
								<< " :	 first " << fracture1 [ fracturesInvolved[0]->getId() ][ 0 ].first
								<< " second " << fracture1 [ fracturesInvolved[0]->getId() ][ 0 ].second << std::endl;
	
	
				std::cout << fracturesInvolved[1]->getId()
								<< " :	 first " << fracture1 [ fracturesInvolved[2]->getId() ][ 0 ].first
								<< " second " << fracture1 [ fracturesInvolved[2]->getId() ][ 0 ].second << std::endl;
	
				std::cout << fracturesInvolved[2]->getId()
								<< " :	 first " << fracture2 [ fracturesInvolved[0]->getId() ][ 0 ].first
								<< " second " << fracture2 [ fracturesInvolved[0]->getId() ][ 0 ].second << std::endl;
	
				std::cout << fracturesInvolved[2]->getId()
								<< " :	 first " << fracture2 [ fracturesInvolved[1]->getId() ][ 0 ].first
								<< " second " << fracture2 [ fracturesInvolved[1]->getId() ][ 0 ].second << std::endl;
	
			}
			*/

        }
	}
    
    return;

} // constructIntesection



FractureIntersect::IntersectionType FractureIntersect::
intersectionType ( getfem::mesh_level_set& meshLevelSet, const size_type& elementID, const sizeVector_Type& levelSets )
{
	const size_type maxIntersection = 2*levelSets.size();
	scalarVector_Type integrationValue ( maxIntersection );
    stringContainer_Type booleanOperation ( maxIntersection );
	stringContainer_Type M_subRegion ( maxIntersection );

	if ( levelSets.size() == 2)
	{
		M_subRegion [ 0 ] = "++";
		M_subRegion [ 1 ] = "--";
		M_subRegion [ 2 ] = "+-";
		M_subRegion [ 3 ] = "-+";

	}

	else if ( levelSets.size() == 3)
	{
		M_subRegion [ 0 ] = "+++";
		M_subRegion [ 1 ] = "+-+";
		M_subRegion [ 2 ] = "+--";
		M_subRegion [ 3 ] = "---";
		M_subRegion [ 4 ] = "-+-";
		M_subRegion [ 5 ] = "-++";
	}

	// Fill the boolean operation between level sets.
	for ( size_type i = 0; i < maxIntersection; ++i )
	{
		booleanOperation [i] = getOperation ( M_subRegion[i], levelSets );

	}

	for ( size_type i = 0; i < maxIntersection; ++i )
	{
		integrationValue [ i ] = integrateWithBooleanOperation ( meshLevelSet, elementID, booleanOperation [i] );
	}

    regionLevelSetPair_Type coppia;
    coppia.first = (size_type) std::count ( integrationValue.begin(), integrationValue.end(), 0 );
    coppia.second = levelSets.size();

    return M_subRegionIntersection [ coppia ];


} // intersectionType



scalar_type FractureIntersect::
integrateWithBooleanOperation ( getfem::mesh_level_set& meshLevelSet, const size_type& elementID,
                                const std::string& operation ) const
{
    // construct the integration method based on the meshLevelSet
    getfem::pintegration_method intTypeIM = getfem::int_method_descriptor ( "IM_TRIANGLE(6)" );
    getfem::mesh_im_level_set meshImLevelSet ( meshLevelSet, getfem::mesh_im_level_set::INTEGRATE_OUTSIDE );
    meshImLevelSet.set_integration_method ( meshLevelSet.linked_mesh().convex_index(), intTypeIM );
    meshImLevelSet.set_simplex_im ( intTypeIM );

    // Set the operation bewteen the level sets
    meshImLevelSet.set_level_set_boolean_operations ( operation );

    meshImLevelSet.adapt();

    getfem::generic_assembly assem;

    assem.set("a=comp(Base(#1));""V(#1)+=a(:);");

    assem.push_mi ( meshImLevelSet );

    // Mesh FEM
    GFMeshFEM_Type meshFEM ( meshLevelSet.linked_mesh() );
    getfem::pfem FEType = getfem::fem_descriptor ( "FEM_PK(2,0)" );
    meshFEM.set_finite_element ( meshLevelSet.linked_mesh().convex_index(), FEType );

    assem.push_mf ( meshFEM );

    // Vector to store the values
    std::vector < scalar_type > V ( meshFEM.nb_dof() );
    gmm::clear ( V );

    assem.push_vec ( V );

    getfem::mesh_region meshElement;
    meshElement.add ( elementID );

    assem.assembly ( meshElement );

    // extract the correct value
    const size_type position = (meshFEM.ind_basic_dof_of_element( elementID ))[0];

    return V[position];

} // integrateWithBooleanOperation


IntersectDataContainer_Type FractureIntersect::getCrossIntersections () const
{
	IntersectDataContainer_Type tmp;
	size_type Num_Cross = getNumberCross();

	if( Num_Cross !=0 )
	{
		for ( size_type i = 0; i< M_intersections.find( Cross )->second.size(); i++ )
			tmp.push_back( M_intersections.find( Cross )->second [ i ] );
	}

    return tmp;

}


IntersectDataContainer_Type FractureIntersect::getBifurcationIntersections () const
{
	IntersectDataContainer_Type tmp;
	
	size_type Num_Bifu = getNumberBifurcation();

	if( Num_Bifu !=0 )
	{
		for ( size_type i = 0; i< M_intersections.find( Bifurcation )->second.size(); i++ )
			tmp.push_back( M_intersections.find( Bifurcation )->second [ i ] );
	}
	
    return tmp;

}


size_type FractureIntersect::getNumberIntersectionOfType ( IntersectionType type ) const
{
	return M_intersections.find( type )->second.size();
} // getNumberIntersectionOfType


size_type FractureIntersect::getNumberCross () const
{
    mapIntersection_Type::const_iterator it;
    size_type numIntersect = 0;
    for ( it = M_intersections.begin(); it != M_intersections.end(); ++it )
    {
	if( it->first == Cross )
        	numIntersect += getNumberIntersectionOfType ( it->first );
    }

    return numIntersect;
} // getNumberCross


size_type FractureIntersect::getNumberBifurcation () const
{
	mapIntersection_Type::const_iterator it;
	size_type numIntersect = 0;
	for ( it = M_intersections.begin(); it != M_intersections.end(); ++it )
	{
	if( it->first == Bifurcation )
			numIntersect += getNumberIntersectionOfType ( it->first );
	}

	return numIntersect;
} // getNumberBifurcation


size_type FractureIntersect::getNumberIntersections () const
{
    mapIntersection_Type::const_iterator it;
    size_type numIntersect = 0;
    for ( it = M_intersections.begin(); it != M_intersections.end(); ++it )
    {
        numIntersect += getNumberIntersectionOfType ( it->first );
    }

    return numIntersect;
} // getNumberIntersections



size_type FractureIntersect::getNumberType () const
{
    mapIntersection_Type::const_iterator it;
    size_type numType = 0;
    for ( it = M_intersections.begin(); it != M_intersections.end(); ++it )
    {
        ++numType;
    }

    return numType;
} // getNumberType



size_type FractureIntersect::getBasisFunctionOfType ( IntersectionType type ) const
{
    return M_basisFunctionOfType.find ( type )->second;
} // getBasisFunctionOfType


void FractureIntersect::setDofFree( const getfem::mesh& M_mesh, FractureHandlerPtr_Type& fracture , size_type i )
{
	bgeot::basic_mesh::ref_mesh_pt_ct nodes = M_mesh.points_of_convex ( i );
	
	if( FindDOF_Intersection( nodes, fracture ) == -1 )
	{
		fracture->pushDofFree( 0 );
	}
	
	else if( FindDOF_Intersection( nodes, fracture ) == 1 )
	{
		size_type NDof = fracture -> getMeshFEMPressure().nb_basic_dof(); 
		fracture->pushDofFree( NDof-1 );
	}	
	
	return;
	
}

size_type FractureIntersect::FindDOF_Intersection( const bgeot::basic_mesh::ref_mesh_pt_ct nodes, FractureHandlerPtr_Type& fracture )
{
	// x e y dei nodi geometrici che definiscono il convesso
	scalar_type x1 = nodes [ 0 ] [ 0 ];
    scalar_type x2 = nodes [ 1 ] [ 0 ];
    scalar_type x3 = nodes [ 2 ] [ 0 ];

    scalar_type y1 = nodes [ 0 ] [ 1 ];
    scalar_type y2 = nodes [ 1 ] [ 1 ];
    scalar_type y3 = nodes [ 2 ] [ 1 ];

    scalar_type minx= fmin(x1, fmin( x2, x3 ) );
    scalar_type maxx= fmax(x1, fmax( x2, x3 ) );

    scalar_type miny= fmin(y1, fmin( y2, y3 ) );
    scalar_type maxy= fmax(y1, fmax( y2, y3 ) );

    // x_massima e x_minima dei punti appartenenti alla frattura
    scalar_type translateAbscissa = fracture->getData().getTranslateAbscissa ();
    scalar_type lengthAbscissa = fracture->getData().getLengthAbscissa ();

    base_node nodo ( 1 );
    nodo [ 0 ] = 0;
    base_node nodo1 ( 1 );
    nodo1 [ 0 ] = 1;

    // y_massima e y_minima dei punti appartenenti alla frattura
    scalar_type translateOrdinata = fmin( fracture->getLevelSet()->getData()->y_map(nodo), fracture->getLevelSet()->getData()->y_map(nodo1) );
    scalar_type lengthOrdinata = fmax( fracture->getLevelSet()->getData()->y_map(nodo), fracture->getLevelSet()->getData()->y_map(nodo1) );

    // controllo che, nel caso in cui la frattura non tagli tutta la mesh di supporto, il convesso sia effettivamente tagliato dalla frattura
    if ( minx >= translateAbscissa && maxx <= translateAbscissa + lengthAbscissa)
    {   
		if ( miny >= translateOrdinata && maxy <= lengthOrdinata )
		{	
			return 0;
		}
		else if ( miny <= lengthOrdinata && maxy >= lengthOrdinata )
		{	
			return 0;
		}
		else if ( miny <= translateOrdinata && maxy >= translateOrdinata )
		{	
			return 0;
		}
		else
		{	
			return 3;
		}

    }

    else if ( minx <= translateAbscissa && maxx >= translateAbscissa )
    {
		if ( miny < translateOrdinata && maxy < translateOrdinata )
		{
			return 3;
		}
		else if ( miny > lengthOrdinata && maxy > lengthOrdinata )
		{
			return 3;
		}

		else
		{
			return -1;
		}
    }

    else if ( minx <= translateAbscissa  + lengthAbscissa && maxx >= translateAbscissa + lengthAbscissa )
    {
    	if ( miny >= lengthOrdinata && maxy >= lengthOrdinata )
		{	
			return 3;
		}
    	else if ( miny <= translateOrdinata && maxy <= translateOrdinata )
		{	
			return 3;
		}
		else
		{	
			return 1;
		}
    }
	
	return 3;

}// compreso