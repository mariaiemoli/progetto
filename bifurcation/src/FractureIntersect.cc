/** FractureIntersect.cc
 *
 */

#include "../include/FractureIntersect.h"


FractureIntersect::FractureIntersect ()
{
    // fill the number of basis functions for each extended element
    M_basisFunctionOfType [ Parallel ] = 1; // TRUCCO!!!!!
    M_basisFunctionOfType [ Cross ] = 1;
    M_basisFunctionOfType [ Bifurcation ] = 1;

    regionLevelSetPair_Type coppia;
    // Fist the number of empty region
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
constructIntesection ( getfem::mesh_level_set& meshLevelSet, const FracturePtrContainer_Type& fractures )
{
    const size_type numFractures = fractures.size();

    size_type globalIndexBifurcation = 0;
    size_type globalIndexCross = 0;

    // Find the intersections
    sizeVector_Type listOfConvex;
    std::vector<dal::bit_vector> listOfLevelSet_bitVector;

    meshLevelSet.find_level_set_potential_intersections ( listOfConvex, listOfLevelSet_bitVector );

    sizeVectorContainer_Type listOfLevelSet ( listOfConvex.size() );

	if( listOfConvex.size() > 0 )
	{
		for ( size_type i = 0; i < listOfConvex.size(); ++i )
		{
	        fromBitVectorToStdVector ( listOfLevelSet_bitVector [ i ], listOfLevelSet [ i ] );

	        IntersectionType type = intersectionType ( meshLevelSet, listOfConvex [ i ], listOfLevelSet [ i ] );
	        
			// Take the pointers of the fractures involved
			FracturePtrContainer_Type fracturesInvolved ( listOfLevelSet[i].size() );

			for ( size_type f = 0; f < fracturesInvolved.size(); ++f )
			{
				 fracturesInvolved [ f ] = fractures [ listOfLevelSet [ i ] [ f ] ];

			}

			// Add a new element for the itersection
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

        // We have to select which is the intersection type.
        for ( size_type i = 0; i < listOfConvex.size(); ++i )
        {
           fromBitVectorToStdVector ( listOfLevelSet_bitVector [ i ], listOfLevelSet [ i ] );
           IntersectionType type = intersectionType ( meshLevelSet, listOfConvex [ i ], listOfLevelSet [ i ] );

           // Take the pointers of the fractures involved
           FracturePtrContainer_Type fracturesInvolved ( listOfLevelSet[i].size() );
           
           for ( size_type f = 0; f < fracturesInvolved.size(); ++f )
           {
                fracturesInvolved [ f ] = fractures [ listOfLevelSet [ i ] [ f ] ];
           }

           if ( type == Cross )
           {
        	   if ( fracturesInvolved.size() == 2 )
                {
                    std::cout << " globalIndexCross " << globalIndexCross << std::endl;

                    size_type indexTmp = globalIndexCross;
                    fracturesInvolved [ 0 ]->setMeshLevelSetFracture ( *fracturesInvolved [ 1 ], indexTmp );

                    indexTmp = globalIndexCross + crossNum;
                    globalIndexCross += fracturesInvolved [ 1 ]->setMeshLevelSetFracture ( *fracturesInvolved [ 0 ], indexTmp );

                }

        	   else if ( fracturesInvolved.size() == 3 )
        	   {
        		   std::cout << " caso cross da sistemare, tre fratture in un triangolo ma due sole si intersecano" << std::endl;
        	   }
           }

           if ( type == Bifurcation )
           {
               std::cout << " globalIndexBifurcation " << globalIndexBifurcation << std::endl;
               size_type indexTmp = 2*crossNum + globalIndexBifurcation;
               fracturesInvolved [ 0 ]->setMeshLevelSetFracture ( *fracturesInvolved [ 1 ], indexTmp );
               indexTmp = 2*crossNum + globalIndexBifurcation + bifuNum;
               fracturesInvolved [ 0 ]->setMeshLevelSetFracture ( *fracturesInvolved [ 2 ], indexTmp );

               indexTmp = 2*crossNum + globalIndexBifurcation + 2*bifuNum;
               fracturesInvolved [ 1 ]->setMeshLevelSetFracture ( *fracturesInvolved [ 0 ], indexTmp );
               indexTmp = 2*crossNum + globalIndexBifurcation + 3*bifuNum ;
               fracturesInvolved [ 1 ]->setMeshLevelSetFracture ( *fracturesInvolved [ 2 ], indexTmp );

               indexTmp = 2*crossNum + globalIndexBifurcation + 4*bifuNum;
               fracturesInvolved [ 2 ]->setMeshLevelSetFracture ( *fracturesInvolved [ 0 ], indexTmp );
               indexTmp = 2*crossNum + globalIndexBifurcation + 5*bifuNum;
               fracturesInvolved [ 2 ]->setMeshLevelSetFracture ( *fracturesInvolved [ 1 ], indexTmp );

               globalIndexBifurcation++;

           }


        if ( type == Cross )
        {
            pairSizeVectorContainer_Type& fracture0 = fracturesInvolved[0]->getFractureIntersectElementsGlobalIndex();
            pairSizeVectorContainer_Type& fracture1 = fracturesInvolved[1]->getFractureIntersectElementsGlobalIndex();

            fracture0 [ fracturesInvolved[1]->getId()][ 0 ].second = fracture1 [ fracturesInvolved[0]->getId()][ 0 ].first;
            fracture1 [ fracturesInvolved[0]->getId()][ 0 ].second = fracture0 [ fracturesInvolved[1]->getId()][ 0 ].first;

            /*
             * funziona solo se due fratture si intersecano una volta sola
             */

            std::cout << "Fratture: " << fracturesInvolved[0]->getId() << "		" << fracturesInvolved[1]->getId() << std::endl;

            std::cout << fracturesInvolved[0]->getId()
            				<< " :	 first " << fracture0 [ fracturesInvolved[1]->getId() ][ 0 ].first
            				<< " second " << fracture0 [ fracturesInvolved[1]->getId() ][ 0 ].second << std::endl;

            std::cout << fracturesInvolved[1]->getId()
            				<< " :	 first " << fracture1 [ fracturesInvolved[0]->getId() ][ 0 ].first
            				<< " second " << fracture1 [ fracturesInvolved[0]->getId() ][ 0 ].second << std::endl;

        }

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

        }
	}

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

    // count the number of zeros
    const size_type numberZero = (size_type) std::count ( integrationValue.begin(), integrationValue.end(), 0 );

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
	
    mapIntersection_Type::const_iterator it;
    
    for ( it = M_intersections.begin(); it != M_intersections.end(); ++it )
    {
    	if ( it->first == Cross)
        {
    		tmp = it->second;
        }
    }

    return tmp;
}


IntersectDataContainer_Type FractureIntersect::getBifurcationIntersections () const
{
	IntersectDataContainer_Type tmp;
	
    mapIntersection_Type::const_iterator it;
    
    for ( it = M_intersections.begin(); it != M_intersections.end(); ++it )
    {
    	if ( it->first == Bifurcation)
        {
    		tmp = it->second;
        }
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
    	if ( it->first == Cross)
        {
    		numIntersect++;
        }
    }

    return numIntersect;

} // getNumberCross


size_type FractureIntersect::getNumberBifurcation () const
{
    mapIntersection_Type::const_iterator it;
    size_type numIntersect = 0;
    for ( it = M_intersections.begin(); it != M_intersections.end(); ++it )
    {
    	if ( it->first == Bifurcation)
        {
    		numIntersect++;
        }
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






