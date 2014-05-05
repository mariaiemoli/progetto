#include "../include/FractureIntersect.h"

FractureIntersect::FractureIntersect ()
{
    // fill the number of basis functions for each extended element
    M_basisFunctionOfType [ Parallel ] = 1; // TRUCCO!!!!!
    M_basisFunctionOfType [ Cross ] = 1;

    regionLevelSetPair_Type coppia;
    // Fist the number of empty region
    coppia.second = 2;
    coppia.first = 1;
    M_subRegionIntersection [ coppia ] = Parallel;
    coppia.first = 0;
    M_subRegionIntersection [ coppia ] = Cross;

    // fill the list of coupling element sub region
    const size_type numMaxSubRegion = 4;
    M_subRegion.resize ( numMaxSubRegion );

    M_subRegion [ 0 ] = "++";
    M_subRegion [ 1 ] = "--";
    M_subRegion [ 2 ] = "+-";
    M_subRegion [ 3 ] = "-+";
}

void FractureIntersect::
constructIntesection ( getfem::mesh_level_set& meshLevelSet, const FracturePtrContainer_Type& fractures )
{
    const size_type numFractures = fractures.size();
    size_type globalIndexFractureIntersect = 0;

    // Find the intersections
    sizeVector_Type listOfConvex;
    std::vector<dal::bit_vector> listOfLevelSet_bitVector;

    meshLevelSet.find_level_set_potential_intersections ( listOfConvex, listOfLevelSet_bitVector );

    sizeVectorContainer_Type listOfLevelSet ( listOfConvex.size() );

    if ( listOfConvex.size() > 0 )
    {

        // We have to select which is the intersection type.
        for ( size_type i = 0; i < listOfConvex.size(); ++i )
        {
           stringContainer_Type regionActive ( M_subRegion );
           fromBitVectorToStdVector ( listOfLevelSet_bitVector [ i ], listOfLevelSet [ i ] );
           IntersectionType type = intersectionType ( meshLevelSet, regionActive,
                                                      listOfConvex [ i ], listOfLevelSet [ i ] );

           // Take the pointers of the fractures involved
           FracturePtrContainer_Type fracturesInvolved ( listOfLevelSet[i].size() );
           for ( size_type f = 0; f < fracturesInvolved.size(); ++f )
           {
                fracturesInvolved [ f ] = fractures [ listOfLevelSet [ i ] [ f ] ];
           }

           // Add a new element for the itersection
           IntersectData intersection;
           intersection.setIntersection ( listOfConvex[i], fracturesInvolved, regionActive, M_subRegion );
           M_intersections [ type ].push_back ( intersection );

        }

        const size_type crossNum = M_intersections [ Cross ].size();
        std::cout << "Num corss " << crossNum << std::endl;

        // We have to select which is the intersection type.
        for ( size_type i = 0; i < listOfConvex.size(); ++i )
        {
           stringContainer_Type regionActive ( M_subRegion );
           fromBitVectorToStdVector ( listOfLevelSet_bitVector [ i ], listOfLevelSet [ i ] );
           IntersectionType type = intersectionType ( meshLevelSet, regionActive,
                                                      listOfConvex [ i ], listOfLevelSet [ i ] );

           // Take the pointers of the fractures involved
           FracturePtrContainer_Type fracturesInvolved ( listOfLevelSet[i].size() );
           for ( size_type f = 0; f < fracturesInvolved.size(); ++f )
           {
                fracturesInvolved [ f ] = fractures [ listOfLevelSet [ i ] [ f ] ];
           }

           if ( type == Cross )
           {
//                for ( size_type f = 0; f < 2; ++f )
                {
                    std::cout << " globalIndexFractureIntersect " << globalIndexFractureIntersect << std::endl;
                    size_type indexTmp = globalIndexFractureIntersect;
                    fracturesInvolved [ 0 ]->setMeshLevelSetFracture ( *fracturesInvolved [ /*(f + 1)%2*/1 ], indexTmp );
                    indexTmp = globalIndexFractureIntersect + crossNum;
                    globalIndexFractureIntersect += fracturesInvolved [ 1 ]->setMeshLevelSetFracture ( *fracturesInvolved [ /*(f + 1)%2*/0 ], indexTmp );
                }
           }

        }

        for ( size_type f = 0; f < numFractures; ++f )
        {
            pairSizeVectorContainer_Type& fracture1 = fractures[f]->getFractureIntersectElementsGlobalIndex();

            for ( size_type otherFracture = 0; otherFracture < numFractures; ++otherFracture )
            {
                pairSizeVectorContainer_Type& fracture2 = fractures[ otherFracture ]->getFractureIntersectElementsGlobalIndex();

                const size_type numIntersections = fracture1 [ otherFracture ].size();
                for ( size_type numInt = 0; numInt < numIntersections; ++numInt )
                {
                    fracture2 [ f ][ numInt ].second = fracture1 [ otherFracture ][ numInt ].first;
                }
            }
        }

    }

    for ( size_type f = 0; f < numFractures; ++f )
    {
        pairSizeVectorContainer_Type& fracture1 = fractures[f]->getFractureIntersectElementsGlobalIndex();
        for ( size_type otherFracture = 0; otherFracture < numFractures; ++otherFracture )
        {
            const size_type numIntersections = fracture1 [ otherFracture ].size();
            for ( size_type numInt = 0; numInt < numIntersections; ++numInt )
            {
                std::cout << "Fratture " << f << " altra " << otherFracture << std::endl;
                std::cout << "NumInt " << numInt << " first " << fracture1 [ otherFracture ][ numInt ].first << " second " << fracture1 [ otherFracture ][ numInt ].second << std::endl;
            }
        }
    }

} // constructIntesection

void FractureIntersect::
findActiveRegion ( const scalarVector_Type& integrationValue, stringContainer_Type& regionActive) const
{
    size_type j = 0;
    for ( size_type i = 0; i < integrationValue.size(); ++i )
    {
        if ( integrationValue[i] == 0 )
        {
            regionActive.erase ( regionActive.begin() + j );
        }
        else
        {
            ++j;
        }
    }

} // findActiveRegion

FractureIntersect::IntersectionType FractureIntersect::
intersectionType ( getfem::mesh_level_set& meshLevelSet, stringContainer_Type& regionActive,
                   const size_type& elementID, const sizeVector_Type& levelSets )
{
    const size_type maxIntersection = regionActive.size();
    scalarVector_Type integrationValue ( maxIntersection );
    stringContainer_Type booleanOperation ( maxIntersection );

    // Fill the boolean operation between level sets. Change to take into account others types.
    for ( size_type i = 0; i < maxIntersection; ++i )
    {
        booleanOperation [i] = getOperation ( M_subRegion[i], levelSets );
    }

    for ( size_type i = 0; i < maxIntersection; ++i )
    {
        integrationValue [ i ] = integrateWithBooleanOperation ( meshLevelSet, elementID, booleanOperation [i] );
    }

    findActiveRegion ( integrationValue, regionActive );

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

size_type FractureIntersect::
getNumberIntersectionOfType ( IntersectionType type ) const
{
        return M_intersections.find( type )->second.size();
} // getNumberIntersectionOfType

size_type FractureIntersect::
getNumberIntersections () const
{
    mapIntersection_Type::const_iterator it;
    size_type numIntersect = 0;
    for ( it = M_intersections.begin(); it != M_intersections.end(); ++it )
    {
        numIntersect += getNumberIntersectionOfType ( it->first );
    }

    return numIntersect;
} // getNumberIntersections

size_type FractureIntersect::
getNumberType () const
{
    mapIntersection_Type::const_iterator it;
    size_type numType = 0;
    for ( it = M_intersections.begin(); it != M_intersections.end(); ++it )
    {
        ++numType;
    }

    return numType;
} // getNumberType

size_type FractureIntersect::
getBasisFunctionOfType ( IntersectionType type ) const
{
    return M_basisFunctionOfType.find ( type )->second;
} // getBasisFunctionOfType
