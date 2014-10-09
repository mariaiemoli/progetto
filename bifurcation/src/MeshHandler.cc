/** MeshHandler.cc
 *
 * Created on: Apr 1, 2011
 *
 * Author: fumagalli
 *
 */

#include "../include/MeshHandler.h"

MeshHandler::MeshHandler ( const GetPot& dataFile, const std::string& sectionDomain ):

			M_meshLevelSet ( M_mesh ),
			M_meshType(dataFile((sectionDomain + "meshTyper").data(),
					"GT_PK(2,1)")),
			M_spaceDimension(dataFile(
					(sectionDomain + "sapceDimension").data(), 2.)),
			M_fEMTypeVector(dataFile(
					(sectionDomain + "fEMTypeVelocity").data(), "FEM_RT0(2)")),
			M_meshFEMVector(M_mesh),
			M_fEMTypeScalar(dataFile(
					(sectionDomain + "fEMPressure").data(), "FEM_PK(2,0)")),
			M_meshFEMScalar(M_mesh),
			M_meshFEMCoefficients(M_mesh),
			M_integrationTypeVector(dataFile((sectionDomain
					+ "integrationTypeVelocity").data(), "IM_TRIANGLE(6)")),
			M_integrationTypeScalar(dataFile((sectionDomain
					+ "integrationTypePressure").data(), "IM_TRIANGLE(1)")),
			M_integrationMethodVector(M_mesh),
			M_integrationMethodScalar(M_mesh),
			M_meshExternal(dataFile((sectionDomain + "meshExternal").data(),
					"none")),
			M_meshFolder(dataFile((sectionDomain + "meshFolder").data(), "./")),
			M_spatialDiscretization(dataFile((sectionDomain
					+ "spatialDiscretization").data(), 10)),
			M_inclination(dataFile(
					(sectionDomain + "spatialInclination").data(), 0.)),
			M_lengthAbscissa(dataFile(
					(sectionDomain + "lengthAbscissa").data(), 1.)),
			M_lengthOrdinate(dataFile(
					(sectionDomain + "lengthOrdinate").data(), 1.)),
			M_lengthQuota(dataFile((sectionDomain + "lengthQuota").data(), 1.))
{
}// costruttore


void MeshHandler::setUpMesh ( )
{
	if ( M_meshExternal == "none" )
	{
		//------------------M_mediumMesh di Omega--------------------------------
		sizeVector_Type numberSubdivision(M_spaceDimension);
		std::fill(numberSubdivision.begin(), numberSubdivision.end(),
				M_spatialDiscretization);

        // Geometric transformation usign primal finite elements type
		M_geometricTransformation = bgeot::geometric_trans_descriptor( M_meshType);

		getfem::regular_unit_mesh(M_mesh, numberSubdivision,
				M_geometricTransformation);

		bgeot::base_matrix transformMatrix(M_spaceDimension, M_spaceDimension);

		scalarVector_Type length(3, 0);
		length [ 0 ] = M_lengthAbscissa;
		length [ 1 ] = M_lengthOrdinate;
		length [ 2 ] = M_lengthQuota;


        for ( size_type i = 0; i < M_spaceDimension; ++i )
        {
            transformMatrix(i, i) = (i < M_spaceDimension) ? length [ i ] : 1.0;
        }
        if ( M_spaceDimension > 1 )
        {
            transformMatrix(0, 1) = M_inclination * M_lengthOrdinate;
        }

        // scale the unit M_mediumMesh to [M_mediumLengthAbscissa,M_mediumLengthOrdinate]
        M_mesh.transformation(transformMatrix);

	}

	else
	      getfem::import_mesh((M_meshFolder + M_meshExternal).data(), M_mesh);

}// setUpMesh


void MeshHandler::setUpRegions ( const FracturesSetPtr_Type& fractures )
{
    // Select cut and uncut elements;
    size_type i_cv = 0;
    dal::bit_vector bv_cv = M_mesh.convex_index();

    // Add all the elements to the uncut region
    for ( i_cv << bv_cv; i_cv != size_type(-1); i_cv << bv_cv )
    {
        M_mesh.region(UNCUT_REGION).add(i_cv);
    }

    const size_type numberFractures = fractures->getNumberFractures ();

    if ( numberFractures > 0 )
    {

        M_extendedDOFScalar.resize ( numberFractures );
        M_extendedDOFVector.resize ( numberFractures );

        for ( size_type f = 0; f < numberFractures; ++f )
        {
            i_cv = 0;
            bv_cv = M_mesh.convex_index();

            for ( i_cv << bv_cv; i_cv != size_type(-1); i_cv << bv_cv )
            {
                if ( fractures->getFracture ( f )->getLevelSet()->getMesh().is_convex_cut(
                        i_cv) )
                {
		    bgeot::basic_mesh::ref_mesh_pt_ct nodes = M_mesh.points_of_convex ( i_cv );

		    if ( compreso( nodes, fractures->getFracture ( f ) ) )
		    {
	                    M_mesh.region(f + FractureData::FRACTURE).add(i_cv);
		    }
                }
            }

            // Check if the triangles are correctly cut, i.e. if some triangles have zero area in or out
            fixCutRegion ( fractures->getFracture ( f ) );
        }

        for ( size_type f = 0; f < numberFractures; ++f )
        {
            i_cv = 0;
            bv_cv = M_mesh.convex_index();
            for ( i_cv << bv_cv; i_cv != size_type(-1); i_cv << bv_cv )
            {
                if ( fractures->getFracture ( f )->getLevelSet()->getMesh().is_convex_cut(
                        i_cv) )
                {
		    bgeot::basic_mesh::ref_mesh_pt_ct nodes = M_mesh.points_of_convex ( i_cv );

		    if ( compreso( nodes, fractures->getFracture ( f ) ) )
		    {
	                    M_mesh.region(UNCUT_REGION).sup(i_cv);
		    }
                }
            }
        }

        // Create the extended dof for each fractures
        for ( size_type f = 0; f < numberFractures; ++f )
        {
            // Number of degree of freedom for the primal variable in the cut region
            dal::bit_vector cuttedRegionNumberDOFPressure =
                    M_meshFEMScalar.basic_dof_on_region(f
                            + FractureData::FRACTURE);

            // Fill the extended dof for the primal
            for ( dal::bv_visitor i(cuttedRegionNumberDOFPressure); !i.finished(); ++i )
            {
            	M_extendedDOFScalar [ f ].push_back(i);
            }

            // Number of degrees of freedom for the dual variable in the cutted region
            dal::bit_vector cuttedRegionNumberDOFVelocity =
                    M_meshFEMVector.dof_on_region(f + FractureData::FRACTURE);

            // Fill the extended dof for the dual
            for ( dal::bv_visitor i(cuttedRegionNumberDOFVelocity); !i.finished(); ++i )
            {
                M_extendedDOFVector [ f ].push_back(i);
            }
        }

        // Create the region mesh for the intersection, now we numerate respect to the element
        // number
        FractureIntersect::mapIntersection_Type::const_iterator itIntersect, itIntersectEnd;
        itIntersect = fractures->getIntersections()->getIntersections().begin();
        itIntersectEnd = fractures->getIntersections()->getIntersections().end();
        for ( ; itIntersect != itIntersectEnd; ++itIntersect )
        {
                const FractureIntersect::IntersectionType intersectionType = itIntersect->first;
                const IntersectDataContainer_Type& fractureIntersect = itIntersect->second;
                const size_type numIntersection = fractureIntersect.size();

                // loop on all the intersections
                for ( size_type k = 0; k < numIntersection; ++k )
                {
                        const IntersectData_Type& currentIntersection = fractureIntersect [ k ];
                        // erase the convex for the corresponding regions
                        const size_type convexId = currentIntersection.getElementID();
                        for ( size_type f = 0; f < currentIntersection.getNumFractures(); ++f )
                        {
                                const size_type fractureID = currentIntersection.getFracture ( f )->getId();
                                M_mesh.region ( FractureData::FRACTURE + fractureID ).sup ( convexId );
                        }

                        // add the convex to the new multiple cut region
                        M_mesh.region ( intersectionType + convexId ).add ( convexId );

                }
        }

        // Fill the extended dof for intersecting fractures
        FractureIntersect::mapIntersection_Type& mapIntersect = fractures->getIntersections()->
                                                                getIntersections();

        FractureIntersect::mapIntersection_Type::const_iterator begin = mapIntersect.begin();
        FractureIntersect::mapIntersection_Type::const_iterator end = mapIntersect.end();
        FractureIntersect::mapIntersection_Type::const_iterator it;

        size_type mapIDIntersectType = 0;
        for ( it = begin; it != end; ++it, ++mapIDIntersectType )
        {

                // resize to take into account element without extra dof
                const IntersectDataContainer_Type& intersection = it->second;
                const FractureIntersect::IntersectionType intersectionType = it->first;

                const size_type basisFunction = fractures->getIntersections()
                                                ->getBasisFunctionOfType ( intersectionType );

                for ( size_type j = 0; j < intersection.size(); ++j )
                {
                        // Current region
                        const size_type region = intersectionType + intersection[j].getElementID();

                        // Add the extended dof for the scalar
                        dal::bit_vector numberDOFScalar = M_meshFEMScalar.basic_dof_on_region ( region );

                        for ( dal::bv_visitor k ( numberDOFScalar ); !k.finished(); ++k )
                        {
                                // add in according to the extra dof needed
                                for ( size_type o = 0; o < basisFunction; ++o )
                                {
                                        M_extendedIntersectDOFScalar.push_back( k );
                                }
                        }

                        // Add the extended dof for the velocity
                        dal::bit_vector numberDOFVector = M_meshFEMVector.dof_on_region( region );

                        for ( dal::bv_visitor k ( numberDOFVector ); !k.finished(); ++k )
                        {
                                // Add in according to the extra dof needed
                                for ( size_type o = 0; o < basisFunction; ++o )
                                {
                                        M_extendedIntersectDOFVector.push_back( k );
                                }
                        }

                }

        }

        std::sort ( M_extendedIntersectDOFVector.begin(), M_extendedIntersectDOFVector.end() );
        sizeVector_Type::iterator itVect;
        itVect = unique ( M_extendedIntersectDOFVector.begin(), M_extendedIntersectDOFVector.end() );
        M_extendedIntersectDOFVector.resize ( itVect - M_extendedIntersectDOFVector.begin() );

    }

}


/*
controlare se queste funzioni non esistono già
attenzione che aglie estremi c'è il rischio che aggiunga due elementi tagliati di troppo, uno per ogni estremo
*/

bool MeshHandler::compreso ( const bgeot::basic_mesh::ref_mesh_pt_ct nodes, FractureHandlerPtr_Type& fracture )
{
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

    scalar_type translateAbscissa = fracture->getData().getTranslateAbscissa ();
    scalar_type lengthAbscissa = fracture->getData().getLengthAbscissa ();

    base_node nodo ( 1 );
    nodo [ 0 ] = 0;
    scalar_type translateOrdinata = fracture->getLevelSet()->getData()->y_map(nodo);
    nodo [ 0 ] = 1;
    scalar_type lengthOrdinata = fracture->getLevelSet()->getData()->y_map(nodo);

    if ( minx > translateAbscissa && maxx < translateAbscissa + lengthAbscissa)
    {
	if ( miny > translateOrdinata && maxy < lengthOrdinata )
        {
		return 1;
	}
	else
	{
		return 0;
	}

    }

    else if ( minx < translateAbscissa && maxx > translateAbscissa )
    {
	if ( miny < translateOrdinata && maxy < translateOrdinata )
        {
		return 0;
	}
	else
	{
		return 1;
	}
    }

    else if ( minx < translateAbscissa  + lengthAbscissa && maxx > translateAbscissa + lengthAbscissa )
    {
	if ( miny > lengthOrdinata && maxy > lengthOrdinata )
        {
		return 0;
	}
	else
	{
		return 1;
	}

    }

    return 0;

}



// Definizione degli elementi finiti
void MeshHandler::setUpFEM ( )
{
	// Dual variable spaces
    getfem::pfem FETypeVelocity = getfem::fem_descriptor(M_fEMTypeVector);

    // Integration type for the dual variable
    getfem::pintegration_method integrationTypeVector =
            getfem::int_method_descriptor(M_integrationTypeVector);

    // Integration method for the dual variable
    M_integrationMethodVector.set_integration_method(M_mesh.convex_index(),
            integrationTypeVector);

    // Finite element space for the dual variable
    M_meshFEMVector.set_qdim(M_spaceDimension);
    M_meshFEMVector.set_finite_element(M_mesh.convex_index(), FETypeVelocity);

    // Finite element type for the primal variable
    getfem::pfem FETypePressure = getfem::fem_descriptor(M_fEMTypeScalar);

    // Integration type for the primal variable
    getfem::pintegration_method integrationTypePressure =
            getfem::int_method_descriptor(M_integrationTypeScalar);

    // Integration method for the primal variable
    M_integrationMethodScalar.set_integration_method(M_mesh.convex_index(),
            integrationTypePressure);

    //  Finite element space for the primal variable
    M_meshFEMScalar.set_finite_element(M_mesh.convex_index(), FETypePressure);

    // Coefficient: P0 FEM

    // Finite element space the coefficients, the same as the primal finite element space
    M_meshFEMCoefficients.set_finite_element(M_mesh.convex_index(),
            FETypePressure);

}// setUpFEM


size_type MeshHandler::getCountExtendedDOFScalar ( const scalar_type& id ) const
{
    size_type total = 0;

    if ( id < 0 )
    {
        return total;
    }

    for ( size_type f = 0; f <= id; ++f )
    {
        total += M_extendedDOFScalar [ f ].size();
    }

    return total;
}// getCountExtendedDOFScalar


size_type MeshHandler::getCountExtendedIntersectDOFScalar () const
{
    return M_extendedIntersectDOFScalar.size();
}// getCountExtendedIntersectDOFScalar


// getCountExtendedDOFVector
size_type MeshHandler::getCountExtendedIntersectDOFVector () const
{
    return M_extendedIntersectDOFVector.size();
}// getCountExtendedIntersectDOFVector


size_type MeshHandler::getCountExtendedDOFVector ( const scalar_type& id ) const
{
    size_type total = 0;

    if ( id < 0 )
    {
        return total;
    }

    for ( size_type f = 0; f <= id; ++f )
    {
        total += M_extendedDOFVector [ f ].size();
    }

    return total;
}// getCountExtendedDOFVector


/*
 questa sistema la cut region creando una lista dei triangoli che non sono veramente tagliati, in pratica quando A1 o A2 sono minori di una certa tolleranza
 */
void MeshHandler::fixCutRegion ( const FractureHandlerPtr_Type& fracture )
{
    const LevelSetHandlerPtr_Type levelSet = fracture->getLevelSet();
    const size_type fractureID = fracture->getId();
    const size_type shiftScalar = M_meshFEMScalar.nb_dof();
    const size_type shiftLevelSet =
            levelSet->getLevelSet().get_mesh_fem().nb_dof();

    scalarVector_Type Va(shiftScalar, 0.), VaIn(shiftScalar, 0.), VaOut(
            shiftScalar, 0.);
    scalarVector_Type uni(shiftScalar, 1);

    scalarVector_Type ls_sec(shiftLevelSet, 0.);

    for ( size_type i = 0; i < shiftLevelSet; ++i )
    {
        ls_sec [ i ] = levelSet->getLevelSet().values(1) [ i ];
    }

    getfem::generic_assembly aree, areeIn, areeOut;
    aree.set("w=data(#1);"
        "t=comp(Base(#1).Base(#1));"
        "V(#1)+=t(:,i).w(i);");

    areeIn.set("w=data(#1);"
        "t=comp(Base(#1).Base(#1));"
        "V(#1)+=t(:,i).w(i);");

    areeOut.set("w=data(#1);"
        "t=comp(Base(#1).Base(#1));"
        "V(#1)+=t(:,i).w(i);");

    aree.push_mi(M_integrationMethodScalar);
    areeIn.push_mi(levelSet->getIntegrationMethodInside());
    areeOut.push_mi(levelSet->getIntegrationMethodOutside());

    aree.push_mf(M_meshFEMScalar);
    areeIn.push_mf(M_meshFEMScalar);
    areeOut.push_mf(M_meshFEMScalar);
    aree.push_mf(M_meshFEMScalar);
    areeIn.push_mf(M_meshFEMScalar);
    areeOut.push_mf(M_meshFEMScalar);

    aree.push_data(uni);
    areeIn.push_data(uni);
    areeOut.push_data(uni);

    aree.push_vec(Va);
    areeIn.push_vec(VaIn);
    areeOut.push_vec(VaOut);

    aree.assembly(M_mesh.region(fractureID + FractureData::FRACTURE));
    areeIn.assembly(M_mesh.region(fractureID + FractureData::FRACTURE));
    areeOut.assembly(M_mesh.region(fractureID + FractureData::FRACTURE));

    for ( size_type i = 0; i < shiftScalar; ++i )
    {
        if ( Va [ i ] > 0 )
        {
            if ( VaIn [ i ] / Va [ i ] <= 0 )
            {
                M_nonCut.push_back(M_meshFEMScalar.first_convex_of_basic_dof(i));
            }

            if ( VaOut [ i ] / Va [ i ] <= 0 )
            {
                M_nonCut.push_back(M_meshFEMScalar.first_convex_of_basic_dof(i));
            }

            VaIn [ i ] = VaIn [ i ] / Va [ i ];
            VaOut [ i ] = VaOut [ i ] / Va [ i ];
        }
        else
        {
            VaIn [ i ] = 2;
            VaOut [ i ] = 2;
        }
    }

} // fix_cut_region